/*
NOTE: For regression checking, check against branch "main_with_indexing_fix"
(specifically, commit: e8f2d9bac4bbf0ad0619930dbdd54d089d9cb351).

This branch fixes a "read one past end of array allocation" bug for offsetarrayA/Ac/B/Bc.
*/
#include "arena.h"
#include <stdlib.h>
#include "cache.h"
#include "estimation.h"
#include "headers.h"
#include "json.hpp"
#include "simulator.h"
#include <chrono>
#include "statistics.h"
#include <cstdlib>
#include <fstream>
#include <filesystem>

struct Arena *global_persist;
struct Arena *global_temp;
struct simulator_state sim;
struct cache cache;

using json = nlohmann::json;

struct matrix {
    b16 transpose;
    b16 dense;

    Coord nrows;
    Coord ncols;
    Coord nzM;    // number of non-zero elements

    Coord *M_backing;
    Coord *Mc_backing;

    Coord **M, **Mc;
    Coord *offsetarrayM, *offsetarrayMc;
};

struct tile_cost_features {
    f64 sparse_cost;
    f64 dense_cost;
    f64 c_nnz_estimate;
    f64 c_density_estimate;
    Coord row_count;
    Coord j_count;
    Coord k_count;
    Coord ti, tj, tk;
};

struct policy_region_features {
    f64 sparse_cost;
    f64 dense_cost;
    f64 c_nnz_estimate;
    f64 c_density_estimate;
    Coord tile_count;
    Coord ti, tj, tk;
};

static enum workload_mode parse_workload_mode(const json &config)
{
    if (!config.contains("workloadMode"))
        return WORKLOAD_MODE_AUTO;

    const std::string mode = config["workloadMode"].get<std::string>();
    if (mode == "auto")
        return WORKLOAD_MODE_AUTO;
    if (mode == "sparse")
        return WORKLOAD_MODE_SPARSE;
    if (mode == "dense")
        return WORKLOAD_MODE_DENSE;

    std::cerr << "Unsupported workloadMode: " << mode << std::endl;
    std::exit(1);
}

static const char *print_workload_mode(enum workload_mode mode)
{
    switch (mode) {
    case WORKLOAD_MODE_AUTO:   return "auto";
    case WORKLOAD_MODE_SPARSE: return "sparse";
    case WORKLOAD_MODE_DENSE:  return "dense";
    case WORKLOAD_MODE_MIXED:  return "mixed";
    }
    return "unknown";
}

static inline Coord tile_linear_idx(Coord major_tile, Coord minor_tile, Coord minor_ntiles)
{
    return major_tile * minor_ntiles + minor_tile;
}

static inline Coord triple_tile_linear_idx(const struct config *cfg, Coord ti, Coord tj, Coord tk)
{
    return (ti * cfg->ttj + tj) * cfg->ttk + tk;
}

static struct matrix_characterization characterize_matrix(const struct matrix *x, Coord tile_rows, Coord tile_cols)
{
    struct matrix_characterization stats = {0};
    if (x->nrows == 0 || x->ncols == 0)
        return stats;

    const f64 total_slots = static_cast<f64>(x->nrows) * static_cast<f64>(x->ncols);
    stats.density = total_slots == 0.0 ? 0.0 : static_cast<f64>(x->nzM) / total_slots;

    f64 sum = 0.0;
    f64 sum_sq = 0.0;
    Coord empty_rows = 0;
    for (Coord row = 0; row < x->nrows; ++row) {
        Coord row_nnz = x->offsetarrayM[row + 1] - x->offsetarrayM[row];
        sum += row_nnz;
        sum_sq += static_cast<f64>(row_nnz) * row_nnz;
        empty_rows += row_nnz == 0;
    }

    stats.mean_row_nnz = sum / x->nrows;
    const f64 mean_sq = sum_sq / x->nrows;
    const f64 variance = max(0.0, mean_sq - stats.mean_row_nnz * stats.mean_row_nnz);
    stats.row_nnz_stddev = std::sqrt(variance);
    stats.empty_row_ratio = static_cast<f64>(empty_rows) / x->nrows;

    if (tile_rows == 0 || tile_cols == 0) {
        stats.tile_fill_ratio = stats.density;
    } else {
        const Coord tile_nrows = div_rup(x->nrows, tile_rows);
        const Coord tile_ncols = div_rup(x->ncols, tile_cols);
        std::vector<Coord> tile_nnz(tile_nrows * tile_ncols, 0);
        for (Coord row = 0; row < x->nrows; ++row) {
            const Coord tile_row = row / tile_rows;
            for (Coord idx = x->offsetarrayM[row]; idx < x->offsetarrayM[row + 1]; ++idx) {
                const Coord col = x->M[row][idx - x->offsetarrayM[row]];
                const Coord tile_col = col / tile_cols;
                ++tile_nnz[tile_row * tile_ncols + tile_col];
            }
        }

        f64 tile_fill_sum = 0.0;
        for (Coord tile_row = 0; tile_row < tile_nrows; ++tile_row) {
            const Coord rows_in_tile = min(tile_rows, x->nrows - tile_row * tile_rows);
            for (Coord tile_col = 0; tile_col < tile_ncols; ++tile_col) {
                const Coord cols_in_tile = min(tile_cols, x->ncols - tile_col * tile_cols);
                const f64 tile_slots = static_cast<f64>(rows_in_tile) * cols_in_tile;
                const Coord nnz = tile_nnz[tile_row * tile_ncols + tile_col];
                tile_fill_sum += tile_slots == 0.0 ? 0.0 : static_cast<f64>(nnz) / tile_slots;
            }
        }
        stats.tile_fill_ratio = tile_fill_sum / (tile_nrows * tile_ncols);
    }

    stats.is_dense_candidate =
        (stats.density >= 0.18) ||
        (stats.tile_fill_ratio >= 0.30 && stats.empty_row_ratio <= 0.10) ||
        (stats.mean_row_nnz >= tile_cols * 0.50 && tile_cols > 0);
    return stats;
}

static std::vector<struct tile_characterization> characterize_matrix_tiles(
    const struct matrix *x,
    Coord tile_rows,
    Coord tile_cols)
{
    if (tile_rows == 0 || tile_cols == 0)
        return {};

    const Coord tile_nrows = div_rup(x->nrows, tile_rows);
    const Coord tile_ncols = div_rup(x->ncols, tile_cols);
    std::vector<struct tile_characterization> tiles(tile_nrows * tile_ncols, {0});

    for (Coord row = 0; row < x->nrows; ++row) {
        const Coord tile_row = row / tile_rows;
        Coord prev_tile_col = (Coord)-1;
        for (Coord idx = x->offsetarrayM[row]; idx < x->offsetarrayM[row + 1]; ++idx) {
            const Coord col = x->M[row][idx - x->offsetarrayM[row]];
            const Coord tile_col = col / tile_cols;
            const Coord tile_idx = tile_linear_idx(tile_row, tile_col, tile_ncols);
            ++tiles[tile_idx].nnz;
            if (tile_col != prev_tile_col) {
                ++tiles[tile_idx].active_major_count;
                prev_tile_col = tile_col;
            }
        }
    }

    for (Coord tile_row = 0; tile_row < tile_nrows; ++tile_row) {
        const Coord rows_in_tile = min(tile_rows, x->nrows - tile_row * tile_rows);
        for (Coord tile_col = 0; tile_col < tile_ncols; ++tile_col) {
            const Coord cols_in_tile = min(tile_cols, x->ncols - tile_col * tile_cols);
            const Coord tile_idx = tile_linear_idx(tile_row, tile_col, tile_ncols);
            const f64 tile_slots = static_cast<f64>(rows_in_tile) * cols_in_tile;
            tiles[tile_idx].fill_ratio =
                tile_slots == 0.0 ? 0.0 : static_cast<f64>(tiles[tile_idx].nnz) / tile_slots;
        }
    }

    return tiles;
}

static inline f64 estimate_sparse_tile_cost(
    const struct tile_characterization &a_tile,
    const struct tile_characterization &b_tile,
    Coord row_count,
    Coord j_count,
    Coord k_count,
    f64 c_est,
    f64 c_density)
{
    const f64 avg_b_fiber =
        static_cast<f64>(b_tile.nnz) / max((Coord)1, b_tile.active_major_count);
    const f64 sparse_mac = static_cast<f64>(a_tile.nnz) * max(1.0, avg_b_fiber);
    const f64 metadata_penalty =
        (1.0 - a_tile.fill_ratio) * row_count +
        (1.0 - b_tile.fill_ratio) * k_count +
        (j_count * 0.25);
    const f64 scatter_penalty =
        c_est * (0.20 + 0.40 * c_density) +
        max(0.0, c_est - static_cast<f64>(a_tile.active_major_count)) * 0.08;
    return sparse_mac + (a_tile.nnz + b_tile.nnz) + scatter_penalty + metadata_penalty;
}

static inline f64 estimate_dense_tile_cost(
    const struct tile_characterization &a_tile,
    const struct tile_characterization &b_tile,
    Coord row_count,
    Coord j_count,
    Coord k_count,
    f64 c_est,
    f64 c_density)
{
    const f64 words_j = div_rup(j_count, (Coord)64);
    const f64 dense_kernel = static_cast<f64>(row_count) * k_count * words_j;
    const f64 bitset_build = a_tile.nnz + b_tile.nnz;
    const f64 bitset_probe =
        c_est * (0.04 + 0.06 * (1.0 - c_density));
    const f64 scan_penalty =
        (row_count + k_count) * words_j * 0.5 +
        (1.0 - 0.5 * (a_tile.fill_ratio + b_tile.fill_ratio)) * j_count;
    return dense_kernel + bitset_build + bitset_probe + scan_penalty;
}

static inline f64 estimate_c_tile_nnz(
    const struct tile_characterization &a_tile,
    const struct tile_characterization &b_tile,
    Coord row_count,
    Coord j_count,
    Coord k_count)
{
    if (row_count == 0 || j_count == 0 || k_count == 0 || a_tile.nnz == 0 || b_tile.nnz == 0)
        return 0.0;

    const f64 a_fill = static_cast<f64>(a_tile.nnz) / max(1.0, static_cast<f64>(row_count) * j_count);
    const f64 b_fill = static_cast<f64>(b_tile.nnz) / max(1.0, static_cast<f64>(j_count) * k_count);
    const f64 expected_overlap = a_fill * b_fill * j_count;
    const f64 row_activity =
        static_cast<f64>(a_tile.active_major_count) / max<Coord>(1, row_count);
    const f64 shared_activity =
        static_cast<f64>(b_tile.active_major_count) / max<Coord>(1, j_count);
    const f64 activation_prob =
        (1.0 - std::exp(-expected_overlap)) * (0.35 + 0.65 * row_activity) * (0.25 + 0.75 * shared_activity);
    const f64 dense_cap = static_cast<f64>(row_count) * k_count;
    const f64 sparse_mac_cap =
        min(static_cast<f64>(a_tile.nnz) * max(1.0, static_cast<f64>(b_tile.nnz) / max<Coord>(1, b_tile.active_major_count)),
            dense_cap);
    return min(dense_cap, min(dense_cap * min(1.0, activation_prob), sparse_mac_cap));
}

static inline f64 transition_penalty(
    enum workload_mode prev_mode,
    enum workload_mode next_mode,
    const struct tile_cost_features &prev_tile,
    const struct tile_cost_features &next_tile,
    f64 base_penalty)
{
    if (prev_mode == next_mode)
        return 0.0;

    const f64 words_j = div_rup(next_tile.j_count, (Coord)64);
    const f64 locality_delta =
        std::abs(prev_tile.c_density_estimate - next_tile.c_density_estimate) * (next_tile.row_count + next_tile.k_count) +
        std::abs(prev_tile.c_nnz_estimate - next_tile.c_nnz_estimate) * 0.02;
    const b16 same_i = prev_tile.ti == next_tile.ti;
    const b16 same_j = prev_tile.tj == next_tile.tj;
    const f64 reuse_disruption = (same_i ? next_tile.row_count * 0.20 : 0.0) + (same_j ? next_tile.j_count * 0.35 : 0.0);

    if (prev_mode == WORKLOAD_MODE_SPARSE && next_mode == WORKLOAD_MODE_DENSE) {
        const f64 dense_setup =
            (next_tile.row_count + next_tile.k_count) * words_j * 0.70 +
            next_tile.c_nnz_estimate * 0.10;
        return base_penalty * 0.55 + dense_setup + locality_delta + reuse_disruption;
    }

    const f64 dense_teardown =
        prev_tile.c_nnz_estimate * 0.12 +
        prev_tile.k_count * 0.45 +
        prev_tile.c_density_estimate * prev_tile.row_count * 0.25;
    return base_penalty * 0.45 + dense_teardown + locality_delta + reuse_disruption;
}

static std::vector<enum workload_mode> build_tile_mode_map(
    const struct config *cfg,
    const struct matrix_characterization &a_stats,
    const struct matrix_characterization &b_stats,
    const std::vector<struct tile_characterization> &a_tiles,
    const std::vector<struct tile_characterization> &b_tiles,
    Coord *policy_span_ti_out,
    Coord *policy_span_tj_out,
    Coord *policy_span_tk_out,
    Coord *sparse_region_count_out,
    Coord *dense_region_count_out,
    f64 *sparse_cost_out,
    f64 *dense_cost_out,
    f64 *mixed_cost_out,
    f64 *oracle_cost_out,
    f64 *policy_overhead_out)
{
    const Coord total_tiles = cfg->tti * cfg->ttj * cfg->ttk;
    std::vector<f64> sparse_cost(total_tiles, 0.0);
    std::vector<f64> dense_cost(total_tiles, 0.0);
    std::vector<struct tile_cost_features> tile_features(total_tiles);

    f64 force_sparse_cost = 0.0;
    f64 force_dense_cost = 0.0;
    f64 oracle_tile_cost = 0.0;
    for (Coord ti = 0; ti < cfg->tti; ++ti) {
        const Coord row_count = min(cfg->iii, cfg->I - ti * cfg->iii);
        for (Coord tj = 0; tj < cfg->ttj; ++tj) {
            const Coord j_count = min(cfg->jjj, cfg->J - tj * cfg->jjj);
            const struct tile_characterization &a_tile =
                a_tiles[tile_linear_idx(ti, tj, cfg->ttj)];
            for (Coord tk = 0; tk < cfg->ttk; ++tk) {
                const Coord k_count = min(cfg->kkk, cfg->K - tk * cfg->kkk);
                const struct tile_characterization &b_tile =
                    b_tiles[tile_linear_idx(tj, tk, cfg->ttk)];
                const Coord tile_idx = triple_tile_linear_idx(cfg, ti, tj, tk);
                const f64 c_est = estimate_c_tile_nnz(a_tile, b_tile, row_count, j_count, k_count);
                const f64 c_density =
                    c_est / max(1.0, static_cast<f64>(row_count) * k_count);
                sparse_cost[tile_idx] =
                    estimate_sparse_tile_cost(a_tile, b_tile, row_count, j_count, k_count, c_est, c_density);
                dense_cost[tile_idx] =
                    estimate_dense_tile_cost(a_tile, b_tile, row_count, j_count, k_count, c_est, c_density);
                oracle_tile_cost += min(sparse_cost[tile_idx], dense_cost[tile_idx]);
                tile_features[tile_idx] = {
                    .sparse_cost = sparse_cost[tile_idx],
                    .dense_cost = dense_cost[tile_idx],
                    .c_nnz_estimate = c_est,
                    .c_density_estimate = c_density,
                    .row_count = row_count,
                    .j_count = j_count,
                    .k_count = k_count,
                    .ti = ti,
                    .tj = tj,
                    .tk = tk,
                };
                force_sparse_cost += sparse_cost[tile_idx];
                force_dense_cost += dense_cost[tile_idx];
            }
        }
    }

    std::vector<struct policy_region_features> regions;
    std::vector<Coord> tile_region_ids(total_tiles, 0);
    const f64 avg_density = 0.5 * (a_stats.density + b_stats.density);
    const Coord max_region_tiles =
        avg_density < 0.05 ? 8 :
        (avg_density < 0.25 ? 4 : 2);
    const f64 preference_similarity_threshold = 0.12;
    const f64 c_density_similarity_threshold = 0.08;
    const f64 strong_preference_discontinuity = 0.28;
    const f64 strong_c_density_discontinuity = 0.14;
    const f64 dense_growth_margin = 0.10;
    const f64 sparse_growth_margin = 0.04;

    auto dense_preference = [](const struct tile_cost_features &f) {
        return (f.sparse_cost - f.dense_cost) / max(1.0, f.sparse_cost + f.dense_cost);
    };

    for (Coord tile_idx = 0; tile_idx < total_tiles; ++tile_idx) {
        const struct tile_cost_features &tile = tile_features[tile_idx];
        if (regions.empty()) {
            regions.push_back({
                .sparse_cost = tile.sparse_cost,
                .dense_cost = tile.dense_cost,
                .c_nnz_estimate = tile.c_nnz_estimate,
                .c_density_estimate = tile.c_density_estimate,
                .tile_count = 1,
                .ti = tile.ti,
                .tj = tile.tj,
                .tk = tile.tk,
            });
            tile_region_ids[tile_idx] = 0;
            continue;
        }

        struct policy_region_features &region = regions.back();
        const f64 region_preference =
            (region.sparse_cost - region.dense_cost) / max(1.0, region.sparse_cost + region.dense_cost);
        const f64 tile_preference = dense_preference(tile);
        const b16 same_mode_preference =
            (region_preference >= 0.0 && tile_preference >= 0.0) ||
            (region_preference < 0.0 && tile_preference < 0.0);
        const b16 same_i = region.ti == tile.ti;
        const f64 region_avg_c_density =
            region.c_density_estimate / max((Coord)1, region.tile_count);
        const b16 close_preference =
            std::abs(region_preference - tile_preference) <= preference_similarity_threshold;
        const b16 close_c_density =
            std::abs(region_avg_c_density - tile.c_density_estimate)
            <= c_density_similarity_threshold;
        const b16 strong_discontinuity =
            std::abs(region_preference - tile_preference) >= strong_preference_discontinuity ||
            std::abs(region_avg_c_density - tile.c_density_estimate) >= strong_c_density_discontinuity;
        const b16 dense_region = region_preference > 0.0;
        const b16 tile_dense = tile_preference > 0.0;
        const b16 dense_growth_allowed =
            dense_region && tile_dense &&
            region_preference >= dense_growth_margin &&
            tile_preference >= dense_growth_margin;
        const b16 sparse_growth_allowed =
            !dense_region && !tile_dense &&
            region_preference <= -sparse_growth_margin &&
            tile_preference <= -sparse_growth_margin;
        const b16 same_tj = region.tj == tile.tj;
        const b16 boundary_friendly =
            same_tj || (!dense_region && !tile_dense);

        if (same_i && same_mode_preference && !strong_discontinuity &&
            close_preference && close_c_density && boundary_friendly &&
            (dense_growth_allowed || sparse_growth_allowed) &&
            region.tile_count < max_region_tiles) {
            region.sparse_cost += tile.sparse_cost;
            region.dense_cost += tile.dense_cost;
            region.c_nnz_estimate += tile.c_nnz_estimate;
            region.c_density_estimate += tile.c_density_estimate;
            region.tile_count += 1;
            tile_region_ids[tile_idx] = regions.size() - 1;
        } else {
            regions.push_back({
                .sparse_cost = tile.sparse_cost,
                .dense_cost = tile.dense_cost,
                .c_nnz_estimate = tile.c_nnz_estimate,
                .c_density_estimate = tile.c_density_estimate,
                .tile_count = 1,
                .ti = tile.ti,
                .tj = tile.tj,
                .tk = tile.tk,
            });
            tile_region_ids[tile_idx] = regions.size() - 1;
        }
    }

    const Coord total_regions = regions.size();
    for (struct policy_region_features &region : regions) {
        if (region.tile_count) {
            region.c_density_estimate /= region.tile_count;
            const f64 reuse_savings =
                (region.tile_count - 1) *
                (region.c_nnz_estimate / max((Coord)1, region.tile_count) * 0.06 +
                 region.c_density_estimate * 18.0);
            region.dense_cost = max(region.dense_cost - reuse_savings, 0.0);
        }
    }

    std::vector<f64> dp_sparse(total_regions, 0.0);
    std::vector<f64> dp_dense(total_regions, 0.0);
    std::vector<u8> prev_sparse(total_regions, WORKLOAD_MODE_SPARSE);
    std::vector<u8> prev_dense(total_regions, WORKLOAD_MODE_DENSE);
    std::vector<enum workload_mode> tile_modes(total_tiles, WORKLOAD_MODE_SPARSE);

    const f64 switch_penalty =
        0.010 * ((force_sparse_cost + force_dense_cost) / max((Coord)1, total_regions));

    dp_sparse[0] = regions[0].sparse_cost;
    dp_dense[0] = regions[0].dense_cost;
    prev_sparse[0] = WORKLOAD_MODE_SPARSE;
    prev_dense[0] = WORKLOAD_MODE_DENSE;

    for (Coord region_idx = 1; region_idx < total_regions; ++region_idx) {
        const f64 dense_to_sparse_penalty = transition_penalty(
            WORKLOAD_MODE_DENSE,
            WORKLOAD_MODE_SPARSE,
            {
                .sparse_cost = regions[region_idx - 1].sparse_cost,
                .dense_cost = regions[region_idx - 1].dense_cost,
                .c_nnz_estimate = regions[region_idx - 1].c_nnz_estimate,
                .c_density_estimate = regions[region_idx - 1].c_density_estimate,
                .row_count = 0,
                .j_count = 0,
                .k_count = 0,
                .ti = regions[region_idx - 1].ti,
                .tj = regions[region_idx - 1].tj,
                .tk = regions[region_idx - 1].tk,
            },
            {
                .sparse_cost = regions[region_idx].sparse_cost,
                .dense_cost = regions[region_idx].dense_cost,
                .c_nnz_estimate = regions[region_idx].c_nnz_estimate,
                .c_density_estimate = regions[region_idx].c_density_estimate,
                .row_count = 0,
                .j_count = 0,
                .k_count = 0,
                .ti = regions[region_idx].ti,
                .tj = regions[region_idx].tj,
                .tk = regions[region_idx].tk,
            },
            switch_penalty);
        const f64 sparse_to_dense_penalty = transition_penalty(
            WORKLOAD_MODE_SPARSE,
            WORKLOAD_MODE_DENSE,
            {
                .sparse_cost = regions[region_idx - 1].sparse_cost,
                .dense_cost = regions[region_idx - 1].dense_cost,
                .c_nnz_estimate = regions[region_idx - 1].c_nnz_estimate,
                .c_density_estimate = regions[region_idx - 1].c_density_estimate,
                .row_count = 0,
                .j_count = 0,
                .k_count = 0,
                .ti = regions[region_idx - 1].ti,
                .tj = regions[region_idx - 1].tj,
                .tk = regions[region_idx - 1].tk,
            },
            {
                .sparse_cost = regions[region_idx].sparse_cost,
                .dense_cost = regions[region_idx].dense_cost,
                .c_nnz_estimate = regions[region_idx].c_nnz_estimate,
                .c_density_estimate = regions[region_idx].c_density_estimate,
                .row_count = 0,
                .j_count = 0,
                .k_count = 0,
                .ti = regions[region_idx].ti,
                .tj = regions[region_idx].tj,
                .tk = regions[region_idx].tk,
            },
            switch_penalty);
        const f64 stay_sparse = dp_sparse[region_idx - 1] + regions[region_idx].sparse_cost;
        const f64 switch_to_sparse =
            dp_dense[region_idx - 1] + dense_to_sparse_penalty + regions[region_idx].sparse_cost;
        if (stay_sparse <= switch_to_sparse) {
            dp_sparse[region_idx] = stay_sparse;
            prev_sparse[region_idx] = WORKLOAD_MODE_SPARSE;
        } else {
            dp_sparse[region_idx] = switch_to_sparse;
            prev_sparse[region_idx] = WORKLOAD_MODE_DENSE;
        }

        const f64 stay_dense = dp_dense[region_idx - 1] + regions[region_idx].dense_cost;
        const f64 switch_to_dense =
            dp_sparse[region_idx - 1] + sparse_to_dense_penalty + regions[region_idx].dense_cost;
        if (stay_dense <= switch_to_dense) {
            dp_dense[region_idx] = stay_dense;
            prev_dense[region_idx] = WORKLOAD_MODE_DENSE;
        } else {
            dp_dense[region_idx] = switch_to_dense;
            prev_dense[region_idx] = WORKLOAD_MODE_SPARSE;
        }
    }

    enum workload_mode mode =
        (dp_dense[total_regions - 1] < dp_sparse[total_regions - 1])
            ? WORKLOAD_MODE_DENSE
            : WORKLOAD_MODE_SPARSE;
    std::vector<enum workload_mode> region_modes(total_regions, WORKLOAD_MODE_SPARSE);
    for (Coord region_idx = total_regions; region_idx-- > 0;) {
        region_modes[region_idx] = mode;
        if (region_idx == 0)
            break;
        mode = (region_modes[region_idx] == WORKLOAD_MODE_SPARSE)
            ? (enum workload_mode)prev_sparse[region_idx]
            : (enum workload_mode)prev_dense[region_idx];
    }

    for (Coord tile_idx = 0; tile_idx < total_tiles; ++tile_idx) {
        tile_modes[tile_idx] = region_modes[tile_region_ids[tile_idx]];
    }

    for (enum workload_mode region_mode : region_modes) {
        *dense_region_count_out += region_mode == WORKLOAD_MODE_DENSE;
        *sparse_region_count_out += region_mode == WORKLOAD_MODE_SPARSE;
    }

    *policy_span_ti_out = 0;
    *policy_span_tj_out = 0;
    *policy_span_tk_out = 0;
    *sparse_cost_out = force_sparse_cost;
    *dense_cost_out = force_dense_cost;
    *mixed_cost_out = min(dp_sparse[total_regions - 1], dp_dense[total_regions - 1]);
    *oracle_cost_out = oracle_tile_cost;
    *policy_overhead_out = max(0.0, *mixed_cost_out - oracle_tile_cost);
    return tile_modes;
}

static struct workload_characterization characterize_workload(
    enum workload_mode requested_mode,
    const struct matrix *matA,
    const struct matrix *matB,
    const struct config *cfg)
{
    const auto t0 = std::chrono::steady_clock::now();
    struct workload_characterization result = {
        .requested_mode = requested_mode,
        .selected_mode = WORKLOAD_MODE_SPARSE,
    };

    result.A = characterize_matrix(matA, cfg->iii, cfg->jjj);
    result.B = characterize_matrix(matB, cfg->jjj, cfg->kkk);
    result.A_tiles = characterize_matrix_tiles(matA, cfg->iii, cfg->jjj);
    result.B_tiles = characterize_matrix_tiles(matB, cfg->jjj, cfg->kkk);
    result.tile_modes.resize(cfg->tti * cfg->ttj * cfg->ttk, WORKLOAD_MODE_SPARSE);

    if (requested_mode == WORKLOAD_MODE_SPARSE) {
        result.selected_mode = WORKLOAD_MODE_SPARSE;
        std::fill(result.tile_modes.begin(), result.tile_modes.end(), WORKLOAD_MODE_SPARSE);
        result.sparse_tile_count = result.tile_modes.size();
        result.decision_reason = "forced sparse mode";
        return result;
    }
    if (requested_mode == WORKLOAD_MODE_DENSE) {
        result.selected_mode = WORKLOAD_MODE_DENSE;
        std::fill(result.tile_modes.begin(), result.tile_modes.end(), WORKLOAD_MODE_DENSE);
        result.dense_tile_count = result.tile_modes.size();
        result.decision_reason = "forced dense mode";
        return result;
    }

    const f64 avg_density = 0.5 * (result.A.density + result.B.density);
    const f64 avg_tile_fill = 0.5 * (result.A.tile_fill_ratio + result.B.tile_fill_ratio);
    const f64 sparse_score =
        (1.0 - result.A.density) +
        (1.0 - result.B.density) +
        result.A.empty_row_ratio +
        result.B.empty_row_ratio +
        avg_density * 1.5 +
        (result.A.row_nnz_stddev / max(1.0, result.A.mean_row_nnz + 1.0)) +
        (result.B.row_nnz_stddev / max(1.0, result.B.mean_row_nnz + 1.0));
    const f64 dense_score =
        (1.0 - result.A.tile_fill_ratio) * 0.5 +
        (1.0 - result.B.tile_fill_ratio) * 0.5 +
        (1.0 - result.A.density) * 0.25 +
        (1.0 - result.B.density) * 0.25;

    const bool dense_candidate = result.A.is_dense_candidate && result.B.is_dense_candidate;
    const bool obvious_dense = avg_density >= 0.35 && avg_tile_fill >= 0.35;
        result.tile_modes = build_tile_mode_map(
        cfg,
        result.A,
        result.B,
        result.A_tiles,
        result.B_tiles,
        &result.policy_region_span_ti,
        &result.policy_region_span_tj,
        &result.policy_region_span_tk,
        &result.sparse_region_count,
        &result.dense_region_count,
        &result.estimated_sparse_cost,
        &result.estimated_dense_cost,
        &result.estimated_mixed_cost,
        &result.estimated_oracle_tile_cost,
        &result.estimated_policy_overhead_cost);
    for (enum workload_mode tile_mode : result.tile_modes) {
        result.dense_tile_count += tile_mode == WORKLOAD_MODE_DENSE;
        result.sparse_tile_count += tile_mode == WORKLOAD_MODE_SPARSE;
    }

    if (result.dense_tile_count == result.tile_modes.size()) {
        result.selected_mode = WORKLOAD_MODE_DENSE;
        result.decision_reason =
            obvious_dense
                ? "auto selected dense mode for clearly dense matrix pair"
                : "auto selected dense mode from per-tile cost map";
    } else if (result.sparse_tile_count == result.tile_modes.size()) {
        result.selected_mode = WORKLOAD_MODE_SPARSE;
        result.decision_reason =
            (!dense_candidate || !(dense_score * 1.15 < sparse_score))
                ? "auto kept sparse mode due to conservative sparse-score margin"
                : "auto cost map still preferred sparse mode for every tile";
    } else {
        result.selected_mode = WORKLOAD_MODE_MIXED;
        result.decision_reason = "auto selected a mixed sparse/dense tile map via offline DP";
    }
    const auto t1 = std::chrono::steady_clock::now();
    result.policy_build_time_us =
        std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    return result;
}

i32 cmp_coord(const void *a, const void *b)
{
    Coord x = *(const Coord *)a;
    Coord y = *(const Coord *)b;
    return (x > y) - (x < y);
}

std::string get_matrix_path(const std::string& matrix_name) {
    std::string relative_roots[] = {
        "./data/",
        "./largedata/",
        "./dense/",
        "./bfs/"
    };

    for (auto& root : relative_roots) {
        // test root/name.mtx
        {
          std::filesystem::path candidate{root + matrix_name + ".mtx"};
          if(std::filesystem::exists(candidate)) return candidate;
        }

        // test root/name/name.mtx
        {
          std::filesystem::path candidate{root + matrix_name + '/' + matrix_name + ".mtx"};
          if(std::filesystem::exists(candidate)) return candidate;
        }
    }
    std::cerr << "Error: " << matrix_name << " not found.\n";
    std::exit(1);
}

void parse_matrix(FILE *f, struct matrix *x)
{
    if (x->dense) {
        Coord longer_dim= max(x->ncols, x->nrows);
        x->M_backing    = (Coord *)arena_push(global_persist, longer_dim*sizeof(*x->M_backing), __alignof__(*x->M_backing), 0);
        x->Mc_backing   = NULL;
        x->M            = (Coord **)arena_push(global_persist, (x->nrows+1)*sizeof(*x->M), __alignof__(*x->M), 0);
        x->Mc           = (Coord **)arena_push(global_persist, (x->ncols+1)*sizeof(*x->Mc), __alignof__(*x->Mc), 0);
        x->offsetarrayM = (Coord *)arena_push(global_persist, (x->nrows+1)*sizeof(x->offsetarrayM[0]), __alignof__(x->offsetarrayM[0]), 0);
        x->offsetarrayMc= (Coord *)arena_push(global_persist, (x->ncols+1)*sizeof(x->offsetarrayMc[0]), __alignof__(x->offsetarrayMc[0]), 0);

        // trick: specify first row/col only (since they are all identical) and let all csr/csc entries alias to the same row/col
        for (Coord i = 0; i < longer_dim; ++i)
            x->M_backing[i] = i;

        for (Coord i = 0; i < x->nrows+1; ++i) {
            x->M[i]             = x->M_backing;
            x->offsetarrayM[i]  = i*x->ncols;
        }
        for (Coord i = 0; i < x->ncols+1; ++i) {
            x->Mc[i]            = x->M_backing;
            x->offsetarrayMc[i] = i*x->nrows;
        }

        return;
    }
    struct Arena_Mark mark = arena_snap(global_temp);

    b16 transpose = x->transpose;

    // Lines limited to 1024 characters by spec
    #define BUFFER_NBYTES (u64)1024
    char readbuffer[BUFFER_NBYTES];

    Coord *raw_rows = (Coord *)arena_push(global_temp, x->nzM*sizeof(*raw_rows), __alignof__(*raw_rows), 0);
    Coord *raw_cols = (Coord *)arena_push(global_temp, x->nzM*sizeof(*raw_cols), __alignof__(*raw_cols), 0);
    Coord *row_lens = (Coord *)arena_push(global_temp, x->nrows*sizeof(*row_lens), __alignof__(*row_lens), 1); // must be zeroed
    Coord *col_lens = (Coord *)arena_push(global_temp, x->ncols*sizeof(*col_lens), __alignof__(*col_lens), 1); // must be zeroed
    Coord *offsets;

    x->M_backing    = (Coord *)arena_push(global_persist, x->nzM*sizeof(*x->M_backing), __alignof__(*x->M_backing), 0);
    x->Mc_backing   = (Coord *)arena_push(global_persist, x->nzM*sizeof(*x->Mc_backing), __alignof__(*x->Mc_backing), 0);
    x->M            = (Coord **)arena_push(global_persist, (x->nrows+1)*sizeof(*x->M), __alignof__(*x->M), 0);
    x->Mc           = (Coord **)arena_push(global_persist, (x->ncols+1)*sizeof(*x->Mc), __alignof__(*x->Mc), 0);
    x->offsetarrayM = (Coord *)arena_push(global_persist, (x->nrows+1)*sizeof(x->offsetarrayM[0]), __alignof__(x->offsetarrayM[0]), 0);
    x->offsetarrayMc= (Coord *)arena_push(global_persist, (x->ncols+1)*sizeof(x->offsetarrayMc[0]), __alignof__(x->offsetarrayMc[0]), 0);

    for (Coord i = 0; i < x->nzM; ++i) {
        /* NOTE(ejs): Each mtx row is a matrix entry.
        xx, yy [zz, lala]
            - xx  = row index
            - yy  = col index
            - zz  = real component of value? (ignored)
            - lala= imag. component of value? (ignored)
        */
        assert(fgets(readbuffer, BUFFER_NBYTES, f));
        Coord xx, yy;
        // f64 zz, lala;
        // assert(sscanf(readbuffer, "%lu %lu %lf %lf", &xx, &yy, &zz, &lala) >= 2);
        static_assert(sizeof(Coord) == sizeof(u64) || sizeof(Coord) == sizeof(u32), "unsupported Coord size");
        switch (sizeof(Coord)) {
        case sizeof(u64): assert(sscanf(readbuffer, "%u %u", &xx, &yy) == 2);   break;
        case sizeof(u32): assert(sscanf(readbuffer, "%u %u", &xx, &yy) == 2);       break;
        }

        // WARNING(ejs): mtx indices are stored 1-based; it converts 0-based representation in code
        Coord row_index = xx-1;
        Coord col_index = yy-1;
        if (transpose)
            swap(row_index, col_index);
        raw_rows[i] = row_index;
        raw_cols[i] = col_index;
        ++row_lens[row_index];
        ++col_lens[col_index];
    }

    // csr
    // prefix sum
    offsets = x->offsetarrayM;
    offsets[0] = 0;
    for (Coord i = 0; i < x->nrows; ++i)
        offsets[i+1] = offsets[i] + row_lens[i];

    // scatter cols into row-major order (reuse row_lens as cursor)
    memset(row_lens, 0, x->nrows * sizeof(*row_lens));
    for (Coord i = 0; i < x->nzM; ++i) {
        Coord r = raw_rows[i];
        Coord pos = offsets[r] + row_lens[r]++;
        x->M_backing[pos] = raw_cols[i];
    }

    // sort within each row by column
    for (Coord r = 0; r < x->nrows; ++r) {
        Coord *base = x->M_backing + offsets[r];
        qsort(base, row_lens[r], sizeof(*base), cmp_coord);
    }
    
    // build csr pointers
    for (Coord i = 0; i < x->nrows; ++i)
        x->M[i] = x->M_backing + offsets[i];
    x->M[x->nrows] = x->M_backing + x->offsetarrayM[x->nrows];

    // csc
    offsets = x->offsetarrayMc;
    offsets[0] = 0;
    for (Coord i = 0; i < x->ncols; ++i)
        offsets[i+1] = offsets[i] + col_lens[i];

    memset(col_lens, 0, x->ncols * sizeof(*col_lens));
    for (Coord i = 0; i < x->nzM; ++i) {
        Coord c = raw_cols[i];
        Coord pos = offsets[c] + col_lens[c]++;
        x->Mc_backing[pos] = raw_rows[i];
    }
    for (Coord c = 0; c < x->ncols; ++c) {
        Coord *base = x->Mc_backing + offsets[c];
        qsort(base, col_lens[c], sizeof(*base), cmp_coord);
    }
    for (Coord i = 0; i < x->ncols; ++i)
        x->Mc[i] = x->Mc_backing + offsets[i];
    x->Mc[x->ncols] = x->Mc_backing + x->offsetarrayMc[x->ncols];

    arena_rewind(mark);
}

// struct matrix_tiling_config {
//     struct Arena *a;
//     Coord major_dim;
//     Coord minor_dim;
//     Coord minor_tile_dim;
//     Coord **map;    // compressed sparse <axis> table
//     Coord *offsets;
// };
// struct matrix_tiling_luts {
//     Coord   **begins;
//     Coord   **sizes;
// };

// struct matrix_tiling_luts
// create_tiling_luts(const struct matrix_tiling_config *cfg)
// {
//     Coord   major_dim, minor_dim, minor_tile_dim, minor_ntiles;
//     Coord   *offsets, **map;
//     major_dim       = cfg->major_dim;
//     minor_dim       = cfg->minor_dim;
//     minor_tile_dim  = cfg->minor_tile_dim;
//     minor_ntiles    = div_rup(minor_dim, minor_tile_dim);
//     offsets         = cfg->offsets;
//     map             = cfg->map;

//     Coord begins[major_dim][minor_ntiles+1] = (Coord **)arena_push(cfg->a, major_dim*(minor_ntiles+1)*sizeof(**begins), __alignof__(**begins), 0);
//     Coord sizes[major_dim][minor_ntiles]    = (Coord **)arena_push(cfg->a, major_dim*(minor_ntiles+1)*sizeof(**begins), __alignof__(**begins), 0);
//     // sizes           = (Coord **)arena_push(cfg->a, major_dim*minor_tile_dim*sizeof(**sizes), __alignof__(**sizes), 1);
//     printf("begins: %p\n", begins);
//     printf("sizes: %p\n", sizes);
//     fflush(stdout);

//     for (Coord ri = 0; ri < major_dim; ++ri) {
//         Coord nzero_elem_cnt= offsets[ri+1] - offsets[ri];
//         Coord wj_tile_base  = 0;
//         Coord wj            = 0;
//         begins[ri][0]       = 0;
//         for (Coord rj = 0; rj < nzero_elem_cnt; ++rj) {
//             while (map[ri][rj] >= (wj_tile_base + minor_tile_dim)) {
//                 wj_tile_base += minor_tile_dim;
//                 ++wj;
//                 begins[ri][wj] = rj;
//             }
//             ++sizes[ri][wj];
//         }

//         while (wj < minor_ntiles) {
//             ++wj;
//             begins[ri][wj] = nzero_elem_cnt;
//         }
//     }

//     return (struct matrix_tiling_luts) {
//         .begins = begins,
//         .sizes  = sizes
//     };
// }

void reset_cursor(struct cursor *c)
{
    c->first= 1;
    c->ti   = 0;
    c->tj   = 0;
    c->tk   = 0;
}

void release_global_persist() { arena_release(global_persist); }
void release_global_temp()    { arena_release(global_temp);    }
void release_cache_backing()  { arena_release(cache.backing);  }

int main(int argc, char *argv[])
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <matrix1 name> <matrix2 name> <config filepath>\n"
                  << "config_filepath is a fully qualified path to the .json config for the run (likely config/config.json).\n"
                  << std::endl;
        return 1;
    }

    global_persist  = arena_alloc(16*GB, MB);
    global_temp     = arena_alloc(16*GB, MB);
    cache.backing   = arena_alloc(16*GB, MB);

    if(std::atexit(release_global_persist)) { std::cerr << "Failed to register dealloc fn" << std::endl; return -1; }
    if(std::atexit(release_global_temp))    { std::cerr << "Failed to register dealloc fn" << std::endl; return -1; }
    if(std::atexit(release_cache_backing))  { std::cerr << "Failed to register dealloc fn" << std::endl; return -1; }

    std::string matrix1_name    = argv[1];
    std::string matrix2_name    = argv[2];
    std::string config_filepath = argv[3];

    std::ifstream file(config_filepath);
    if (!file.is_open()) {
        std::cerr << "Error opening config file." << std::endl;
        return 1;
    }

    json config;
    file >> config;

    int transpose           = config["transpose"].get<int>();
    f32 tmpsram             = config["cachesize"].get<f32>();
    f32 tmpbandw            = config["memorybandwidth"].get<f32>();
    int baselinetest        = config["baselinetest"].get<int>();
    bool condensedOP        = config["condensedOP"].get<bool>();
    std::string tile_dir    = config["tileDir"].get<std::string>();
    std::string output_dir  = config["outputDir"].get<std::string>();
    enum workload_mode requested_workload_mode = parse_workload_mode(config);
    std::string dense_matrices = config["denseMatrix"].get<std::string>(); // "A", "B", "both", "neither", default neither

    cachesize = tmpsram * 262144 * 0.9;
    inputcachesize = cachesize;
    HBMbandwidth = (tmpbandw / 4.0) * 0.6;
    int tmpPE = config["PEcnt"].get<int>();
    PEcnt = tmpPE;
    mergecnt = tmpPE;
    HBMbandwidthperPE = HBMbandwidth / PEcnt;
    int tmpbank = config["srambank"].get<int>();
    sramBank = tmpbank;

    std::string matrix1_filepath = get_matrix_path(matrix1_name);
    std::string matrix2_filepath = get_matrix_path(matrix2_name);
    std::string output_filepath = output_dir
        +   (1 ? "C" : "_") // +   (ISCACHE ? "C" : "_")
        +   printDataFlow[dataflow]
        +   (baselinetest ? "Base_" : "570Cache_")
        +   std::to_string(tmpsram)
        +   "MB_" + std::to_string(tmpbandw)
        +   "GBs_" + std::to_string(tmpPE)
        +   "PEs_" + std::to_string(tmpbank) + "sbanks_"
        +   "_" + matrix1_name + "_" + matrix2_name + "_"
        +   printFormat[format] + "_" + (transpose ? "1" : "0") 
	+   "_dense_" + dense_matrices
	+ ".txt";

    FILE *matrix1_file  = fopen(matrix1_filepath.c_str(), "r");
    FILE *matrix2_file  = fopen(matrix2_filepath.c_str(), "r");
    assert(matrix1_file);
    assert(matrix2_file);
    assert(freopen(output_filepath.c_str(), "w", stdout));

    struct matrix matA = {0};
    struct matrix matB = {0};
    {
        #define TEMP_BUFFER_NBYTES 1024
        char buf[TEMP_BUFFER_NBYTES];

        // read and ignore annotation '%' lines
        while (fgets(buf, TEMP_BUFFER_NBYTES, matrix1_file)) {
            if (buf[0] != '%')
                break;
        }
        sscanf(buf, "%u%u%u", &matA.nrows, &matA.ncols, &matA.nzM);

        // read and ignore annotation '%' lines
        while (fgets(buf, TEMP_BUFFER_NBYTES, matrix2_file)) {
            if (buf[0] != '%')
                break;
        }
        sscanf(buf, "%u%u%u", &matB.nrows, &matB.ncols, &matB.nzM);
    }

    matA.dense = dense_matrices == "A" || dense_matrices == "both";
    matB.dense = dense_matrices == "B" || dense_matrices == "both";
    if(matA.dense) matA.nzM = matA.nrows * matA.ncols;
    if(matB.dense) matB.nzM = matB.nrows * matB.ncols;
    // matA.dense = (b16)((matA.nrows*matA.ncols) == matA.nzM);
    // matB.dense = (b16)((matB.nrows*matB.ncols) == matB.nzM);
    matA.transpose = (b16)transpose;
    matB.transpose = (b16)((matB.nrows == matB.ncols) ? transpose : !transpose);
    if (matA.transpose)
        swap(matA.ncols, matA.nrows);
    if (matB.transpose)
        swap(matB.ncols, matB.nrows);
    assert(matA.ncols == matB.nrows);

    FILE *tile_file = fopen((tile_dir + matrix1_name).c_str(), "r");
    assert(tile_file);
    Coord t_i, t_j, t_k;
    {
        #define TEMP_BUFFER_NBYTES 1024
        char buf[TEMP_BUFFER_NBYTES];

        assert(fgets(buf, TEMP_BUFFER_NBYTES, tile_file));
        assert(sscanf(buf, "%u%u%u", &t_i, &t_j, &t_k) == 3);
    }

    const struct config cfg = {
        .dataflow   = Gust,
        .interorder = IJK,
        .format     = RR,
        .workload_mode = requested_workload_mode,

        .I          = matA.nrows,
        .J          = matA.ncols,
        .K          = matB.ncols,
        .tti        = div_rup(matA.nrows, t_i),
        .ttj        = div_rup(matA.ncols, t_j),
        .ttk        = div_rup(matB.ncols, t_k),
        .iii        = t_i,
        .jjj        = t_j,
        .kkk        = t_k,
    };
    sim = initialize_simulator(&cfg);
    
    parse_matrix(matrix1_file, &matA);
    A   = matA.M;
    Ac  = matA.Mc;
    offsetarrayA = matA.offsetarrayM;
    offsetarrayAc= matA.offsetarrayMc;

    if (condensedOP) {
        assert(0); // this path is not maintained

        // // memory management for sparchA
        // sparchA = new std::vector<int>[sim.cfg.J]();
        // if (sparchA == nullptr) {
        //     if (sparchA != nullptr)
        //         delete[] sparchA;
        //     std::cerr << "Error allocating memory for sparchA" << std::endl;
        //     std::exit(1);
        // }

        // // if use the condensed OP dataflow, need to preprocess the A matrix into
        // // the condensed format first. first put the data into sparchA[], then put
        // // it back to A[], and call gust dataflow
        // for (int j = 0; j < sim.cfg.J; j++) {
        //     for (int i = 0; i < sim.cfg.I; i++) {
        //         if (static_cast<int>(A[i].size()) > j) {
        //             sparchA[j].push_back(A[i][j]);
        //         }
        //     }
        // }

        // for (int j = 0; j < sim.cfg.J; j++) {
        //     A[j].clear();
        //     for (int i = 0; i < static_cast<int>(sparchA[j].size()); i++) {
        //         A[j].push_back(sparchA[j][i]);
        //     }
        // }

        // delete[] sparchA;
    }

    long long totalempty = 0;

    long long totalincache = 0;

    long long totaltagmatch48 = 0;
    long long totaltagmatch16 = 0;

    for (uint32_t i = 1; i < sim.cfg.I; i++) {
        Coord size_im1 = offsetarrayA[i] - offsetarrayA[i-1];
        if (size_im1 < 48) {
            totalempty += (48 - size_im1);
        }
        totalincache    += min(48u, size_im1);
        totaltagmatch48 += div_rup(size_im1, 48);
        totaltagmatch16 += div_rup(size_im1, 16);
    }

    /////////////////////////// input B /////////////////////////////////
    parse_matrix(matrix2_file, &matB);
    B   = matB.M;
    Bc  = matB.Mc;
    offsetarrayB = matB.offsetarrayM;
    offsetarrayBc= matB.offsetarrayMc;
    sim.workload = characterize_workload(requested_workload_mode, &matA, &matB, &sim.cfg);


    // {
    //     struct matrix_tiling_config tile_cfg = {.a = global_persist};
    //     struct matrix_tiling_luts   luts;
    //     b16 A_row_major = (dataflow == Inner) || (dataflow == Gust);
    //     tile_cfg        = A_row_major
    //         ? (struct matrix_tiling_config) {
    //             .major_dim      = sim.cfg.I,
    //             .minor_dim      = sim.cfg.J,
    //             .minor_tile_dim = sim.cfg.jjj,
    //             .map            = A,
    //             .offsets        = offsetarrayA,
    //         }
    //         : (struct matrix_tiling_config) {
    //             .major_dim      = sim.cfg.J,
    //             .minor_dim      = sim.cfg.I,
    //             .minor_tile_dim = sim.cfg.iii,
    //             .map            = Ac,
    //             .offsets        = offsetarrayAc,
    //         };
    //     luts = create_tiling_luts(&tile_cfg);
    //     sim.beginsA = luts.begins;
    //     sim.sizesA  = luts.sizes;

    //     b16 B_row_major = (dataflow == Outer) || (dataflow == Gust);
    //     tile_cfg        = B_row_major
    //         ? (struct matrix_tiling_config) {
    //             .major_dim      = sim.cfg.J,
    //             .minor_dim      = sim.cfg.K,
    //             .minor_tile_dim = sim.cfg.kkk,
    //             .map            = B,
    //             .offsets        = offsetarrayB,
    //         }
    //         : (struct matrix_tiling_config) {
    //             .major_dim      = sim.cfg.K,
    //             .minor_dim      = sim.cfg.J,
    //             .minor_tile_dim = sim.cfg.jjj,
    //             .map            = Bc,
    //             .offsets        = offsetarrayBc,
    //         };
    //     luts = create_tiling_luts(&tile_cfg);
    //     sim.beginsB = luts.begins;
    //     sim.sizesB  = luts.sizes;
    // }

    // WARNING: hardcoded cursor to IJK (tile-level), Gust (in-tile)
    {
        sim.cursor = (struct cursor) {
            .outer_dim      = sim.cfg.I,
            .middle_dim     = sim.cfg.J,
            .inner_dim      = sim.cfg.K,
            .outer_tile_dim = sim.cfg.iii,
            .middle_tile_dim= sim.cfg.jjj,
            .inner_tile_dim = sim.cfg.kkk,
            .outer_ntiles   = sim.cfg.tti,
            .middle_ntiles  = sim.cfg.ttj,
            .inner_ntiles   = sim.cfg.ttk,

            // .ti             = 0,
            // .tj             = 0,
            // .tk             = 0,

            .outer_tile_idx = &sim.cursor.ti,
            .middle_tile_idx= &sim.cursor.tj,
            .inner_tile_idx = &sim.cursor.tk,

            .outer_wrap     = &sim.cursor.wrap_ti, // outer should never wrap
            .middle_wrap    = &sim.cursor.wrap_tj,
            .inner_wrap     = &sim.cursor.wrap_tk
        };

        sim.cursor.A    = (struct tile) {
            .major_dim      = sim.cfg.iii,
            .minor_dim      = sim.cfg.jjj,
            .map            = (const Coord **)A,
            .offsets        = (const Coord *)offsetarrayA,
            .minor_wrap     = &sim.cursor.wrap_tj,
            .major_tile_idx = &sim.cursor.ti,
            .minor_tile_idx = &sim.cursor.tj
        };

        sim.cursor.B    = (struct tile) {
            .major_dim      = sim.cfg.jjj,
            .minor_dim      = sim.cfg.kkk,
            .map            = (const Coord **)B,
            .offsets        = (const Coord *)offsetarrayB,
            .minor_wrap     = &sim.cursor.wrap_tk,
            .major_tile_idx = &sim.cursor.tj,
            .minor_tile_idx = &sim.cursor.tk
        };

        sim.cursor.A.begins = (Coord *)arena_push(global_persist, sim.cursor.A.major_dim*sizeof(Coord), __alignof__(Coord), 1);
        sim.cursor.A.sizes  = (Coord *)arena_push(global_persist, sim.cursor.A.major_dim*sizeof(Coord), __alignof__(Coord), 1);
        sim.cursor.B.begins = (Coord *)arena_push(global_persist, sim.cursor.B.major_dim*sizeof(Coord), __alignof__(Coord), 1);
        sim.cursor.B.sizes  = (Coord *)arena_push(global_persist, sim.cursor.B.major_dim*sizeof(Coord), __alignof__(Coord), 1);
    }



    /******************Config************************************/

    printf("Matrix A: %u x %u, number of non-zeros = %u\n", matA.nrows, matA.ncols, matA.nzM);
    printf("*** ratio of empty %lf, ratio of not empty %lf\n", totalempty / (sim.cfg.I * 48.0), 1 - (totalempty / (sim.cfg.I * 48.0)));
    printf("*** ratio of in cache %lf\n", totalincache / ((f64)matA.nzM));
    printf("** ratio tag access 48 %lf\n", sim.cfg.I / ((f64)sim.cfg.I + totaltagmatch48));
    printf("** ratio tag access 16 %lf\n", sim.cfg.I / ((f64)sim.cfg.I + totaltagmatch16));
    printf("Matrix B: %u x %u, number of non-zeros = %u\n", matB.nrows, matB.ncols, matB.nzM);
    printf("transpose: %d\n", matB.transpose);
    printf("I = %llu, K = %llu, J = %llu\n", sim.cfg.I, sim.cfg.K, sim.cfg.J);
    printf("workload mode requested = %s, selected = %s\n",
        print_workload_mode(sim.workload.requested_mode),
        print_workload_mode(sim.workload.selected_mode));
    printf("workload decision = %s\n", sim.workload.decision_reason.c_str());
    printf("A density = %.6lf, mean_row_nnz = %.2lf, row_nnz_stddev = %.2lf, empty_row_ratio = %.6lf, tile_fill_ratio = %.6lf\n",
        sim.workload.A.density,
        sim.workload.A.mean_row_nnz,
        sim.workload.A.row_nnz_stddev,
        sim.workload.A.empty_row_ratio,
        sim.workload.A.tile_fill_ratio);
    printf("B density = %.6lf, mean_row_nnz = %.2lf, row_nnz_stddev = %.2lf, empty_row_ratio = %.6lf, tile_fill_ratio = %.6lf\n",
        sim.workload.B.density,
        sim.workload.B.mean_row_nnz,
        sim.workload.B.row_nnz_stddev,
        sim.workload.B.empty_row_ratio,
        sim.workload.B.tile_fill_ratio);
    printf("tile mode counts: sparse = %llu, dense = %llu\n",
        sim.workload.sparse_tile_count,
        sim.workload.dense_tile_count);
    if (sim.workload.policy_region_span_ti || sim.workload.policy_region_span_tj || sim.workload.policy_region_span_tk) {
        printf("policy region span (tiles): ti = %u, tj = %u, tk = %u; region counts: sparse = %u, dense = %u\n",
            sim.workload.policy_region_span_ti,
            sim.workload.policy_region_span_tj,
            sim.workload.policy_region_span_tk,
            sim.workload.sparse_region_count,
            sim.workload.dense_region_count);
    } else if (sim.workload.sparse_region_count || sim.workload.dense_region_count) {
        printf("policy region grouping = adaptive similarity-based; region counts: sparse = %u, dense = %u\n",
            sim.workload.sparse_region_count,
            sim.workload.dense_region_count);
    }
    if (sim.workload.estimated_sparse_cost || sim.workload.estimated_dense_cost || sim.workload.estimated_mixed_cost) {
        printf("estimated tile costs: sparse = %.2lf, dense = %.2lf, mixed = %.2lf\n",
            sim.workload.estimated_sparse_cost,
            sim.workload.estimated_dense_cost,
            sim.workload.estimated_mixed_cost);
        printf("estimated oracle tile cost = %.2lf, estimated policy overhead = %.2lf, policy build time = %.2lf us\n",
            sim.workload.estimated_oracle_tile_cost,
            sim.workload.estimated_policy_overhead_cost,
            sim.workload.policy_build_time_us);
    }
    /************************************************************/

    getParameter(); // sets estEffMAC

    configPartial(0.05, 0.5, 0.45);

    /////////////// Baseline configurations

    if (baselinetest) {
        // EWH
        // Incorporate SeaCache into baseline
        puts("***************** SeaCache *******************");
        printf("nnzB:%u  K:%u  J/TJ:%u  nzlB:%u\n",
            matB.nzM, sim.cfg.K, (sim.cfg.J + sim.cfg.jjj - 1) / sim.cfg.jjj,
            matB.nzM / (sim.cfg.K * ((sim.cfg.J + sim.cfg.jjj - 1) / sim.cfg.jjj))
        );

        ISCACHE = 1;
        cachesize               = 262144;
        CACHE_BLOCK_NELEMS      = 16;
        CACHE_BLOCK_NELEMS_LOG2 = getlog(CACHE_BLOCK_NELEMS);
        setSET();

        cache.cfg = {
            .block_nelems       = 1,
            .block_nelems_log2  = 1,
            .scheme             = CACHE_SCHEME_FLFU,
        };

        adaptive_prefetch = 1;
        useVirtualTag = 1;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;

        reset_cursor(&sim.cursor);
        runTile(sim.cfg.kkk);

        adaptive_prefetch = 0;
        useVirtualTag = 0;

        adaptive_prefetch = 0;

        ////////////  InnserSP
        // static FLRU + 16 words scheme0
        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   test InnerSP   "
             "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        prefetchSize = inputcachesize / 6;
        cacheScheme = CACHE_SCHEME_INNER_SP;
        cachesize = inputcachesize;
        CACHE_BLOCK_NELEMS = 16;
        CACHE_BLOCK_NELEMS_LOG2 = 4;
        setSET();
        reset_cursor(&sim.cursor);
        runTile(sim.cfg.kkk);

        ////////////  Sparch
        // dynamic FLRU + 128KB prefetch size + 144 words scheme0
        puts("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   test Sparch   "
             "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        ISCACHE = 1;
        cacheScheme = CACHE_SCHEME_SPARCH;
        prefetchSize = inputcachesize / 6;
        cachesize = inputcachesize - prefetchSize;
        CACHE_BLOCK_NELEMS = 144;
        CACHE_BLOCK_NELEMS_LOG2 = 8;
        setSET();
        // calculate metadata overhead.
        // if metadata overflow, choose smaller tile
        int newkkk = sim.cfg.kkk;
        int newttk = sim.cfg.ttk;
        // if can keep, just use current kkk
        if (cachesize > sim.cfg.kkk * 2) {
            cachesize -= sim.cfg.kkk * 2;
        } else {
            // if can't keep, use smaller kkk
            // (make kkk*2 to be half cachesize)
            newkkk = cachesize / 4;
            newttk = (sim.cfg.K + sim.cfg.kkk - 1) / sim.cfg.kkk;
            cachesize -= sim.cfg.kkk * 2;
        }
        reset_cursor(&sim.cursor);
        runTile(newkkk);
        // return to the default setting
        CACHE_BLOCK_NELEMS = 16;
        CACHE_BLOCK_NELEMS_LOG2 = 4;
        cachesize = inputcachesize;
        setSET();

        ////////////  X-cache
        puts("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   test X-cache   "
             "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        // LRU + 4 words scheme0
        // just same as using scheme0 with cacheline = 4
        ISCACHE = 1;
        cacheScheme = CACHE_SCHEME_BASE;
        cachesize = inputcachesize;
        CACHE_BLOCK_NELEMS = 4;
        CACHE_BLOCK_NELEMS_LOG2 = 2;
        setSET();
        reset_cursor(&sim.cursor);
        runTile(sim.cfg.kkk);

        // return to the default setting
        CACHE_BLOCK_NELEMS = 16;
        CACHE_BLOCK_NELEMS_LOG2 = 4;
        setSET();

        puts("!!!!!!!!!!!!!!!!!!!!  Scratchpad   !!!!!!!!!!!!!!!!!!!!!!!");
        ISCACHE = 0;

        configPartial(0.05, 0.9, 0.05);

        reinitialize();

        reset_cursor(&sim.cursor);
        run();
    }

    if (!baselinetest) {
        puts("\n!!!!!!!!!!!!!!!!!!!! EECS570 !!!!!!!!!!!!!!!!!!!!");
        // cribbed settings from SCACHE
        ISCACHE = 1;
        cachesize               = 262144;
        CACHE_BLOCK_NELEMS      = 16;
        CACHE_BLOCK_NELEMS_LOG2 = getlog(CACHE_BLOCK_NELEMS);
        setSET();

        cache.cfg = {
            .block_nelems       = 1,
            .block_nelems_log2  = 1,
            .scheme             = CACHE_SCHEME_FLFU,
        };

        adaptive_prefetch = 1;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;

        // adaptive sparse-dense scheme, also uses virtual tag in sparse mode. Overloading this
        // a little but it's ok.
        useVirtualTag = 2;

        reset_cursor(&sim.cursor);
        runTile(sim.cfg.kkk);

        printf("Dense installations: %lld\n", totalDenseInstalls);
        printf("Dense hits: %lld\n", totalDenseHits);
    }

    bool ablationtest = 0;
    if (ablationtest) {

        adaptive_prefetch = 0;

        /////////////// ablation test

        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!! scheme0 (base)   "
             "!!!!!!!!!!!!!!!!!!!!!!!!");
        puts("CacheScheme 0");
        ISCACHE = 1;
        cacheScheme = CACHE_SCHEME_BASE;
        cachesize = inputcachesize;
        setSET();
        runTile(sim.cfg.kkk);

        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!! scheme1 (mapping)   "
             "!!!!!!!!!!!!!!!!!!!!!!!!");

        puts("CacheScheme 1");
        ISCACHE = 1;
        cacheScheme = CACHE_SCHEME_MAPPING;
        cachesize = inputcachesize;
        setSET();
        runTile(sim.cfg.kkk);

        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!! scheme88 without virtue   "
             "!!!!!!!!!!!!!!!!!!!!!!!!");

        useVirtualTag = 0;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;
        prefetchSize = cachesize / 6;
        runTile(sim.cfg.kkk);

        puts("\n!!!!!!!!!!!!!!!!!!!!!!!!!! scheme88 with virtue   "
             "!!!!!!!!!!!!!!!!!!!!!!!!");

        puts("CacheScheme 88 practical FLFU  with virtual tag 1/6");
        useVirtualTag = 1;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;
        prefetchSize = cachesize / 6;
        runTile(sim.cfg.kkk);
        useVirtualTag = 0;

        puts("CacheScheme 88 practical FLFU  with virtual tag 1/16");
        useVirtualTag = 1;
        cacheScheme = CACHE_SCHEME_FLFU;
        cachesize = inputcachesize;
        prefetchSize = cachesize / 16;
        runTile(sim.cfg.kkk);
        useVirtualTag = 0;
    }

    return 0;
}

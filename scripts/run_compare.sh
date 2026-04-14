#!/bin/bash
# Run matrices under sparse / dense / auto modes and collect results.
# Dense mode is skipped for very large matrices (> NNZ_DENSE_LIMIT) to avoid
# impractical runtimes.
set -e

BINARY=./scache
BASE_CONFIG=config/config.json
CONFIG_DIR=config

NNZ_DENSE_LIMIT=500000   # skip forced-dense above this nnz

mkdir -p "$CONFIG_DIR"

# Generate per-mode config files
for mode in sparse dense auto; do
    python3 -c "
import json
with open('$BASE_CONFIG') as f:
    cfg = json.load(f)
cfg['workloadMode'] = '$mode'
with open('$CONFIG_DIR/config_$mode.json', 'w') as f:
    json.dump(cfg, f, indent=4)
"
done

run_pair() {
    local mat1="$1"
    local mat2="$2"
    local label="$3"

    [ -f "data/${mat1}.mtx" ] || { echo "  SKIP (no mtx): $mat1"; return; }
    [ -f "tiles/${mat1}" ]    || { echo "  SKIP (no tile): $mat1"; return; }
    [ -f "data/${mat2}.mtx" ] || { echo "  SKIP (no mtx): $mat2"; return; }

    # Get nnz from header line
    local nnz
    nnz=$(grep -v "^%" "data/${mat1}.mtx" | head -1 | awk '{print $3}')

    echo "  [sparse] $label"
    "$BINARY" "$mat1" "$mat2" "$CONFIG_DIR/config_sparse.json" 2>/dev/null \
        || echo "  ERROR sparse: $label"

    if [ "$nnz" -le "$NNZ_DENSE_LIMIT" ]; then
        echo "  [dense]  $label"
        "$BINARY" "$mat1" "$mat2" "$CONFIG_DIR/config_dense.json" 2>/dev/null \
            || echo "  ERROR dense: $label"
    else
        echo "  [dense]  SKIPPED (nnz=$nnz > limit)"
    fi

    echo "  [auto]   $label"
    "$BINARY" "$mat1" "$mat2" "$CONFIG_DIR/config_auto.json" 2>/dev/null \
        || echo "  ERROR auto: $label"
}

# ── Synthetic matrices (256×256, designed to stress mixed-mode) ──────────────
SYNTH=(
    synth_dense_full
    synth_dense_band
    synth_mixed_checker
    synth_mixed_quadrant
    synth_mixed_stripes
    synth_clustered_blocks
    synth_sparse_random
    synth_sparse_diag
    dense256_tmp
    mixed256_tmp
)

echo "=== Synthetic matrices ==="
for name in "${SYNTH[@]}"; do
    echo "--- $name ---"
    run_pair "$name" "$name" "$name"
done

# ── Real sparse matrices (self-mult) ─────────────────────────────────────────
REAL_SQUARE=(
    filter3D
    ship_001
    amazon0601
    web-Google
)

echo ""
echo "=== Real sparse matrices ==="
for name in "${REAL_SQUARE[@]}"; do
    echo "--- $name ---"
    run_pair "$name" "$name" "$name"
done

# ── LLM weight matrices ───────────────────────────────────────────────────────
LLM_SQUARE=(
    gpt2_l0_attn_out_sp90
    bert_l0_q_sp90
    bert_l0_out_sp90
)
LLM_RECT_PAIRS=(
    "gpt2_l0_mlp_fc_sp90  gpt2_l0_mlp_proj_sp90  gpt2_l0_mlp_fcXproj"
    "bert_l0_ffn1_sp90    bert_l0_ffn2_sp90       bert_l0_ffn1Xffn2"
)

echo ""
echo "=== LLM weight matrices ==="
for name in "${LLM_SQUARE[@]}"; do
    echo "--- $name ---"
    run_pair "$name" "$name" "$name"
done
for pair in "${LLM_RECT_PAIRS[@]}"; do
    read -r m1 m2 lbl <<< "$pair"
    echo "--- $lbl ---"
    run_pair "$m1" "$m2" "$lbl"
done

echo ""
echo "All runs complete. Output files in: output/"

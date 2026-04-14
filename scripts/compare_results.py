"""
Parse SeaCache-sim output files and produce a comparison table:
  - Pure Sparse mode (forced)
  - Pure Dense mode (forced)
  - Auto (tile-level detect-and-switch)

Metrics extracted per run:
  - SeaCache total cycle (our algorithm)
  - workload decision / selected mode
  - tile mode counts (sparse/dense)
  - density
  - hitrate
"""

import os
import re
import sys
from collections import defaultdict

OUTPUT_DIR = "output"

# ------------------------------------------------------------------ #
# Parsing helpers
# ------------------------------------------------------------------ #

def parse_output_file(path):
    """Return a dict of fields parsed from one simulator output file."""
    result = {}
    with open(path) as f:
        text = f.read()

    def first(pattern, conv=str):
        m = re.search(pattern, text)
        return conv(m.group(1)) if m else None

    result["nrows"]          = first(r"Matrix A: (\d+) x", int)
    result["ncols"]          = first(r"Matrix A: \d+ x (\d+)", int)
    result["nnzA"]           = first(r"number of non-zeros = (\d+)", int)
    result["density_A"]      = first(r"A density = ([\d.]+)", float)
    result["selected_mode"]  = first(r"workload mode requested = \w+, selected = (\w+)")
    result["decision"]       = first(r"workload decision = (.+)")
    result["sparse_tiles"]   = first(r"tile mode counts: sparse = (\d+)", int)
    result["dense_tiles"]    = first(r"tile mode counts: sparse = \d+, dense = (\d+)", int)
    result["est_sparse"]     = first(r"estimated tile costs: sparse = ([\d.]+)", float)
    result["est_dense"]      = first(r"estimated tile costs: sparse = [\d.]+, dense = ([\d.]+)", float)
    result["est_mixed"]      = first(r"estimated tile costs: .+ mixed = ([\d.]+)", float)

    # SeaCache block: the FIRST "total cycle" after the SeaCache header
    seacache_block = re.split(r"!!!|Scratchpad", text)[0]
    seacache_block = re.split(r"\*+ SeaCache \*+", seacache_block, maxsplit=1)
    if len(seacache_block) > 1:
        m = re.search(r"total cycle = (\d+)", seacache_block[1])
        result["seacache_cycles"] = int(m.group(1)) if m else None
        m = re.search(r"hitrate = ([\d.]+)", seacache_block[1])
        result["seacache_hitrate"] = float(m.group(1)) if m else None

    return result


def find_output_files():
    """
    Returns dict keyed by (matrix_label, mode) → filepath.
    Filename pattern:
      CGustBase_<hw>__<mat1>_<mat2>_RR_<transpose>.txt
    We encode mode in the config used, but the filename doesn't include mode.
    We work around this by re-running outputs into mode-specific subdirs.
    """
    files = {}
    for fname in os.listdir(OUTPUT_DIR):
        if not fname.endswith(".txt"):
            continue
        path = os.path.join(OUTPUT_DIR, fname)
        # Determine mode from inside the file
        with open(path) as f:
            text = f.read()
        m = re.search(r"workload mode requested = (\w+)", text)
        mode = m.group(1) if m else "unknown"

        # Extract matrix label from filename:  __<mat1>_<mat2>_RR
        m = re.search(r"__(.+)_RR_\d+\.txt$", fname)
        label = m.group(1) if m else fname

        key = (label, mode)
        # If duplicate (re-run), keep latest
        if key not in files or os.path.getmtime(path) > os.path.getmtime(files[key]):
            files[key] = path
    return files


# ------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------ #

def main():
    files = find_output_files()

    # Group by matrix label
    labels = sorted(set(label for label, mode in files))

    MODES = ["sparse", "dense", "auto"]
    COL_W = 16

    # Header
    header = f"{'Matrix':<45}" + "".join(f"{'['+m+'] cycles':>{COL_W}}" for m in MODES) + \
             f"{'auto/sparse':>{COL_W}}" + f"{'auto/dense':>{COL_W}}" + \
             f"{'selected':>{COL_W}}" + f"{'sp_tiles':>{COL_W}}" + f"{'dn_tiles':>{COL_W}}" + \
             f"{'density':>{COL_W}}"
    print(header)
    print("-" * len(header))

    for label in labels:
        row_data = {}
        for mode in MODES:
            key = (label, mode)
            if key in files:
                row_data[mode] = parse_output_file(files[key])

        if not row_data:
            continue

        cycles = {m: row_data[m].get("seacache_cycles") for m in MODES if m in row_data}
        sp = cycles.get("sparse")
        dn = cycles.get("dense")
        au = cycles.get("auto")

        auto_vs_sparse = f"{au/sp:.3f}x" if (au and sp) else "N/A"
        auto_vs_dense  = f"{au/dn:.3f}x" if (au and dn) else "N/A"

        ref = row_data.get("auto") or row_data.get("sparse") or {}
        selected   = ref.get("selected_mode", "?")
        sp_tiles   = ref.get("sparse_tiles", "?")
        dn_tiles   = ref.get("dense_tiles", "?")
        density    = ref.get("density_A")
        density_s  = f"{density:.5f}" if density else "?"

        def fmt(v):
            return f"{v:>{COL_W},}" if v else f"{'N/A':>{COL_W}}"

        row = f"{label:<45}" + fmt(sp) + fmt(dn) + fmt(au) + \
              f"{auto_vs_sparse:>{COL_W}}" + f"{auto_vs_dense:>{COL_W}}" + \
              f"{selected:>{COL_W}}" + f"{sp_tiles!s:>{COL_W}}" + \
              f"{dn_tiles!s:>{COL_W}}" + f"{density_s:>{COL_W}}"
        print(row)

    print()
    print("Notes:")
    print("  auto/sparse < 1.0x  =>  auto is FASTER than pure sparse (switch to dense helps)")
    print("  auto/dense  < 1.0x  =>  auto is FASTER than pure dense  (switch to sparse helps)")
    print("  auto/sparse > 1.0x  =>  auto is slower than pure sparse  (overhead or wrong switch)")


if __name__ == "__main__":
    main()

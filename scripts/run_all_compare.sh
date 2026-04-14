#!/bin/bash
# Run all runnable matrices under sparse/dense/auto workload modes
# and collect output files for comparison.
set -e

BINARY=./scache
BASE_CONFIG=config/config.json
CONFIG_DIR=config
OUT_DIR=output

mkdir -p "$CONFIG_DIR" "$OUT_DIR"

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

# -----------------------------------------------------------------
# Square matrices: self-multiplication A x A
# -----------------------------------------------------------------
SQUARE_MATRICES=(
    amazon0601
    filter3D
    ship_001
    web-Google
    # LLM 768x768 square weights
    gpt2_l0_attn_out_sp90
    gpt2_l1_attn_out_sp90
    bert_l0_q_sp90
    bert_l0_k_sp90
    bert_l0_v_sp90
    bert_l0_out_sp90
    bert_l1_q_sp90
    bert_l1_k_sp90
    bert_l1_v_sp90
    bert_l1_out_sp90
)

# -----------------------------------------------------------------
# Non-square compatible pairs: A x B where A.ncols == B.nrows
# tile file must exist for matrix1; mtx must exist for both
# -----------------------------------------------------------------
# Format: "mat1 mat2 label"
PAIRS=(
    "gpt2_l0_mlp_fc_sp90   gpt2_l0_mlp_proj_sp90  gpt2_l0_mlp_fcXproj"
    "gpt2_l1_mlp_fc_sp90   gpt2_l1_mlp_proj_sp90  gpt2_l1_mlp_fcXproj"
    "bert_l0_ffn1_sp90     bert_l0_ffn2_sp90       bert_l0_ffn1Xffn2"
    "bert_l1_ffn1_sp90     bert_l1_ffn2_sp90       bert_l1_ffn1Xffn2"
)

run_matrix() {
    local mat1="$1"
    local mat2="$2"
    local mode="$3"
    local config="$CONFIG_DIR/config_${mode}.json"

    if [ ! -f "data/${mat1}.mtx" ]; then
        echo "  SKIP (no mtx): $mat1"
        return
    fi
    if [ ! -f "tiles/${mat1}" ]; then
        echo "  SKIP (no tile): $mat1"
        return
    fi
    if [ ! -f "data/${mat2}.mtx" ]; then
        echo "  SKIP (no mtx for mat2): $mat2"
        return
    fi

    echo "  Running [$mode] $mat1 x $mat2 ..."
    "$BINARY" "$mat1" "$mat2" "$config" 2>/dev/null || echo "  ERROR: $mat1 x $mat2 [$mode]"
}

echo "=== Running square matrices (A x A) ==="
for name in "${SQUARE_MATRICES[@]}"; do
    echo "--- $name ---"
    for mode in sparse dense auto; do
        run_matrix "$name" "$name" "$mode"
    done
done

echo ""
echo "=== Running compatible matrix pairs ==="
for pair in "${PAIRS[@]}"; do
    read -r mat1 mat2 label <<< "$pair"
    echo "--- $label ---"
    for mode in sparse dense auto; do
        run_matrix "$mat1" "$mat2" "$mode"
    done
done

echo ""
echo "Done. Output files in: $OUT_DIR/"

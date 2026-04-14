"""
Extract weight matrices from LLMs, apply magnitude pruning to create
sparse matrices, and export as .mtx + tile files for SeaCache-sim.

Usage: python3 scripts/extract_llm_matrices.py
"""

import os
import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite
import torch
from transformers import AutoModelForCausalLM, AutoConfig

DATA_DIR = "data"
TILES_DIR = "tiles"

# Tile sizes used by existing matrices (from tiles/amazon0601: 403394 12607 403394)
# We'll use a reasonable tile size for LLM matrices
TILE_I = 32
TILE_J = 32
TILE_K = 32

SPARSITY_LEVELS = [0.5, 0.9, 0.95]  # fraction of weights zeroed out


def to_sparse_mtx(dense_weight: np.ndarray, name: str, sparsity: float):
    """
    Magnitude-prune a dense weight matrix to target sparsity,
    write .mtx and tile config files.
    Returns (nnz, nrows, ncols)
    """
    # Flatten to 2D
    if dense_weight.ndim > 2:
        rows = dense_weight.shape[0]
        cols = int(np.prod(dense_weight.shape[1:]))
        dense_weight = dense_weight.reshape(rows, cols)

    # Magnitude pruning: zero out smallest |w|
    flat = np.abs(dense_weight).flatten()
    threshold = np.percentile(flat, sparsity * 100)
    mask = np.abs(dense_weight) >= threshold
    sparse = dense_weight * mask

    csr = sp.csr_matrix(sparse)
    nnz = csr.nnz
    nrows, ncols = csr.shape
    print(f"  {name}: {nrows}x{ncols}, nnz={nnz}, density={nnz/(nrows*ncols):.4f}")

    mtx_path = os.path.join(DATA_DIR, f"{name}.mtx")
    with open(mtx_path, "w") as f:
        f.write(f"%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{nrows} {ncols} {nnz}\n")
        cx = csr.tocoo()
        for r, c in zip(cx.row, cx.col):
            # MTX is 1-indexed
            f.write(f"{r+1} {c+1}\n")

    # Tile config: nrows nzM ncols  (same format as tiles/amazon0601)
    tile_path = os.path.join(TILES_DIR, name)
    with open(tile_path, "w") as f:
        f.write(f"{TILE_I} {TILE_J} {TILE_K}")

    return nnz, nrows, ncols


def extract_from_gpt2():
    print("Loading GPT-2 (small)...")
    config = AutoConfig.from_pretrained("gpt2")
    model = AutoModelForCausalLM.from_config(config)  # random weights, no download needed

    matrices = []
    for layer_idx in range(min(2, config.n_layer)):  # first 2 layers
        block = model.transformer.h[layer_idx]

        # Attention projections: Q/K/V combined (n_embd x 3*n_embd) and output (n_embd x n_embd)
        attn_qkv = block.attn.c_attn.weight.detach().numpy()   # (n_embd, 3*n_embd)
        attn_out = block.attn.c_proj.weight.detach().numpy()   # (n_embd, n_embd)

        # MLP projections: (n_embd, 4*n_embd) and (4*n_embd, n_embd)
        mlp_fc   = block.mlp.c_fc.weight.detach().numpy()     # (n_embd, 4*n_embd)
        mlp_proj = block.mlp.c_proj.weight.detach().numpy()   # (4*n_embd, n_embd)

        matrices.append((f"gpt2_l{layer_idx}_attn_qkv", attn_qkv))
        matrices.append((f"gpt2_l{layer_idx}_attn_out", attn_out))
        matrices.append((f"gpt2_l{layer_idx}_mlp_fc",   mlp_fc))
        matrices.append((f"gpt2_l{layer_idx}_mlp_proj", mlp_proj))

    return matrices


def extract_from_bert():
    print("Loading BERT-base...")
    from transformers import AutoModel
    config = AutoConfig.from_pretrained("bert-base-uncased")
    model = AutoModel.from_config(config)

    matrices = []
    for layer_idx in range(min(2, config.num_hidden_layers)):
        layer = model.encoder.layer[layer_idx]

        q = layer.attention.self.query.weight.detach().numpy()
        k = layer.attention.self.key.weight.detach().numpy()
        v = layer.attention.self.value.weight.detach().numpy()
        out = layer.attention.output.dense.weight.detach().numpy()
        ffn1 = layer.intermediate.dense.weight.detach().numpy()
        ffn2 = layer.output.dense.weight.detach().numpy()

        matrices.append((f"bert_l{layer_idx}_q",    q))
        matrices.append((f"bert_l{layer_idx}_k",    k))
        matrices.append((f"bert_l{layer_idx}_v",    v))
        matrices.append((f"bert_l{layer_idx}_out",  out))
        matrices.append((f"bert_l{layer_idx}_ffn1", ffn1))
        matrices.append((f"bert_l{layer_idx}_ffn2", ffn2))

    return matrices


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(TILES_DIR, exist_ok=True)

    all_matrices = []
    all_matrices += extract_from_gpt2()
    all_matrices += extract_from_bert()

    # For SpGEMM (A x B), we pair each weight matrix with itself or another
    # Use sparsity=0.9 as the primary test (sparse region detection most relevant here)
    sparsity = 0.9
    names = []
    print(f"\nExporting matrices at sparsity={sparsity}:")
    for base_name, weight in all_matrices:
        name = f"{base_name}_sp{int(sparsity*100)}"
        nnz, nrows, ncols = to_sparse_mtx(weight, name, sparsity)
        names.append((name, nrows, ncols, nnz))

    print("\nGenerated matrix pairs for simulation:")
    print("=" * 60)
    for name, nrows, ncols, nnz in names:
        print(f"  ./scache {name} {name} config/config.json")

    # Write a run script
    with open("scripts/run_llm_sim.sh", "w") as f:
        f.write("#!/bin/bash\nset -e\n\n")
        f.write("echo '=== SeaCache LLM Matrix Simulation ==='\n")
        for name, nrows, ncols, nnz in names:
            # Only run square or compatible-dim matrices (A.ncols == B.nrows)
            f.write(f"echo '--- {name} ---'\n")
            f.write(f"./scache {name} {name} config/config.json\n")
        f.write("echo 'Done.'\n")
    os.chmod("scripts/run_llm_sim.sh", 0o755)
    print("\nRun script written to scripts/run_llm_sim.sh")


if __name__ == "__main__":
    main()

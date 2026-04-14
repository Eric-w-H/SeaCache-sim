#!/usr/bin/env python3
"""
Build real-workload matrices for SeaCache-sim.

Supported sources:
1) SNAP graph edge lists (.txt or .txt.gz) -> square adjacency .mtx
2) Longformer-style attention mask (window + global tokens) -> square .mtx

This script writes both:
- data/<name>.mtx
- tiles/<name>   (single-line tile sizes: "ti tj tk")
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
from pathlib import Path
from typing import Iterable, List, Set, Tuple


REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
TILES_DIR = REPO_ROOT / "tiles"


def _write_mtx_pattern(path: Path, nrows: int, ncols: int, coords: List[Tuple[int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{nrows} {ncols} {len(coords)}\n")
        for r, c in coords:
            # MatrixMarket is 1-based.
            f.write(f"{r + 1} {c + 1}\n")


def _write_tile(path: Path, ti: int, tj: int, tk: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{ti} {tj} {tk}\n")


def _pick_tile_size(n: int, target_tiles_per_dim: int = 32) -> int:
    # Keep enough tiles for mixed-mode decisions while avoiding tiny-tile overhead.
    val = max(16, n // max(1, target_tiles_per_dim))
    # Round to power-of-two for stable comparisons.
    p2 = 1 << max(4, int(math.log2(val)))
    return min(n, p2)


def ingest_snap_edge_list(
    input_path: Path,
    output_name: str,
    undirected: bool,
    target_tiles_per_dim: int,
    labels_path: Path | None,
) -> None:
    opener = gzip.open if input_path.suffix == ".gz" else open
    edges_raw: List[Tuple[int, int]] = []
    node_ids: Set[int] = set()

    with opener(input_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            u = int(parts[0])
            v = int(parts[1])
            edges_raw.append((u, v))
            node_ids.add(u)
            node_ids.add(v)

    if not node_ids:
        raise RuntimeError(f"no edges parsed from {input_path}")

    # Reindex to dense [0..N-1] for compact matrices.
    # If labels are provided, group nodes by label first to expose realistic
    # community blocks (useful for GNN adjacency experiments).
    label_map = {}
    if labels_path is not None:
        lopener = gzip.open if labels_path.suffix == ".gz" else open
        with lopener(labels_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                nid = int(parts[0])
                lbl = int(parts[1])
                label_map[nid] = lbl

    if label_map:
        sorted_nodes = sorted(node_ids, key=lambda nid: (label_map.get(nid, 10**12), nid))
    else:
        sorted_nodes = sorted(node_ids)
    id_map = {old: new for new, old in enumerate(sorted_nodes)}
    coords_set: Set[Tuple[int, int]] = set()
    for u_old, v_old in edges_raw:
        u = id_map[u_old]
        v = id_map[v_old]
        coords_set.add((u, v))
        if undirected and u != v:
            coords_set.add((v, u))

    n = len(id_map)
    coords = sorted(coords_set)
    density = len(coords) / float(n * n)

    mtx_path = DATA_DIR / f"{output_name}.mtx"
    tile_path = TILES_DIR / output_name
    _write_mtx_pattern(mtx_path, n, n, coords)

    tile = _pick_tile_size(n, target_tiles_per_dim=target_tiles_per_dim)
    _write_tile(tile_path, tile, tile, tile)

    print(
        f"[snap] wrote {output_name}: {n}x{n}, nnz={len(coords)}, "
        f"density={density:.6f}, tile={tile}"
    )
    print(f"       mtx:  {mtx_path}")
    print(f"       tile: {tile_path}")


def ingest_longformer_mask(
    output_name: str,
    seq_len: int,
    window: int,
    global_tokens: int,
    target_tiles_per_dim: int,
) -> None:
    if global_tokens > seq_len:
        raise ValueError("global_tokens must be <= seq_len")
    if window < 0:
        raise ValueError("window must be non-negative")

    gset = set(range(global_tokens))
    half = window // 2
    coords_set: Set[Tuple[int, int]] = set()

    # Attention layout (connectivity mask):
    # - local sliding window around each token
    # - global tokens attend to all and are attended by all
    for i in range(seq_len):
        lo = max(0, i - half)
        hi = min(seq_len - 1, i + half)
        for j in range(lo, hi + 1):
            coords_set.add((i, j))
        for g in gset:
            coords_set.add((i, g))
            coords_set.add((g, i))

    coords = sorted(coords_set)
    density = len(coords) / float(seq_len * seq_len)

    mtx_path = DATA_DIR / f"{output_name}.mtx"
    tile_path = TILES_DIR / output_name
    _write_mtx_pattern(mtx_path, seq_len, seq_len, coords)
    tile = _pick_tile_size(seq_len, target_tiles_per_dim=target_tiles_per_dim)
    _write_tile(tile_path, tile, tile, tile)

    print(
        f"[longformer] wrote {output_name}: {seq_len}x{seq_len}, nnz={len(coords)}, "
        f"density={density:.6f}, window={window}, global={global_tokens}, tile={tile}"
    )
    print(f"             mtx:  {mtx_path}")
    print(f"             tile: {tile_path}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Ingest real workloads into SeaCache format.")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_snap = sub.add_parser("snap", help="Convert SNAP edge list to .mtx + tile")
    p_snap.add_argument("--input", required=True, type=Path, help="Path to .txt or .txt.gz edge list")
    p_snap.add_argument("--name", required=True, help="Output matrix/tile base name")
    p_snap.add_argument("--undirected", action="store_true", help="Mirror edges")
    p_snap.add_argument(
        "--labels",
        type=Path,
        default=None,
        help="Optional node-label file (node label). Reorders nodes by label.",
    )
    p_snap.add_argument("--target-tiles", type=int, default=32, help="Approx tiles per dimension")

    p_long = sub.add_parser("longformer", help="Generate Longformer attention mask matrix")
    p_long.add_argument("--name", required=True, help="Output matrix/tile base name")
    p_long.add_argument("--seq-len", type=int, default=4096)
    p_long.add_argument("--window", type=int, default=512)
    p_long.add_argument("--global-tokens", type=int, default=32)
    p_long.add_argument("--target-tiles", type=int, default=32, help="Approx tiles per dimension")

    return p


def main() -> None:
    args = build_parser().parse_args()
    if args.cmd == "snap":
        ingest_snap_edge_list(
            input_path=args.input,
            output_name=args.name,
            undirected=args.undirected,
            target_tiles_per_dim=args.target_tiles,
            labels_path=args.labels,
        )
        return
    if args.cmd == "longformer":
        ingest_longformer_mask(
            output_name=args.name,
            seq_len=args.seq_len,
            window=args.window,
            global_tokens=args.global_tokens,
            target_tiles_per_dim=args.target_tiles,
        )
        return
    raise AssertionError("unreachable")


if __name__ == "__main__":
    main()

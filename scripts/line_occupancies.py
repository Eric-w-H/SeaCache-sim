import sys
import numpy as np

LINE_NBYTES         = 64
ELEM_DATA_NBYTES    = 8
ELEM_COORD_NBYTES   = 4 # assume each dimension is at most 2^32 (~4 billion) elements long

filepath    = sys.argv[1]
col_major   = '--col-major' in sys.argv
header_seen = False
row_dim = col_dim = None
rows, cols  = [], []

with open(filepath) as f:
    for line in f:
        if line.startswith('%'):
            continue
        if not header_seen:
            header_seen = True
            parts = line.split()
            row_dim, col_dim = int(parts[0]), int(parts[1])
            continue
        parts = line.split()
        rows.append(int(parts[0]))
        cols.append(int(parts[1]))

rows, cols = np.array(rows), np.array(cols)
if col_major:
    cl = (cols * row_dim + rows) * ELEM_DATA_NBYTES // LINE_NBYTES
else:
    cl = (rows * col_dim + cols) * ELEM_DATA_NBYTES // LINE_NBYTES
_, counts   = np.unique(cl, return_counts=True)
occ, freq   = np.unique(counts, return_counts=True)
total_lines = np.sum(freq) # nonempty lines only
breakeven   = LINE_NBYTES/(ELEM_DATA_NBYTES + ELEM_COORD_NBYTES)

print(f"Line size:          {LINE_NBYTES:4} bytes")
print(f"Elem data size:     {ELEM_DATA_NBYTES:4} bytes")
print(f"Elem coord size:    {ELEM_COORD_NBYTES:4} bytes")
print(f"Break-even occupancy for dense line: > {breakeven:.2f} elems per line\n")
for o, f in zip(occ, freq):
    print(f"occupancy {o}: {f/total_lines*100:6.2f}% ({f} lines)")

dense_lines = np.sum(freq[occ > breakeven])
print(f"Lines above break-even: {dense_lines}/{total_lines} ({dense_lines/total_lines*100:.2f}%)")
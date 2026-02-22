#!/usr/bin/env python3
import glob, os, sys
import numpy as np
from collections import defaultdict

# ----------------------------------------
# 1) Parse arguments
# ----------------------------------------
# python 03_countOverlapped.py <folder_path> <depth_threshold> <output_file>
folder_path, depth_thresh_str, output_file_path = sys.argv[1:]
depth_threshold = float(depth_thresh_str)

# ----------------------------------------
# 2) Find all .pos_dep.txt files
# ----------------------------------------
file_paths = glob.glob(os.path.join(folder_path, "*.pos_dep.txt"))
if not file_paths:
    raise ValueError(f"No .pos_dep.txt files found in {folder_path}")

# ----------------------------------------
# 3) Count lines per file (after depth filter)
# ----------------------------------------
file_row_counts = {}
files_filtered_lines = {}

for fp in file_paths:
    kept_lines = []
    with open(fp) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            cols = line.split()        # split on whitespace
            if len(cols) < 3:
                continue              # skip malformed rows
            depth = float(cols[2])
            if depth < depth_threshold:
                continue              # skip below‐threshold
            kept_lines.append(line)
    file_row_counts[fp] = len(kept_lines)
    files_filtered_lines[fp] = kept_lines

if not any(file_row_counts.values()):
    raise ValueError("No lines passed the depth filter in any file.")

# ----------------------------------------
# 4) Outlier filtering on file sizes
# ----------------------------------------
counts = np.array(list(file_row_counts.values()), dtype=float)
mean, std = counts.mean(), counts.std()
low, high = mean - 2*std, mean + 2*std

filtered_files = [
    fp for fp, cnt in file_row_counts.items()
    if low <= cnt <= high
]
if not filtered_files:
    raise ValueError("All files were filtered out as outliers.")

# ----------------------------------------
# 5) Build full‐row & prefix sets for each kept file
# ----------------------------------------
files_seq = {}
files_seq_pos = {}
for fp in filtered_files:
    rows = files_filtered_lines[fp]
    seq = set("{}".format(*line.split()[:1]) for line in rows)
    files_seq[fp] = seq
    seq_pos = set("{},{}".format(*line.split()[:2]) for line in rows)
    files_seq_pos[fp] = seq_pos

# ----------------------------------------
# 6) Count common items at various % thresholds
# ----------------------------------------
def collect_common_items(files_dict, thresholds):
    freq = defaultdict(int)
    for items in files_dict.values():
        for x in items:
            freq[x] += 1
    n = len(files_dict)
    counts = {}
    shared_items = {}
    for p in thresholds:
        need = int((p / 100) * n + 0.9999)
        selected = {x for x, cnt in freq.items() if cnt >= need}
        counts[p] = len(selected)
        shared_items[p] = selected
    return counts, shared_items

thresholds = [100, 90, 80, 70, 60, 50]
common_seq_counts, common_seq_sets = collect_common_items(files_seq, thresholds)
common_seq_pos_counts, common_seq_pos_sets = collect_common_items(files_seq_pos, thresholds)

# ----------------------------------------
# 7) Write the summary file (original output)
# ----------------------------------------
with open(output_file_path, "w") as out:
    out.write(f"Depth threshold: {depth_threshold}\n")
    out.write("Files after outlier filter:\n")
    for fp in filtered_files:
        out.write(os.path.basename(fp) + "\n")

    out.write("\nCommon Transcript:\n")
    for p in thresholds:
        out.write(f"{p}% {common_seq_counts[p]}\n")

    out.write("\nCommon Position (seq,pos):\n")
    for p in thresholds:
        out.write(f"{p}% {common_seq_pos_counts[p]}\n")

# ----------------------------------------
# 8) Write transcript IDs and transcript-position pairs to files
# ----------------------------------------
output_dir = os.path.dirname(output_file_path)
base_prefix = os.path.splitext(os.path.basename(output_file_path))[0]

for p in thresholds:
    out_seq = os.path.join(output_dir, f"{base_prefix}_shared_transcripts_{p}pct.txt")
    out_seq_pos = os.path.join(output_dir, f"{base_prefix}_shared_positions_{p}pct.txt")

    with open(out_seq, "w") as f_seq:
        for tid in sorted(common_seq_sets[p]):
            f_seq.write(tid + "\n")

    with open(out_seq_pos, "w") as f_pos:
        for tid_pos in sorted(common_seq_pos_sets[p]):
            f_pos.write(tid_pos + "\n")



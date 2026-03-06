#!/software/skylake/software/Python/3.10.8-GCCcore-12.2.0/bin/python
import os
import sys
import re
import pysam
from collections import OrderedDict

# === Defaults ===
DELIM_PATTERN = r"[_\-.]"   # split sample IDs on _, -, or .
GROUP_INDEX = 1             # "second part" of the sample ID
COVERAGE_CUTOFFS = [0.7, 0.8, 0.9]
COVERAGE_LABELS = {0.7: "70", 0.8: "80", 0.9: "90"}  # for output

# Four independent criteria
CRITERIA = OrderedDict({
    "DP>=20": lambda sdat: (sdat.get("GT") is not None and None not in sdat["GT"]
                            and sdat.get("DP") is not None and sdat["DP"] >= 20),
    "DP>=30": lambda sdat: (sdat.get("GT") is not None and None not in sdat["GT"]
                            and sdat.get("DP") is not None and sdat["DP"] >= 30),
    "DP>=40": lambda sdat: (sdat.get("GT") is not None and None not in sdat["GT"]
                            and sdat.get("DP") is not None and sdat["DP"] >= 40),
    "GQ>=20": lambda sdat: (sdat.get("GT") is not None and None not in sdat["GT"]
                            and sdat.get("GQ") is not None and sdat["GQ"] >= 20),
})

def prefix_candidates(base_name_first_token: str):

    cands = [base_name_first_token]
    # strip ploidy level
    stripped = re.sub(r'_[0-9]+x$', '', base_name_first_token)
    if stripped and stripped != base_name_first_token:
        cands.append(stripped)
    # de-duplicate while preserving order
    seen, out = set(), []
    for c in cands:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out

def read_skip_samples(replicates_path):

    skip = set()
    with open(replicates_path, "r") as fh:
        header = fh.readline()  # consume header
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\s+", line)
            if len(parts) >= 3:
                skip.add(parts[2])
    return skip

def read_diversity_map(diversity_path):

    d = {}
    with open(diversity_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\s+", line)
            if len(parts) >= 2:
                key = parts[0]
                try:
                    d[key] = int(parts[1])
                except ValueError:
                    # Keep it out if not an int
                    continue
    return d

def group_from_sample(sample, delim_pattern, group_index):
    toks = re.split(delim_pattern, sample)
    if len(toks) <= group_index:
        raise ValueError(
            f"Sample '{sample}' has no token at index {group_index} after split by '{delim_pattern}'. Tokens={toks}"
        )
    return toks[group_index]

def build_groups(sample_names, delim_pattern, group_index, skip_samples):

    # Filter out skipped samples
    kept_indices = [i for i, s in enumerate(sample_names) if s not in skip_samples]
    kept_samples = [sample_names[i] for i in kept_indices]

    groups = OrderedDict()
    for local_i, s in enumerate(kept_samples):
        g = group_from_sample(s, delim_pattern, group_index)
        groups.setdefault(g, []).append(local_i)

    return groups, kept_samples, kept_indices

def ensure_header(path):
    need_header = (not os.path.exists(path)) or os.path.getsize(path) == 0
    if need_header:
        with open(path, "w") as f:
            f.write("base_name,group_name,criteria,cutoff,Expected_H,diversity\n")

def process_vcf(vcf_path, replicates_path, diversity_path, output_file):
    base_name = os.path.basename(vcf_path)
    # base prefix for mapping diversity (substring before first dot)
    # base_prefix = base_name.split(".")[0]
    base_first_token = base_name.split(".")[0]
    prefixes = prefix_candidates(base_first_token)

    # Read auxiliary files
    skip_samples = read_skip_samples(replicates_path)
    diversity_map = read_diversity_map(diversity_path)

    vcf = pysam.VariantFile(vcf_path)

    # Build groups from samples (after skipping)
    original_samples = list(vcf.header.samples)
    groups, kept_samples, kept_indices = build_groups(original_samples, DELIM_PATTERN, GROUP_INDEX, skip_samples)

    # Map kept local indices to actual record sample names for fast access
    kept_sample_names = [original_samples[i] for i in kept_indices]

    # Accumulators: group -> criterion -> cutoff -> list of Expected H per site (only biallelic)
    het_values = {g: {crit: {c: [] for c in COVERAGE_CUTOFFS} for crit in CRITERIA} for g in groups}

    for rec in vcf.fetch():  # full VCF
        is_biallelic = rec.alts is not None and len(rec.alts) == 1

        # Preload per-sample data for kept samples only, in kept order
        sdat_list = [rec.samples[name] for name in kept_sample_names]

        # For each group & criterion, check coverage and (if biallelic) compute H
        for g, group_local_indices in groups.items():
            if not group_local_indices:
                continue

            for crit_name, crit_fn in CRITERIA.items():
                # Which samples in this group pass this criterion?
                passing_local = [li for li in group_local_indices if crit_fn(sdat_list[li])]
                if not passing_local:
                    continue

                valid_fraction = len(passing_local) / len(group_local_indices)

                for cutoff in COVERAGE_CUTOFFS:
                    if valid_fraction >= cutoff and is_biallelic:
                        # Count alleles using only passing samples
                        ref_count = 0
                        alt_count = 0
                        for li in passing_local:
                            gt = sdat_list[li].get("GT")
                            if gt is None or None in gt:
                                continue
                            ref_count += sum(1 for a in gt if a == 0)
                            alt_count += sum(1 for a in gt if a == 1)

                        total = ref_count + alt_count
                        if total > 0:
                            p = alt_count / total
                            het = 2 * p * (1 - p)
                            het_values[g][crit_name][cutoff].append(het)

    # Write CSV rows
    ensure_header(output_file)
    with open(output_file, "a") as f:
        for g in groups.keys():
            # Resolve diversity value for this base+group
            # key = f"{base_prefix}_{g}"
            # diversity_val = diversity_map.get(key, "NA")
            diversity_val = "NA"
            for pref in prefixes:
                k = f"{pref}_{g}"
                if k in diversity_map:
                    diversity_val = diversity_map[k]
                    break

            for crit_name in CRITERIA.keys():
                for cutoff in COVERAGE_CUTOFFS:
                    vals = het_values[g][crit_name][cutoff]
                    if len(vals) == 0:
                        f.write(f"{base_name},{g},{crit_name},{COVERAGE_LABELS[cutoff]},NA,{diversity_val}\n")
                    else:
                        avg = sum(vals) / len(vals)
                        f.write(f"{base_name},{g},{crit_name},{COVERAGE_LABELS[cutoff]},{avg:.4f},{diversity_val}\n")

# === Run ===
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <vcf_path> <fileA> <fileB> <output_file>")
        sys.exit(1)

    vcf_path = sys.argv[1]
    replicates_path = sys.argv[2]
    diversity_path = sys.argv[3]
    output_file = sys.argv[4]

    for pth, label in [(vcf_path, "VCF"), (replicates_path, "fileA"), (diversity_path, "fileB")]:
        if not os.path.exists(pth):
            print(f"Error: {label} not found: {pth}")
            sys.exit(1)

    process_vcf(vcf_path, replicates_path, diversity_path, output_file)


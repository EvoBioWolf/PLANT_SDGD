#!/software/skylake/software/Python/3.10.8-GCCcore-12.2.0/bin/python
import os
import sys
import re
import pysam
from collections import OrderedDict

# === Defaults ===
DELIM_PATTERN = r"[_\-.]"   # split sample IDs on _, -, or .
GROUP_INDEX = 1             # "second part" of the sample ID

# we now only use 80% coverage for Fst
COVERAGE_CUTOFF = 0.8
COVERAGE_LABEL = "80"

# Four independent criteria (same structure as your previous script)
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
    """
    Return possible prefixes for diversity lookup.
    1) the raw first token (e.g., 'VicCra_4x')
    2) the token with a trailing '_<digits>x' removed (e.g., 'VicCra')
    """
    cands = [base_name_first_token]
    stripped = re.sub(r'_[0-9]+x$', '', base_name_first_token)
    if stripped and stripped != base_name_first_token:
        cands.append(stripped)
    seen, out = set(), []
    for c in cands:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out

## exclude replicates
def read_skip_samples(replicates_path):
    """Read replicates files (with header). Return a set of sample IDs from the 3rd column to skip."""
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
    """Read plot species richness (no header). Return dict key -> int diversity."""
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
    """Return OrderedDict(group -> list of sample indices) and filtered sample list."""
    kept_indices = [i for i, s in enumerate(sample_names) if s not in skip_samples]
    kept_samples = [sample_names[i] for i in kept_indices]

    groups = OrderedDict() # group is plot
    for local_i, s in enumerate(kept_samples):
        g = group_from_sample(s, delim_pattern, group_index)
        groups.setdefault(g, []).append(local_i)

    return groups, kept_samples, kept_indices

def ensure_header(path):
    need_header = (not os.path.exists(path)) or os.path.getsize(path) == 0
    if need_header:
        with open(path, "w") as f:
            f.write(
                "base_name,group1,group2,criteria,cutoff,"
                "n_loci,H_T_sum,H_S_sum,Fst,group1_diversity,group2_diversity\n"
            )

def process_vcf(vcf_path, replicates_path, diversity_path, output_file):
    base_name = os.path.basename(vcf_path)
    base_first_token = base_name.split(".")[0]
    prefixes = prefix_candidates(base_first_token)

    # Read replicates and richness files
    skip_samples = read_skip_samples(replicates_path)
    diversity_map = read_diversity_map(diversity_path)

    vcf = pysam.VariantFile(vcf_path)

    # Build groups from samples (after skipping)
    original_samples = list(vcf.header.samples)
    groups, kept_samples, kept_indices = build_groups(
        original_samples, DELIM_PATTERN, GROUP_INDEX, skip_samples
    )

    # Names in the same order as local indices used in `groups`
    kept_sample_names = [original_samples[i] for i in kept_indices]

    # Map group -> diversity value (species richness)
    group_diversity = {}
    for g in groups.keys():
        div_val = "NA"
        for pref in prefixes:
            k = f"{pref}_{g}"
            if k in diversity_map:
                div_val = diversity_map[k]
                break
        group_diversity[g] = div_val

    # Prepare pairwise structures
    group_names = list(groups.keys())
    pair_keys = []
    for i in range(len(group_names)):
        for j in range(i + 1, len(group_names)):
            pair_keys.append((group_names[i], group_names[j]))

    # pair_stats[(g1,g2)][crit] = {"H_T_sum":..., "H_S_sum":..., "n_loci":...}
    pair_stats = {
        (g1, g2): {
            crit_name: {"H_T_sum": 0.0, "H_S_sum": 0.0, "n_loci": 0}
            for crit_name in CRITERIA.keys()
        }
        for (g1, g2) in pair_keys
    }

    # Iterate over variants
    for rec in vcf.fetch():
        # Only biallelic SNPs
        if rec.alts is None or len(rec.alts) != 1:
            continue

        # Preload per-sample data for kept samples only, in kept order
        sdat_list = [rec.samples[name] for name in kept_sample_names]

        # For each criterion, compute per-group allele counts + coverage
        for crit_name, crit_fn in CRITERIA.items():
            # group_info[g] = None (if group fails coverage or has no alleles)
            # or dict with ref, alt, n_alleles
            group_info = {}

            for g, group_local_indices in groups.items():
                group_size = len(group_local_indices)
                if group_size == 0:
                    group_info[g] = None
                    continue

                ref_count = 0
                alt_count = 0
                n_pass_geno = 0

                for li in group_local_indices:
                    sdat = sdat_list[li]
                    if not crit_fn(sdat):
                        continue

                    gt = sdat.get("GT")
                    if gt is None or None in gt:
                        continue

                    n_pass_geno += 1
                    for a in gt:
                        if a == 0:
                            ref_count += 1
                        elif a == 1:
                            alt_count += 1
                        # ignore other allele codes if present

                valid_fraction = n_pass_geno / group_size
                total_alleles = ref_count + alt_count

                if valid_fraction >= COVERAGE_CUTOFF and total_alleles > 0:
                    group_info[g] = {
                        "ref": ref_count,
                        "alt": alt_count,
                        "n_alleles": total_alleles,
                        "valid_fraction": valid_fraction,
                    }
                else:
                    group_info[g] = None

            # Now compute Nei-style components for each pair where both groups pass
            for (g1, g2) in pair_keys:
                info1 = group_info.get(g1)
                info2 = group_info.get(g2)
                if info1 is None or info2 is None:
                    continue

                nA = info1["n_alleles"]
                nB = info2["n_alleles"]
                if nA == 0 or nB == 0:
                    continue

                pA = info1["alt"] / nA
                pB = info2["alt"] / nB

                # Expected heterozygosity within each group
                HA = 2.0 * pA * (1.0 - pA)
                HB = 2.0 * pB * (1.0 - pB)

                # Within-population component (H_S) weighted by allele counts
                HS = (nA * HA + nB * HB) / (nA + nB)

                # Pooled allele frequency and total heterozygosity (H_T)
                p_bar = (nA * pA + nB * pB) / (nA + nB)
                HT = 2.0 * p_bar * (1.0 - p_bar)

                # Accumulate sums (we still count the locus even if HT==0)
                stats = pair_stats[(g1, g2)][crit_name]
                stats["H_T_sum"] += HT
                stats["H_S_sum"] += HS
                stats["n_loci"] += 1

    # Write CSV rows
    ensure_header(output_file)
    with open(output_file, "a") as f:
        for (g1, g2) in pair_keys:
            div1 = group_diversity.get(g1, "NA")
            div2 = group_diversity.get(g2, "NA")

            for crit_name in CRITERIA.keys():
                stats = pair_stats[(g1, g2)][crit_name]
                n_loci = stats["n_loci"]
                H_T_sum = stats["H_T_sum"]
                H_S_sum = stats["H_S_sum"]

                if n_loci == 0 or H_T_sum == 0.0:
                    fst_str = "NA"
                else:
                    fst = (H_T_sum - H_S_sum) / H_T_sum
                    fst_str = f"{fst:.6f}"

                # H_T_sum and H_S_sum can be useful to keep; format nicely
                if n_loci == 0:
                    H_T_str = "NA"
                    H_S_str = "NA"
                else:
                    H_T_str = f"{H_T_sum:.6f}"
                    H_S_str = f"{H_S_sum:.6f}"

                f.write(
                    f"{base_name},{g1},{g2},{crit_name},{COVERAGE_LABEL},"
                    f"{n_loci},{H_T_str},{H_S_str},{fst_str},{div1},{div2}\n"
                )

# === Run ===
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script_pairwise_fst.py <vcf_path> <fileA> <fileB> <output_file>")
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


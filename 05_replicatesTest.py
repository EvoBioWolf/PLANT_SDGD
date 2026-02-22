#!/software/skylake/software/Python/3.10.8-GCCcore-12.2.0/bin/python
import sys
import pysam
from collections import defaultdict

if len(sys.argv) < 4:
    sys.exit(f"Usage: {sys.argv[0]} pairs.txt input.vcf.gz output.tsv")

pairs_file = sys.argv[1]
vcf_file   = sys.argv[2]
out_file   = sys.argv[3]

# ----------------------------
# Load sample pairs
# ----------------------------
pairs = []
with open(pairs_file) as f:
    for line in f:
        if line.strip():
            if line.startswith('ID'):
                continue
            pid, s1, s2 = line.strip().split()
            pairs.append((pid, s1, s2))

# ID	Sample1	Sample2
# random1	28-H3-04	28-M2-05
# random2	28-H3-12	28-H3-04
# random3	28-H4-01	28-M1-09

# ----------------------------
# Define stats structure
# ----------------------------
def make_stats():
    return {
        "total_sites": 0,
        "identical_total_exact": 0,
        "identical_nonref_exact": 0,
        "identical_het_exact": 0,
        "different_exact": 0,
        "identical_total_collapsed": 0,
        "identical_nonref_collapsed": 0,
        "identical_het_collapsed": 0,
        "different_collapsed": 0,
    }

# Filtering criteria (independent)
criteria = {
    "DP20": lambda s1, s2: s1.get("DP", 0) >= 20 and s2.get("DP", 0) >= 20,
    "DP30": lambda s1, s2: s1.get("DP", 0) >= 30 and s2.get("DP", 0) >= 30,
    "DP40": lambda s1, s2: s1.get("DP", 0) >= 40 and s2.get("DP", 0) >= 40,
    "GQ20": lambda s1, s2: s1.get("GQ", 0) >= 20 and s2.get("GQ", 0) >= 20,
}

# Results storage
results = defaultdict(lambda: {crit: make_stats() for crit in criteria})

# ----------------------------
# Helper functions
# ----------------------------
def gt_to_str(gt):
    """Convert GT tuple to string '0011'."""
    return "".join(str(a) for a in sorted(gt) if a is not None)

def collapse_gt(gt_str):
    """Collapse ploidy to diploid-style 00, 01, 11."""
    if not gt_str:
        return None
    alleles = set(gt_str)
    if alleles == {"0"}:
        return "00"
    elif alleles == {"1"}:
        return "11"
    elif "0" in alleles and "1" in alleles:
        return "01"
    else:
        return None

# ----------------------------
# Process VCF
# ----------------------------
vcf = pysam.VariantFile(vcf_file)

for rec in vcf.fetch():
    # only bi-allelic SNPs
    if not rec.alts or len(rec.alts) != 1:
        continue

    for pid, s1, s2 in pairs:
        if s1 not in rec.samples or s2 not in rec.samples:
            continue

        d1, d2 = rec.samples[s1], rec.samples[s2]
        gt1, gt2 = d1.get("GT"), d2.get("GT")

        # skip missing genotypes
        if gt1 is None or gt2 is None:
            continue
        if any(a is None for a in gt1) or any(a is None for a in gt2):
            continue

        gt1_str, gt2_str = gt_to_str(gt1), gt_to_str(gt2)
        gt1_c, gt2_c     = collapse_gt(gt1_str), collapse_gt(gt2_str)

        for crit_name, crit_fn in criteria.items():
            if not crit_fn(d1, d2):
                continue

            stats = results[(pid, s1, s2)][crit_name]
            stats["total_sites"] += 1

            # exact ploidy
            if gt1_str == gt2_str:
                stats["identical_total_exact"] += 1
                if gt1_str in ["01", "11", "0001", "0011", "0111", "1111"]:
                    stats["identical_nonref_exact"] += 1     # most important
                if gt1_str in ["01", "0001", "0011", "0111"]:
                    stats["identical_het_exact"] += 1
            else:
                stats["different_exact"] += 1

            # collapsed diploid-style
            if gt1_c and gt2_c:    # 0/0/0/1 VS. 0/0/1/1
                if gt1_c == gt2_c:
                    stats["identical_total_collapsed"] += 1
                    if gt1_c in ["01", "11"]:
                        stats["identical_nonref_collapsed"] += 1
                    if gt1_c == "01":
                        stats["identical_het_collapsed"] += 1
                else:
                    stats["different_collapsed"] += 1

# ----------------------------
# Write output
# ----------------------------
with open(out_file, "w") as out:
    header = [
        "pairID", "sample1", "sample2", "criteria",
        "total_sites",
        "identical_total_exact", "identical_nonref_exact", "identical_het_exact", "different_exact",
        "identical_total_collapsed", "identical_nonref_collapsed", "identical_het_collapsed", "different_collapsed",
    ]
    out.write("\t".join(header) + "\n")

    for (pid, s1, s2), crit_dict in results.items():
        for crit_name, stats in crit_dict.items():
            row = [pid, s1, s2, crit_name] + [str(stats[h]) for h in header[4:]]
            out.write("\t".join(row) + "\n")


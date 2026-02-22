#! /software/skylake/software/Python/3.10.8-GCCcore-12.2.0/bin/python
import sys
import math
import re
import pysam
from collections import defaultdict

vcf_file, diversity_file, expr_file, filter_file, species_name, output_file = sys.argv[1:]

# ===== diversity data (species richness) ===== 
diversity_log2, diversity = {}, {}
with open(diversity_file) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) != 2:
            continue
        key, val = parts
        try:
            v = float(val)
            # protect against non-positive values
            if v > 0:
                diversity_log2[key] = math.log2(v)
                diversity[key] = v
            else:
                diversity_log2[key] = ""
                diversity[key] = ""
        except ValueError:
            diversity_log2[key] = ""
            diversity[key] = ""

# ===== Load expression data (from isoSeq data) ===== 
expr_counts = {} 
with open(expr_file) as f: 
    for line in f: 
        parts = line.strip().split() 
        if len(parts) != 2:
            continue 
        expr, gene_id = parts
        expr_counts[gene_id] = float(expr)

# ===== Load replicates (exclude replicate sequence) ===== 
replicate_samples = set()
with open(filter_file, "r") as f:
    header = True
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        cols = line.split()
        if header:
            header = False
            continue
        if len(cols) >= 3:
            replicate_samples.add(cols[2])

# ===== helper =====
def extract_plot_id(sample_id: str):
    """Split on '-' or '_' and return the second token (e.g., '28-H3-14' -> 'H3', '28-M2_13' -> 'M2')."""
    tokens = re.split(r"[-_]", sample_id)
    return tokens[1] if len(tokens) >= 2 else ""

def gt_to_str(gt):
    """Convert GT tuple to string"""
    if gt is None or any(a is None for a in gt):
        return None
    return "".join(map(str, sorted(gt)))

def is_heterozygous(gt_str):
    """Return True if genotype contains both 0 and 1"""
    alleles = set(gt_str)
    return "0" in alleles and "1" in alleles

def safe_ln(x):
    try:
        if x is None or x <= 0:
            return ""
        return math.log(x)
    except Exception:
        return ""


# Criteria definitions
CRITERIA = [
    ("GQ20", lambda gq, dp: (gq is not None and gq >= 20)),
    ("DP20", lambda gq, dp: (dp is not None and dp >= 20)),
    ("DP30", lambda gq, dp: (dp is not None and dp >= 30)),
    ("DP40", lambda gq, dp: (dp is not None and dp >= 40)),
]


# Open VCF
vcf = pysam.VariantFile(vcf_file)

# results[(chrom, sample, criterion)] = [genotyped_count, polymorphic_count]
results = defaultdict(lambda: [0, 0])

# Iterate through VCF
for rec in vcf.fetch():

    chrom = rec.chrom

    for sample, data in rec.samples.items():
        if sample in replicate_samples:
            continue

        gt = data.get("GT")
        gt_str = gt_to_str(gt)
        if gt_str is None:
            continue  # skip missing genotypes

        # Pull per-sample FORMAT fields
        gq = data.get("GQ")
        dp = data.get("DP")

        # Coerce missing to 0 (so the comparisons work)
        gq = 0 if gq is None else gq
        dp = 0 if dp is None else dp

        for crit_name, crit_func in CRITERIA:
            if crit_func(gq, dp):
                key = (chrom, sample, crit_name)
                # genotyped site: passed the criterion and GT is present
                results[key][0] += 1
                # polymorphic site: heterozygous 0/1
                if is_heterozygous(gt_str):
                    results[key][1] += 1

# Write output
with open(output_file, "w") as out:
    out.write("\t".join([
        "geneID", "sampleID", "genotyped_site", "polymorphic_site",
        "Expr", "plotID", "diversity", "diversity_log2", "genotyped_site_ln",
        "Expr_ln", "criteria"
    ]) + "\n")

    # stable order: by chrom, sample, our CRITERIA order
    crit_order = {name: i for i, (name, _) in enumerate(CRITERIA)}
    def sort_key(item):
        (chrom, sample, crit), _ = item
        return (chrom, sample, crit_order.get(crit, 999))

    for (chrom, sample, crit), (genotyped, polymorphic) in sorted(results.items(), key=sort_key):
        # join NumReads by geneID (chrom)
        expression = expr_counts.get(chrom, 1.0)

        # plotID and diversity_log2
        plot_id = extract_plot_id(sample)
        div_key = f"{species_name}_{plot_id}" if plot_id else ""
        div_log2 = diversity_log2.get(div_key, "")
        div = diversity.get(div_key, "")

        # logs
        genotyped_ln = safe_ln(genotyped)
        expression_ln = safe_ln(expression)

        out.write("\t".join(map(str, [
            chrom,
            sample,
            genotyped,
            polymorphic,
            expression,
            plot_id,
            div,
            div_log2 if div_log2 != "" else "",
            genotyped_ln if genotyped_ln != "" else "",
            expression_ln if expression_ln != "" else "No_expression",
            crit
        ])) + "\n")



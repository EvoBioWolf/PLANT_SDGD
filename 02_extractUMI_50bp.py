import gzip
from Bio import SeqIO
import sys

# Input and output gzipped files for R1 and R2
input_r1, input_r2, output_r1, output_r2 = sys.argv[1:]

# Open gzipped input and output files
with gzip.open(input_r1, "rt") as r1_in, gzip.open(input_r2, "rt") as r2_in, \
     gzip.open(output_r1, "wt") as r1_out, gzip.open(output_r2, "wt") as r2_out:

    # Iterate over paired-end reads
    for r1_record, r2_record in zip(SeqIO.parse(r1_in, "fastq"), SeqIO.parse(r2_in, "fastq")):
        # Trim R1 to first 16 bases
        trimmed_r1 = r1_record[:16]  # Slicing SeqRecord directly
        trimmed_r1.seq = r1_record.seq[:16]
        trimmed_r1.letter_annotations["phred_quality"] = r1_record.letter_annotations["phred_quality"][:16]

        # Trim R2 to first 50 bases
        trimmed_r2 = r2_record[:50]  # Slicing SeqRecord directly
        trimmed_r2.seq = r2_record.seq[:50]
        trimmed_r2.letter_annotations["phred_quality"] = r2_record.letter_annotations["phred_quality"][:50]

        # Write the modified reads to output
        SeqIO.write(trimmed_r1, r1_out, "fastq")
        SeqIO.write(trimmed_r2, r2_out, "fastq")

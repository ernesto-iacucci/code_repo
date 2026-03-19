# Language: Python
# RNA-Seq pipeline using Python
# This example integrates subprocess calls for command-line tools (FastQC, Trim Galore, HISAT2, featureCounts)
# and uses Python for processing, file management, and minimal analysis.

import os
import subprocess
import glob
import pandas as pd

# -------------------------
# Configuration
# -------------------------

# Paths to raw fastq files (paired-end)
RAW_DATA_DIR = "raw_data"
FASTQ_FILES = glob.glob(os.path.join(RAW_DATA_DIR, "*.fastq.gz"))

# Output directories
QC_DIR = "qc_results"
TRIMMED_DIR = "trimmed_data"
ALIGN_DIR = "alignment"
COUNT_DIR = "counts"

# Reference genome files (replace with actual paths)
REFERENCE_GENOME = "reference/genome.fa"
GTF_FILE = "reference/annotation.gtf"

# Ensure output directories exist
for folder in [QC_DIR, TRIMMED_DIR, ALIGN_DIR, COUNT_DIR]:
    os.makedirs(folder, exist_ok=True)

# -------------------------
# Step 1: Quality Control with FastQC
# -------------------------
for fq in FASTQ_FILES:
    # Run FastQC on each FASTQ file
    subprocess.run(["fastqc", fq, "-o", QC_DIR])

# -------------------------
# Step 2: Trimming adapters using Trim Galore
# -------------------------
for fq in FASTQ_FILES:
    # Create output filename
    out_fq = os.path.join(TRIMMED_DIR, os.path.basename(fq))
    # Run Trim Galore
    subprocess.run(["trim_galore", "--paired", "--output_dir", TRIMMED_DIR, fq])

# Gather trimmed files for alignment (assuming paired-end: *_val_1.fq.gz and *_val_2.fq.gz)
TRIMMED_R1 = sorted(glob.glob(os.path.join(TRIMMED_DIR, "*_val_1.fq.gz")))
TRIMMED_R2 = sorted(glob.glob(os.path.join(TRIMMED_DIR, "*_val_2.fq.gz")))

# -------------------------
# Step 3: Alignment using HISAT2
# -------------------------
for r1, r2 in zip(TRIMMED_R1, TRIMMED_R2):
    sample_name = os.path.basename(r1).split("_val_1")[0]
    # Output BAM file path
    bam_file = os.path.join(ALIGN_DIR, f"{sample_name}.bam")
    # HISAT2 alignment and conversion to BAM
    hisat2_cmd = f"hisat2 -x {REFERENCE_GENOME} -1 {r1} -2 {r2} | samtools view -bS - > {bam_file}"
    subprocess.run(hisat2_cmd, shell=True)

# -------------------------
# Step 4: Count reads with featureCounts
# -------------------------
bam_files = glob.glob(os.path.join(ALIGN_DIR, "*.bam"))
# featureCounts command to count reads per gene
featurecounts_cmd = ["featureCounts", "-a", GTF_FILE, "-o", os.path.join(COUNT_DIR, "gene_counts.txt")] + bam_files
subprocess.run(featurecounts_cmd)

# -------------------------
# Step 5: Load counts and basic analysis with Python (pandas)
# -------------------------
counts_file = os.path.join(COUNT_DIR, "gene_counts.txt")
counts_df = pd.read_csv(counts_file, sep="\t", comment="#", index_col=0)
# Keep only count columns after first annotation columns (like Chr, Start, End)
counts_only = counts_df.iloc[:, 5:]
print("Top 5 genes by counts:")
print(counts_only.sum(axis=1).sort_values(ascending=False).head())

# -------------------------
# Note / Next Steps
# -------------------------
# - For differential expression, consider using DESeq2 or edgeR in R by exporting counts_only.
# - This pipeline assumes paired-end fastq files and that tools like FastQC, Trim Galore,
#   HISAT2, samtools, and featureCounts are installed and in PATH.
# - Extensive logging and multi-threading can be added for production use.
# - Python can also use subprocess.Popen with better error handling for robustness.
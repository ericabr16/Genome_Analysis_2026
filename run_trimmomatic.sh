#!/bin/bash
#SBATCH -A uppmax2026-1-61
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J trim_qc_embr16
#SBATCH --output=trim_qc_%j.out
#SBATCH --error=trim_qc_%j.err

# Define environment directories
RAW_DIR="/proj/uppmax2026-1-61/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/raw"
OUT_DIR="/proj/uppmax2026-1-61/nobackup/work/embr16/transcriptomics"

cd $OUT_DIR

# Load specific Pelle modules
module load Trimmomatic/0.39-Java-17
module load FastQC/0.12.1-Java-17

# Run Trimmomatic PE
trimmomatic PE -phred33 \
  $RAW_DIR/ERR1797972_1.fastq.gz \
  $RAW_DIR/ERR1797972_2.fastq.gz \
  ERR1797972_trimmed_R1_paired.fq.gz ERR1797972_trimmed_R1_unpaired.fq.gz \
  ERR1797972_trimmed_R2_paired.fq.gz ERR1797972_trimmed_R2_unpaired.fq.gz \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Run FastQC on the outputs
fastqc ERR1797972_trimmed_R1_paired.fq.gz ERR1797972_trimmed_R2_paired.fq.gz

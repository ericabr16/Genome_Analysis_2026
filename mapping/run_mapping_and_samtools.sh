#!/bin/bash
#SBATCH -A uppmax2026-1-61
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -J mapping_and_samtools
#SBATCH -M pelle

# Load required bioinformatics modules on UPPMAX
module load bioinfo-tools
module load hisat2
module load samtools
module load htseq

# Define workspace directories
REFERENCE="/proj/uppmax2026-1-61/nobackup/work/embr16/assembly/output/E745_assembly.contigs.fasta"
INDEX_PREFIX="E745_index"
DATA_DIR="/proj/uppmax2026-1-61/nobackup/work/embr16/transcriptomics"
OUTPUT_DIR="/proj/uppmax2026-1-61/nobackup/work/embr16/transcriptomics/mapping"

cd $OUTPUT_DIR

# Step 1: Build the Hisat2 reference genome index
hisat2-build $REFERENCE $INDEX_PREFIX

# List of sample run IDs matching your biological replicates
SAMPLES=("ERR1797972" "ERR1797973" "ERR1797974" "ERR1797969" "ERR1797970" "ERR1797971")

# Step 2: Run loop for Alignment, SAMtools conversion/sorting/indexing, and HTSeq counting
for SAMPLE in "${SAMPLES[@]}"
do
    echo "Processing sample: ${SAMPLE}"

    # 2a. Align paired-end reads to indexed reference
    hisat2 -x $INDEX_PREFIX \
           -1 ${DATA_DIR}/${SAMPLE}_trimmed_R1_paired.fastq.gz \
           -2 ${DATA_DIR}/${SAMPLE}_trimmed_R2_paired.fastq.gz \
           -S ${SAMPLE}.sam --threads 4

    # 2b. SAMtools: Convert raw SAM to compressed binary BAM
    samtools view -bS ${SAMPLE}.sam > ${SAMPLE}.bam

    # 2c. SAMtools: Sort alignment by genomic coordinates
    samtools sort ${SAMPLE}.bam -o ${SAMPLE}.sorted.bam

    # 2d. SAMtools: Index the sorted BAM (generates the .bai file)
    samtools index ${SAMPLE}.sorted.bam

    # 2e. SAMtools QC: Calculate mapping statistics (used for your Wiki table!)
    samtools flagstat ${SAMPLE}.sorted.bam > ${SAMPLE}_flagstat.txt

    # 2f. Quantification: Run HTSeq-count to generate final expression tables
    htseq-count -f bam -r pos -s no -t CDS -i ID \
                ${SAMPLE}.sorted.bam ${OUTPUT_DIR}/E745_cleaned.gff > ${OUTPUT_DIR}/counts/${SAMPLE}.counts

    # Clean up massive raw intermediate alignment files to save storage space
    rm ${SAMPLE}.sam ${SAMPLE}.bam
done

echo "Mapping pipeline and post-mapping SAMtools processing completed successfully."

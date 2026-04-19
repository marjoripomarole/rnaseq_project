#!/bin/bash
# scripts/01_pipeline.sh
set -euo pipefail

SAMPLES="SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517 SRR1039520 SRR1039521"
THREADS=8
STAR_IDX="data/reference/star_index"
GTF="data/reference/Homo_sapiens.GRCh38.110.gtf"

for sample in ${SAMPLES}; do
    echo "========== Processing ${sample} =========="

    # QC
    fastqc data/raw/${sample}_1.fastq.gz data/raw/${sample}_2.fastq.gz \
        -o results/qc/ -t 2

    # Trim
    fastp -i data/raw/${sample}_1.fastq.gz -I data/raw/${sample}_2.fastq.gz \
          -o results/trimmed/${sample}_R1.fastq.gz -O results/trimmed/${sample}_R2.fastq.gz \
          --qualified_quality_phred 20 --length_required 50 \
          --html results/trimmed/${sample}_fastp.html \
          --json results/trimmed/${sample}_fastp.json --thread 4

    # Align with STAR (using FIFOs to stream compressed input on macOS)
    R1_FIFO=$(mktemp -u /tmp/${sample}_R1.XXXXXX)
    R2_FIFO=$(mktemp -u /tmp/${sample}_R2.XXXXXX)
    mkfifo "${R1_FIFO}" "${R2_FIFO}"
    gunzip -c results/trimmed/${sample}_R1.fastq.gz > "${R1_FIFO}" &
    gunzip -c results/trimmed/${sample}_R2.fastq.gz > "${R2_FIFO}" &
    STAR --genomeDir ${STAR_IDX} \
         --readFilesIn "${R1_FIFO}" "${R2_FIFO}" \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix results/aligned/${sample}_ \
         --runThreadN ${THREADS}
    rm -f "${R1_FIFO}" "${R2_FIFO}"

    # Index BAM
    samtools index results/aligned/${sample}_Aligned.sortedByCoord.out.bam
done

# Count reads per gene
featureCounts \
    -a ${GTF} \
    -o results/counts/gene_counts.txt \
    -T ${THREADS} \
    -p --countReadPairs \
    -s 2 \
    results/aligned/*_Aligned.sortedByCoord.out.bam

# Aggregate QC reports
multiqc results/ -o results/ --force

echo "Pipeline complete! Check results/multiqc_report.html"

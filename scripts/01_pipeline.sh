#!/bin/bash
# scripts/01_pipeline.sh
set -euo pipefail

SAMPLES="SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517 SRR1039520 SRR1039521"
THREADS=8
HISAT2_IDX="data/reference/hisat2_index/GRCh38"
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

    # Align with HISAT2 (STAR 2.7.11b has a known ARM64 macOS read-input bug)
    hisat2 -p ${THREADS} \
           -x ${HISAT2_IDX} \
           -1 results/trimmed/${sample}_R1.fastq.gz \
           -2 results/trimmed/${sample}_R2.fastq.gz \
           --dta \
           2>results/aligned/${sample}_hisat2.log \
        | samtools sort -@ ${THREADS} -o results/aligned/${sample}.bam

    # Index BAM
    samtools index results/aligned/${sample}.bam
done

# Count reads per gene
featureCounts \
    -a ${GTF} \
    -o results/counts/gene_counts.txt \
    -T ${THREADS} \
    -p --countReadPairs \
    -s 2 \
    results/aligned/*.bam

# Aggregate QC reports
multiqc results/ -o results/ --force

echo "Pipeline complete! Check results/multiqc_report.html"

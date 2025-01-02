# ChIP-seq data processing 워크플로우

## Trimming

trimming.sh

```bash
#!/bin/bash

# Set environment variables
export bdir="/data-blog/bi5/res"
export hg38_bowtieidx="/data-blog/bi5/hg38_bowtie_idx/hg38"
export hg38_bwaidx="/data-blog/bi5/hg38_bwa_index/hg38.fa"
export ebv_bowtie2idx="/data-blog/bi5/EBV_bowtie2_idx/NC_007605.1"
export ebv_bwaidx="/data-blog/bi5/EBV_bwa_index/NC_007605.1.fa"

# Set working directory
cd /data-blog/bi5/samples

# Create output directory
mkdir -p trimmed

# Set sample list and paths
samplist=("Input" "p65" "RIgG")
trimmomatic_jar="/data/packages/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapter_path="/data/packages/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# Trimming loop
for sampname in "${samplist[@]}"; do
    mkdir -p "trimmed/${sampname}"
    java -jar $trimmomatic_jar PE -threads 40 -trimlog "trimmed/${sampname}/log.txt" \
        "${sampname}_1.fastq" "${sampname}_2.fastq" \
        "trimmed/${sampname}/${sampname}_1.trimmed.fastq" "trimmed/${sampname}/${sampname}_1.unpaired.fastq" \
        "trimmed/${sampname}/${sampname}_2.trimmed.fastq" "trimmed/${sampname}/${sampname}_2.unpaired.fastq" \
        ILLUMINACLIP:"${adapter_path}":2:30:10 \
        LEADING:3 TRAILING:3 MINLEN:36
done

```

## Align

align.sh

```bash
#!/bin/bash

# Set environment variables
export bdir="/data-blog/bi5/res"
export hg38_bowtieidx="/data-blog/bi5/hg38_bowtie_idx/hg38"
export hg38_bwaidx="/data-blog/bi5/hg38_bwa_index/hg38.fa"
export ebv_bowtie2idx="/data-blog/bi5/EBV_bowtie2_idx/NC_007605.1"
export ebv_bwaidx="/data-blog/bi5/EBV_bwa_index/NC_007605.1.fa"

# Set working directory
cd /data-blog/bi5/samples

# Create directories for alignment results
mkdir -p $bdir/aln/bwa
mkdir -p $bdir/aln/bwt2

# Set sample list
samplist=("Input" "p65" "RIgG")

# EBV alignment loop
for sampname in "${samplist[@]}"; do
    bwa mem -t 20 -v 2 $ebv_bwaidx \
        "trimmed/${sampname}/${sampname}_1.trimmed.fastq" "trimmed/${sampname}/${sampname}_2.trimmed.fastq" \
        > "$bdir/aln/bwa/${sampname}_PE.bwa.sam"
    samtools view -hf 2 "$bdir/aln/bwa/${sampname}_PE.bwa.sam" | \
        samtools sort -o "$bdir/aln/bwa/${sampname}_PE.bwa.bam" -O BAM -@ 20 -
    samtools index -@ 20 "$bdir/aln/bwa/${sampname}_PE.bwa.bam"
    bamCoverage -b "$bdir/aln/bwa/${sampname}_PE.bwa.bam" \
        -o "$bdir/aln/bwa/${sampname}_PE.bwa.bam.bigwig"
    macs2 callpeak -t "$bdir/aln/bwa/${sampname}_PE.bwa.bam" -f BAMPE \
        -n "${sampname}" --outdir "$bdir/peaks_ebv"
done

# hg38 alignment loop
for sampname in "${samplist[@]}"; do
    bowtie2 -k1 --no-unal -p 40 --qc-filter -x $hg38_bowtieidx \
        -1 "trimmed/${sampname}/${sampname}_1.trimmed.fastq" -2 "trimmed/${sampname}/${sampname}_2.trimmed.fastq" \
        > "$bdir/aln/bwt2/${sampname}_PE_hg38.sam"
    samtools view -hf 2 "$bdir/aln/bwt2/${sampname}_PE_hg38.sam" | \
        samtools sort -o "$bdir/aln/bwt2/${sampname}_PE_hg38.bam" -O BAM -@ 20 -
    sambamba view -h -t 20 -f bam -p -F "[XS] == null and not unmapped and not duplicate" \
        "$bdir/aln/bwt2/${sampname}_PE_hg38.bam" > "$bdir/aln/bwt2/${sampname}_PE_hg38_u.bam"
    samtools index -@ 20 "$bdir/aln/bwt2/${sampname}_PE_hg38_u.bam"
    macs2 callpeak -t "$bdir/aln/bwt2/${sampname}_PE_hg38_u.bam" -f BAMPE \
        -n "${sampname}" --outdir "$bdir/peaks"
done
```

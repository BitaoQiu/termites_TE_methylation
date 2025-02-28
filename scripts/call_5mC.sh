#!/bin/sh
# Required files: 
# *_5mc.bam (HiFi reads with 5mC information)
# Genome assembly with minimap index. e.g., Mbel.mmi
# Genome assembly. e.g., Mbel.fasta

# Step 1. Map HiFi reads (with 5mC information) to the genome assembly.
genome_prefix=Mbel
genome=${genome_prefix}.mmi
mc_file=${genome_prefix}_5mc.fofn
ls *_5mc.bam > ${mc_file}
output=${genome_prefix}.5mc_mapped.bam
pbmm2 index $genome_prefix.fasta $genome
pbmm2 align $genome $mc_file --log-level INFO --sort -j 64 $output

# Step 2. Filter mapped reads with clipped aligment (See: https://github.com/tseemann/samclip).
genome=${genome_prefix}.fasta
bam_input=${genome_prefix}.5mc_mapped.bam
bam_output=${bam_input%.bam}.filtered.bam
samtools faidx ${genome_prefix}.fasta
samtools view -@ 64 -h $bam_input |  ~/bin/samclip-0.4.0/samclip --ref $genome --max 100 | samtools sort -@ 64 > $bam_output
samtools index -@ 64 $bam_output

# Step 3: Call genomic CpG 5mC information using pb-CpG-tools (see: https://github.com/PacificBiosciences/pb-CpG-tools)
out_prefix=${bam_output%.bam}.model.pbmm2
model=~/bin/pb-CpG-tools/models/pileup_calling_model.v1.tflite
~/bin/pb-CpG-tools/bin/aligned_bam_to_cpg_scores --bam $bam_output --output-prefix $out_prefix --threads 64 --model $model

###########################
##    WILDCARDS SET-UP   ##
###########################

import os
import glob
import subprocess
import pandas as pd
import yaml
import re
import itertools
from pathlib import Path


#### Functions to set up samples list in config from different sequencing batches

# Rule to update the configuration file
rule update_config:
  output:
    "snakeprofile/config.yaml"
  shell:
    "python scripts/extract_lcWGS_samples.py"

# Function to get samples from the config file
def get_samples(batch):
    with open("snakeprofile/config.yaml", 'r') as file:
        config = yaml.safe_load(file)
    return config['batches'][batch]

# Function to get batches from the config file
def get_batches():
    with open("snakeprofile/config.yaml", 'r') as file:
        config = yaml.safe_load(file)
    return list(config['batches'].keys())

# Ensure that the config is updated before any other rule runs
# Get all batches from the config
batches = get_batches()

# Get samples for each batch
samples_batch_1 = get_samples('batch_1')
samples_batch_2 = get_samples('batch_2')

# Try using this for expand in rule all: 
  # ("",
  # batch=lambda wildcards: get_batches(),
  # sample=lambda wildcards, batch: get_samples(batch))

##### Functions to set up lists of sample prefixes for merging .BAM files from different sequencing batches

# Extract the prefix of a sample name (e.g LCTGP-2)
def extract_sample_prefix(sample_name):
    match = re.match(r"^(.*?-\d+)", sample_name)
    if match:
        return match.group(1)
    return sample_name

# Find replicates across BAM files as inputs into merge_replicates rule
def find_replicates(sample_prefix, hap):
    pattern = f"results/bam_raw/hap{hap}/batch*/{sample_prefix}_*_hap{hap}.bam"
    files = list(Path().glob(pattern))
    
    # Debugging: Write the matched files to a debug file
    with open(f"debug/debug_find_replicates_{sample_prefix}_hap{hap}.txt", "w") as f:
        f.write(f"Pattern: {pattern}\n")
        if files:
            f.write("Found files:\n")
            for file in files:
                f.write(f"{file}\n")
        else:
            f.write("No files found.\n")
    
    return files

# Get all unique sample prefixes by extracting from samples listed in config file
def list_sample_prefixes():
    samples_batch_1 = get_samples('batch_1')
    samples_batch_2 = get_samples('batch_2')
    all_samples = set(samples_batch_1 + samples_batch_2)
    sample_prefixes = {extract_sample_prefix(sample) for sample in all_samples}
    return list(sample_prefixes)


# List of haplotypes
haps = ["1", "2"]

# List of sample prefixes
sample_prefixes = list_sample_prefixes()

# Debugging: Write sample_prefixes and haps to files
with open("debug_sample_prefixes.txt", "w") as f:
    f.write("Sample Prefixes:\n")
    for sp in sample_prefixes:
        f.write(f"{sp}\n")

with open("debug_haps.txt", "w") as f:
    f.write("Haplotypes:\n")
    for hap in haps:
        f.write(f"{hap}\n")

# List of populations
POPULATIONS = ("HPW", "IDNP-MW", "LCTGP", "MFNP", "PPP", "RLPLV", "SWCP")


#### Functions to create scaffold names list per scaffold
# Get first 24 scaffold names to later estimate LD per scaffold
# Not using other scaffolds as they don't quite align to the expected 24 chromosomes

rule get_scaffold_names:
  input:
    "data/reference/lupinehap1.fasta",
    "data/reference/lupinehap2.fasta"
  output:
    "results/scaffolds/hap1_scaffolds.txt",
    "results/scaffolds/hap2_scaffolds.txt"
  shell:
    """
    python scripts/extract_24_scaffold_names_by_hap.py
    """

with open("results/scaffolds/hap1_scaffolds.txt", "r") as file:
    HAP1SCAFFOLDS = [line.strip() for line in file.readlines()]

with open("results/scaffolds/hap2_scaffolds.txt", "r") as file:
    HAP2SCAFFOLDS = [line.strip() for line in file.readlines()]

# Split scaffold names by comma and create a list
HAP1SCAFFOLDS = [name.strip() for name in HAP1SCAFFOLDS]
HAP2SCAFFOLDS = [name.strip() for name in HAP2SCAFFOLDS]
HAP1SCAFFOLD_PREFIXES = [s.split("__")[0] for s in HAP1SCAFFOLDS]
HAP2SCAFFOLD_PREFIXES = [s.split("__")[0] for s in HAP2SCAFFOLDS]

# Function give full scaffold name given prefix of scaffold
def map_prefix_to_full_scaffold(prefix, hap_type):
    scaffold_list = HAP1SCAFFOLDS if hap_type == 1 else HAP2SCAFFOLDS
    for scaffold in scaffold_list:
        if scaffold.startswith(prefix):
            return scaffold
    return None  
# Return None or an appropriate default if not found

# Varying levels of Bonferroni-Hotchberg correction for ngsParalog
BH_VARS=[50,40,30,20,10,5]

# Varying levels of site depth to check filter levels
depth_params = [(50, 1500), (50, 2000), (50, 2500), (50, 3000), (100, 1500), (100, 2000), (100, 2500), (100, 3000)]

# Varying levels of minor allele frequency cutoffs
minMAF_params = [0.01, 0.001, 0.0001, 0.00001]

# For identifying populatio pairs for Fst analysis
POP_COMBINATIONS = list(itertools.combinations(POPULATIONS, 2))

###########################
## BEGINNING OF WORKFLOW ##
###########################


# Expand the final files using the updated configuration
rule all:
  input:
    expand("results/realSFS/hap2/fst/{pop1}_{pop2}_fst_global.txt", pop1=[x[0] for x in POP_COMBINATIONS], pop2=[x[1] for x in POP_COMBINATIONS])
    #expand("results/theta/hap2/{population}_log_scale.out", population=POPULATIONS)
    #expand("results/ngsRelate/hap2/{population}_inbreeding.out", population=POPULATIONS),
    #"results/plots/hap2/PCAngsd/all_popln_canonical_SNP_PCAngsd.png"
    #expand("results/ngsF/hap2/by_popln/{population}_ngsF-HMM_inbreeding.indF", population=POPULATIONS)
    #expand("results/theta/hap2/{population}_out.thetasWindow.gz.pestPG", population=POPULATIONS)
    #expand("results/ngsRelate/hap2/{population}_clone_summary.txt", population=POPULATIONS),
    #expand("results/theta/hap2/{population}_out.thetasWindow.gz.pestPG", population=POPULATIONS)
    #expand("results/plots/hap2/SFS/{population}_globalSFS_control_folded_minMAF{min_MAF}.png",
           #population=POPULATIONS, min_MAF = minMAF_params),
    #expand("results/plots/hap2/HWE/{population}_hwe_raw_SNPs_control_min{min_depth}_max{max_depth}.png",
           #population=POPULATIONS, min_depth=[50, 100], max_depth=[3050, 3100]),
    #expand("results/bed/hap{hap}/deviant_SNPs/{population}_devian.pestPGt_SNPs_calcLR_BH_corrected.BED",
            #hap=haps, population=POPULATIONS),
    #expand("results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM.BED",
            #hap=haps, population=POPULATIONS)
    

# Trim adapter ends off each sequence file using Trimmomatic
rule trim_reads:
  input:
    r1="data/{batch}/{sample}_L001_R1_001.fastq.gz",
    r2="data/{batch}/{sample}_L001_R2_001.fastq.gz"
  output:
    r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    r1_unp="results/fastq_trimmed/{batch}/{sample}_R1_unpaired.fastq.gz",
    r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz",
    r2_unp="results/fastq_trimmed/{batch}/{sample}_R2_unpaired.fastq.gz"
  log:
    "results/logs/trim_reads/{batch}/{sample}.log"
  envmodules:
    "trimmomatic/0.39"
  params:
    adapters="$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa"
  shell:
    "java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE {input.r1} {input.r2} "
    "{output.r1} {output.r1_unp} {output.r2} {output.r2_unp} "
    "ILLUMINACLIP:{params.adapters}:2:30:10:True "
    "LEADING:3 "
    "TRAILING:3 "
    "SLIDINGWINDOW:4:15 " # removes low quality bases
    "MINLEN:36 2> {log}"
  
# Unzip trimmed files as fastqc 0.12.0 cannot properly read compressed files. 
# Tried with 0.12.1 and it works!
rule unzip_files:
  input:
    zipped_r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    zipped_r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output: 
    unzipped_r1="results/fastq_trimmed_unzip/{batch}/{sample}_R1.fastq",
    unzipped_r2="results/fastq_trimmed_unzip/{batch}/{sample}_R2.fastq"
  shell:
    "gunzip -c {input.zipped_r1} > {output.unzipped_r1} && gunzip -c {input.zipped_r2} > {output.unzipped_r2}"

# Run FastQC per each trimmed sequence file
# Attempting zipped files since using fastqc/0.12.1
rule fastqc:
  input:
    fastq_r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    fastq_r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output:
    html_report_r1="results/fastqc/{batch}/{sample}_R1_fastqc.html",
    zip_report_r1="results/fastqc/{batch}/{sample}_R1_fastqc.zip",
    html_report_r2="results/fastqc/{batch}/{sample}_R2_fastqc.html",
    zip_report_r2="results/fastqc/{batch}/{sample}_R2_fastqc.zip"
  log:
    path="results/logs/fastQC/{batch}/{sample}.log"
  envmodules:
    "fastqc/0.12.1"
  shell:
    "fastqc {input.fastq_r1} {input.fastq_r2} --outdir results/fastqc/{wildcards.batch} 2> {log.path}"

# Create an aggregated FastQC report using MultiQC.
# Note that we create separate MultiQC reports for both batch 1 and 2
rule multiqc_raw:
  input:
    fastqc_dir="results/fastqc/{batch}"
  output:
    html_report="results/multiqc/multiqc_report_trimmed_{batch}.html"
  log:
    path="results/logs/multiqc/{batch}.log"
  params:
    fastqc_dir="results/fastqc/{batch}"
  shell:
    "multiqc -n {output.html_report} {input.fastqc_dir} 2> {log.path}"

# Creating faidx index for reference genome
rule faidx_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    "data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  log:
    "results/logs/refgen/lupinehap{hap}_faidx.log"
  envmodules:
    "samtools/1.20"
  shell:
    "samtools faidx {input} 2> {log}"

# Rules for indexing reference genomes (haplotypes 1 and 2)
rule index_reference:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  log:
    "results/logs/refgen/lupinehap{hap}_bwa_index.log"
  envmodules:
    "bwa/0.7.18"
  shell:
    "bwa index {input} 2> {log}"

# Rules for creating dictionary files
rule create_dict:
  input:
    "data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    "data/reference/hap{hap}/lupinehap{hap}.dict"
  log:
    "results/logs/refgen/hap{hap}_dict.log"
  envmodules:
    "samtools/1.20"
  shell:
    "samtools dict {input} > {output} 2> {log}"

# Removing reads smaller than 70bp so 'bwa mem' works better
# Note for a read to be kept, it must be greater than 70bp in BOTH read 1 and 2.
# There are some tools available in envmodules 'seqkit', but I found none that were compatible
rule filter_short_reads_with_seqkit:
  input:
    r1="results/fastq_trimmed/{batch}/{sample}_R1.fastq.gz",
    r2="results/fastq_trimmed/{batch}/{sample}_R2.fastq.gz"
  output:
    r1_filtered="results/fastq_filtered/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_filtered="results/fastq_filtered/{batch}/{sample}_R2_filtered.fastq.gz"
  envmodules:
    "StdEnv/2020",
    "seqkit/2.3.1"
  log:
    "results/logs/filter_fastq_70bp/{batch}/{sample}_filter_fastq.log"
  shell:
    """
    seqkit seq -m 70 {input.r1} | gzip > {output.r1_filtered}
    seqkit seq -m 70 {input.r2} | gzip > {output.r2_filtered}
    """


# Pair filtered reads back together
# Note: naming convention stays the same! Just different folder
rule pair_filtered_reads:
  input:
    r1_filtered="results/fastq_filtered/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_filtered="results/fastq_filtered/{batch}/{sample}_R2_filtered.fastq.gz"
  output:
    r1_paired="results/fastq_paired/{batch}/{sample}_R1_filtered.fastq.gz",
    r2_paired="results/fastq_paired/{batch}/{sample}_R2_filtered.fastq.gz"
  envmodules:
    "StdEnv/2020",
    "seqkit/2.3.1"
  log:
    "results/logs/pair_fastq_70bp/{batch}/{sample}_pair_fastq.log"
  shell:
    """
    seqkit pair \
    -1 {input.r1_filtered} \
    -2 {input.r2_filtered} \
    -O results/fastq_paired/{wildcards.batch}
    """

# NOTE: Prior to mapping, some people like to merge fastqs from the same individual/library and remove PCR duplicates prior to mapping using SuperDeduper from HTStream (last used 2020)
# This might be helpful in reducing heterozygote excess
# I have opted not to do this yet as it may be outdated.

# Mapping/Aligning reads to reference haplotypes
rule map_reads:
  input:
    r1="results/fastq_paired/{batch}/{sample}_R1_filtered.fastq.gz",
    r2="results/fastq_paired/{batch}/{sample}_R2_filtered.fastq.gz",
    genome="data/reference/hap{hap}/lupinehap{hap}.fasta",
    idx=multiext("data/reference/hap{hap}/lupinehap{hap}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
  output:
    "results/bam_raw/hap{hap}/{batch}/{sample}_hap{hap}.bam"
  log:
    "results/logs/map_reads/hap{hap}/{batch}/{sample}_hap{hap}.log"
  envmodules:
    "bwa/0.7.18",
    "samtools/1.20"
  threads: 8 
  params:
    RG="-R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA' "
  shell:
    """
    (bwa mem {params.RG} -t {threads} {input.genome} {input.r1} {input.r2} |\
    samtools view -u |\
    samtools sort - > {output}) 2> {log}
    """

##############################
#### DATA QUALITY CONTROL ####
##############################

#### DATA QUALITY CONTROL ####

# Steps in QC :
  # 1. Mark PCR and optical duplicates
  # 2. Produce HWE and SFS plots to evaluate mapping problems
  # 3. Identify paralogous regions causing mapping problems
      # - Requires ngsParalog and 1st run of ANGSD
      # - ANGSD: SNP call and input SNPs (aka polymorphic sites) into ngsParalog
      # - ngsParalog: probablistic call of paralogous regions
      # - After ngsParalog, calculate p-values based on chi-sq df = 1 with Benjamini-Hochberg correction
  # 4. Verify heterozygote excess is reduced
      # - Re-run ANGSD excluding paralogous regions and PCR duplicates
      # - Visualize SFS

## STEP 1: MERGE REPLICATES & MARK PCR & OPTICAL DUPLICATES 


# merge replicates if found, otherwise rename and move
rule merge_replicates:
  input:
    lambda wildcards: find_replicates(wildcards.sample_prefix, wildcards.hap)
  output:
    "results/bam_raw/hap{hap}/merged/{sample_prefix}_hap{hap}.bam"
  log:
    "results/logs/merge_replicates/hap{hap}/{sample_prefix}_hap{hap}.log"
  envmodules:
    "samtools/1.17"
  shell:
    """
    module load StdEnv/2020
    module load samtools/1.17
    echo "Input files: {input}" >> {log}
    echo "Output file: {output}" >> {log}
    
    if [ -z "{input}" ]; then
      echo "No files found for {wildcards.sample_prefix} in hap{wildcards.hap}" >> {log}
      exit 1
    elif [ $(echo {input} | wc -w) -eq 1 ]; then
      echo "Single file found for {wildcards.sample_prefix}, {wildcards.hap}, moving it." >> {log}
      mv {input} {output}
    else
      echo "Multiple files found for {wildcards.sample_prefix}, {wildcards.hap}. Merging..." >> {log}
      samtools merge -o {output} {input} 2>> {log}
    fi
    sync
    """




# Marking and removing PCR duplicates + index
# Note that: /bam_mkdup/ is where marked duplicates are marked but not removed. This is backed up in the projects folder.
# marked duplicates will be removed downstream in ANGSD
rule mark_remove_duplicates:
  input:
    raw_bam="results/bam_raw/hap{hap}/merged/{sample_prefix}_hap{hap}.bam"
  output:
    bam="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bam",
    bai="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bai",
    metrics="results/qc/mkdup_metrics/{sample_prefix}_hap{hap}.metrics"
  log:
    "results/logs/mark_remove_duplicates/hap{hap}/{sample_prefix}_hap{hap}.log"
  envmodules:
    "gatk/4.4.0.0"
  shell:
    """
    gatk MarkDuplicates \
    --CREATE_INDEX \
    -I {input.raw_bam} \
    -O {output.bam} \
    -M {output.metrics} \
    --REMOVE_DUPLICATES true \
    2> {log}
    """


# Clip overlapping reads
rule clip_overlapping_reads:
  input:
    bam="results/bam_mkdup/hap{hap}/{sample_prefix}_hap{hap}_mkdup.bam"
  output:
    clipped_bam="results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bam",
    clipped_index="results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bai"
  log:
    "results/logs/clip_overlap/hap{hap}/{sample_prefix}_clip_overlapping_reads.log"
  shell:
    """
    module --force purge
    module load nixpkgs/16.09 intel/2018.3
    module load bamutil/1.0.14
    bam clipOverlap --in {input.bam} --out {output.clipped_bam} 2> {log}
    module --force purge
    module load StdEnv/2023 samtools/1.18
    samtools index -b {output.clipped_bam} -o {output.clipped_index} --threads 2
    module --force purge
    """


# Create bam list per population for entry into realign indels
rule generate_clipped_bam_list_per_population:
  input:
    expand("results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bam", sample_prefix=sample_prefixes, hap=haps),
  output:
    "data/lists/hap{hap}/{population}_clipped_hap{hap}.list"
  wildcard_constraints:
    population="|".join(POPULATIONS)
  run:
    bam_files = input
    output_file = output[0]
    population = wildcards.population
    hap = wildcards.hap

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if population in bam_file and f"_hap{hap}_" in bam_file:
                output.write(f"{bam_file}\n")


# Create interval list of indels
rule indel_list:
  input:
    bam_list="data/lists/hap{hap}/{population}_clipped_hap{hap}.list",
    reference="data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    intervals="data/lists/hap{hap}/{population}_hap{hap}_indels.intervals"
  log:    
    "results/logs/indel_list/hap{hap}/{population}_hap{hap}_indel_list.log"
  shell:
    """
    module --force purge
    module load nixpkgs/16.09 gatk/3.8
    java -Xmx16g \
    -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R {input.reference} \
    -I {input.bam_list} \
    -o {output.intervals} \
    -drf BadMate \
    2> {log}
    module --force purge
    """


# Realign reads around indels - not entirely necessary if using ANGSD
rule realign_indels:
  input:
    bam="results/bam_clipped/hap{hap}/{sample_prefix}_hap{hap}_clipped.bam",
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    intervals=lambda wildcards: "data/lists/hap" + wildcards.hap + "/" + next(pop for pop in POPULATIONS if pop in wildcards.sample_prefix) + f"_hap{wildcards.hap}_indels.intervals"
  output:
    realigned_bam="results/bam_realign/hap{hap}/{sample_prefix}_hap{hap}_realign.bam"
  log:
    "results/logs/realign_indels/hap{hap}/{sample_prefix}_hap{hap}_realign_indels.log"
  shell:
    """
    module --force purge
    module load nixpkgs/16.09 gatk/3.8
    java -Xmx16g \
    -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -I {input.bam} \
    -R {input.ref} \
    -targetIntervals {input.intervals} \
    -o {output.realigned_bam}\
    &> {log}
    module --force purge
    """


# Rule to create .txt file of BAM files for each haplotype
rule generate_bam_list_per_haplotype:
  input:
    expand("results/bam_realign/hap{hap}/{sample_prefix}_hap{hap}_realign.bam", sample_prefix=sample_prefixes, hap=haps)
  output:
    "data/lists/hap{hap}/all_realign_hap{hap}.txt"
  run:
    bam_files = input
    output_file = output[0]
    hap = wildcards.hap

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if f"_hap{hap}_" in bam_file:
                output.write(f"{bam_file}\n")


# Merge all bams to estimate per sample depth
rule merge_bams:
  input:
    bam_list="data/lists/hap{hap}/all_realign_hap{hap}.txt"
  output:
    merged_bam="results/bam_merge/hap{hap}/merged_hap{hap}.bam"
  log:
    "results/logs/merge_bams/merge_bam_{hap}.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    samtools merge -b {input.bam_list} {output.merged_bam} 2> {log}
    """

rule index_bams:
  input:
    merged_bam="results/bam_merge/hap{hap}/merged_hap{hap}.bam"
  output:
    merged_bai="results/bam_merge/hap{hap}/merged_hap{hap}.bam.bai"
  log:
    "results/logs/merged_bams/index_bai_{hap}.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    samtools index {input.merged_bam} {output.merged_bai} 2> {log}
    """


# Estimate depth per bp per sample
# -d flag : set max-depth to infinity
# -A flag : don't skip anamalous read pairs marked w/ FLAG but missing properly paired flag set
# -q flag : minimum mapping quality
# -Q flag : minimum base quality
# -r flag : specify regions for pileup; needs indexed BAM files
# -a flag : all positions including unused reference sequences
# -f flag : faidx-indexed reference fasta file
# -ff flag : exclude SECONDARY, QCFAIL, DUP
# Note: previously used mpileup with all 189 .bam files, but now merged them as one file because the positions are different across samples with flag -s
# Using ngsQC tools 'bamstats', created by Dr. Tyler Linderoth

# run for only haplotype 2
rule estimate_depth_hap1:
  input:
    bam="results/bam_merge/hap1/merged_hap1.bam",
    bai="results/bam_merge/hap1/merged_hap1.bam.bai",
    ref="data/reference/hap1/lupinehap1.fasta"
  output:
    depth="results/depths/hap1/{hap1scaffold_prefix}_depth_est.txt.gz"
  params:
    hap1scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap1scaffold_prefix, 1)
  envmodules:
    "samtools/1.20"
  shell:
    """
    bamstats {input.bam} -r {params.hap1scaffold}\
    -A -d 77000000 -q 0 -Q 0 --ff UNMAP,SECONDARY,QCFAIL,DUP\
    -s -aa -f {input.ref} | gzip > {output.depth}
    """

rule estimate_depth_hap2:
  input:
    bam="results/bam_merge/hap2/merged_hap2.bam",
    bai="results/bam_merge/hap2/merged_hap2.bam.bai",
    ref="data/reference/hap2/lupinehap2.fasta"
  output:
    depth="results/depths/hap2/{hap2scaffold_prefix}_depth_est.txt.gz"
  params:
    hap2scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap2scaffold_prefix, 2)
  envmodules:
    "samtools/1.20"
  shell:
    """
    bamstats {input.bam} -r {params.hap2scaffold}\
    -A -d 77000000 -q 0 -Q 0 --ff UNMAP,SECONDARY,QCFAIL,DUP\
    -s -aa -f {input.ref} | gzip > {output.depth}
    """


# Plot aggregate depths per population
rule plot_depth_per_scaffold_hap1:
  input:
    depth_file="results/depths/hap1/{hap1scaffold_prefix}_depth_est.txt.gz"
  output:
    plot="results/plots/hap1/depths/{hap1scaffold_prefix}_depth_histogram.png"
  shell:
    """
    python scripts/plot_depths.py {input.depth_file} {output.plot}
    """

rule plot_depth_per_scaffold_hap2:
  input:
    depth_file="results/depths/hap2/{hap2scaffold_prefix}_depth_est.txt.gz"
  output:
    plot="results/plots/hap2/depths/{hap2scaffold_prefix}_depth_histogram.png"
  shell:
    """
    python scripts/plot_depths.py {input.depth_file} {output.plot}
    """

rule plot_aggregate_depths_hap1:
  input:
    aggregated_depth=expand("results/depths/hap1/{hap1scaffold_prefix}_depth_est.txt.gz", hap1scaffold_prefix=HAP1SCAFFOLD_PREFIXES)
  output:
    plot="results/plots/hap1/depths/hap1_aggregate_depth_histogram.png"
  shell:
    """
    python scripts/plot_depths_WG.py {input.aggregated_depth} {output.plot}
    """

rule plot_aggregate_depths_hap2:
  input:
    aggregated_depth=expand("results/depths/hap2/{hap2scaffold_prefix}_depth_est.txt.gz", hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES)
  output:
    plot="results/plots/hap2/depths/hap2_aggregate_depth_histogram.png"
  shell:
    """
    python scripts/plot_depths_WG.py {input.aggregated_depth} {output.plot}
    """










## STEP 2: Quality control check of HWE and SFS

# Rule to create .txt file of BAM files per population
rule generate_bam_list_per_population:
  input:
    expand("results/bam_realign/hap{hap}/{sample_prefix}_hap{hap}_realign.bam", sample_prefix=sample_prefixes, hap=haps),
  output:
    "data/lists/hap{hap}/{population}_realign_hap{hap}.txt"
  wildcard_constraints:
    population="|".join(POPULATIONS)
  run:
    bam_files = input
    output_file = output[0]
    population = wildcards.population
    hap = wildcards.hap

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if population in bam_file and f"_hap{hap}_" in bam_file:
                output.write(f"{bam_file}\n")

# ANGSD by population:  To calculate SFS without ngsParalog filtering
# We try different depth levels to get a better distribution of HWE. 

rule angsd_SFS_control_by_population_on_all_sites:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt"
  output:
    arg_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.arg",
    mafs_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.mafs.gz",
    depth_sample="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.depthSample",
    depth_global="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.depthGlobal",
    saf_1="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.saf.idx",
    saf_2="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.saf.pos.gz",
    saf_3="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.saf.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt",
    min_MAF=lambda wildcards: wildcards.min_MAF
    #min_depth=lambda wildcards: wildcards.min_depth,
    #max_depth=lambda wildcards: wildcards.max_depth
  log:
    "results/logs/angsd/hap{hap}/raw/by_popln/angsd_raw_sites_control_hap{hap}_{population}_minMAF{min_MAF}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 2\
    -C 50\
    -minMapQ 30\
    -minQ 20\
    -minInd 20\
    -minMAF {params.min_MAF}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doSaf 1\
    -anc {params.ref}\
    &> {log}
    """

# -minMaf 0.01 makes everything with a steep SFS

# Rule for generating global SFS with different depth settings
rule global_SFS_control_by_population:
  input:
    saf_idx="results/angsd/hap{hap}/raw/by_popln/{population}_raw_sites_control_minMAF{min_MAF}.saf.idx"
  output:
    sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_minMAF{min_MAF}.sfs"
  log:
    sfs1="results/logs/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_minMAF{min_MAF}.log"
  threads: 40
  shell:
    """
    winsfs {input.saf_idx} -t {threads} --seed 1 -v > {output.sfs} 2> {log.sfs1}
    """


rule fold_global_SFS_control_by_population:
  input:
    sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_minMAF{min_MAF}.sfs"
  output:
    sfs_folded="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_folded_minMAF{min_MAF}.sfs"
  log:
    sfs2="results/logs/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_folded_minMAF{min_MAF}.log"
  threads: 40
  shell:
    """
    winsfs view --fold {input.sfs} -v > {output.sfs_folded} 2> {log.sfs2}
    """


# Rule for plotting the folded SFS with different depth settings
rule global_SFS_control_by_population_plots:
  input:
    unfolded_sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_minMAF{min_MAF}.sfs",
    folded_sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_control_folded_minMAF{min_MAF}.sfs"
  output:
    unfolded_plot="results/plots/hap{hap}/SFS/{population}_globalSFS_control_unfolded_minMAF{min_MAF}.png",
    folded_plot="results/plots/hap{hap}/SFS/{population}_globalSFS_control_folded_minMAF{min_MAF}.png"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    """
    Rscript scripts/SFS_1D_graph.R {input.folded_sfs} {output.folded_plot}
    Rscript scripts/SFS_1D_graph.R {input.unfolded_sfs} {output.unfolded_plot}
    """


# Running ANGSD HWE analysis with different depth settings
rule angsd_HWE_by_population_on_control_SNPs:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.arg",
    mafs_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.mafs.gz",
    hwe_file="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.hwe.gz",
    depth_sample="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.depthSample",
    depth_global="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.depthGlobal"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt",
    min_depth=lambda wildcards: wildcards.min_depth,
    max_depth=lambda wildcards: wildcards.max_depth
  log:
    "results/logs/angsd/hap{hap}/raw/by_popln/angsd_raw_SNPs_control_hap{hap}_{population}_min{min_depth}_max{max_depth}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 2\
    -C 50\
    -minMapQ 30\
    -minQ 20 \
    -minMaf 0.01\
    -minHWEpval 0.01\
    -setMinDepth {params.min_depth}\
    -setMaxDepth {params.max_depth}\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -doHWE 1\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    &> {log}
    """


# Unzipping the HWE output with different depth settings
rule unzip_hwe_control:
  input:
    zipped_hwe="results/angsd/hap{hap}/raw/by_popln/{population}_raw_SNPs_control_min{min_depth}_max{max_depth}.hwe.gz"
  output:
    unzip_hwe="results/angsd/hap{hap}/raw/by_popln/{population}_hwe_raw_SNPs_control_min{min_depth}_max{max_depth}.lr"
  shell:
    """
    zcat {input.zipped_hwe} > {output.unzip_hwe}
    """

# Generating HWE histograms with different depth settings
rule hwe_histogram_control:
  input:
    lr_file="results/angsd/hap{hap}/raw/by_popln/{population}_hwe_raw_SNPs_control_min{min_depth}_max{max_depth}.lr"
  output:
    plot="results/plots/hap{hap}/HWE/{population}_hwe_raw_SNPs_control_min{min_depth}_max{max_depth}.png"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/hetero_excess_header.R {input.lr_file} {output.plot}"



## STEP 3: IDENTIFY PARALOGOUS REGIONS 

# Call SNPs with liberal rules for input to ngsParalog
# NOTE: -SNP_pvalue 0.05 so very liberal SNP calls
# NOTE: NO filters applied to gather as much as data as possible as input for ngsParalog
rule angsd_raw_SNP:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt"
  output:
    arg_file="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_realign.arg",
    mafs_file="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_realign.mafs.gz",
    hwe_file="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_realign.hwe.gz",
    depth_sample="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_realign.depthSample",
    depth_global="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_realign.depthGlobal"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_realign",
    hap="{hap}",
    population="{population}"
  log:
    "results/logs/angsd/hap{hap}/raw/ngsParalog_input/angsd_SNP_raw_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -C 50\
    -GL 2\
    -SNP_pval 0.05\
    -minMapQ 30\
    -minQ 20\
    -minInd 20 \
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -doHWE 1\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    &> {log}
    """


# Create BED files so including only SNPs into ngsParalog
# NOTE: BED files indexed at 0bp because SAMtools to create pileup requires 0bp index
rule convert_mafs_to_bed:
  input:
    mafs_gz="results/angsd/hap{hap}/raw/ngsParalog_input/{population}_raw_SNPs_realign.mafs.gz"
  output:
    bed_file="results/bed/hap{hap}/raw_SNPs/{population}_raw_SNPs_realign.BED"
  shell:
   """
   gunzip -c {input.mafs_gz} | awk 'NR>1 {{print $1, $2 - 1, $2}}' > {output.bed_file}
   dos2unix {output.bed_file}  # Add this line to convert line endings
   """


# Run ngsParalog on all SNPs across the genome (parallelized by running 1 job per scaffold)
rule ngsParalog_hap1: 
  input:
    bam_ngsPara=lambda wildcards: expand("results/bam_realign/hap1/{sample}_hap1_realign.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)]),
    ref="data/reference/hap1/lupinehap1.fasta",
    bed_file="results/bed/hap1/raw_SNPs/{population}_raw_SNPs_realign.BED"
  output:
    paralog_output="results/ngs_paralog/hap1/by_popln/{population}_scaffolds/{population}-{hap1scaffold_prefix}.lr"
  log:
    "results/logs/ngs_paralog/hap1/by_popln/{population}_scaffolds/{population}-{hap1scaffold_prefix}.log"
  params:
    hap1scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap1scaffold_prefix, 1)
  envmodules:
    "samtools/1.20"
  shell:
    """
    rm -f {output.paralog_output} #remove existing output file if it exists
    touch {output.paralog_output}
    samtools mpileup {input.bam_ngsPara} -A -d 77000000 -q 0 -Q 0 --ff UNMAP,QCFAIL,DUP \
    -l {input.bed_file} -r {params.hap1scaffold} -f {input.ref} 2>> {log} | \
    /home/socamero/ngsParalog/ngsParalog calcLR -infile - -outfile {output.paralog_output} -allow_overwrite 1 \
    -minQ 20 -minind 20 -mincov 1 \
    -runinfo 1 \
    2>> {log} || true
    """


rule ngsParalog_hap2:
  input:
    bam_ngsPara=lambda wildcards: expand("results/bam_realign/hap2/{sample}_hap2_realign.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)]),
    ref="data/reference/hap2/lupinehap2.fasta",
    bed_file="results/bed/hap2/raw_SNPs/{population}_raw_SNPs_realign.BED"
  output:
    paralog_output="results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.lr"
  log:
    "results/logs/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.log"
  params:
    hap2scaffold=lambda wildcards: map_prefix_to_full_scaffold(wildcards.hap2scaffold_prefix, 2)
  envmodules:
    "samtools/1.20"
  shell:
    """
    rm -f {output.paralog_output} #remove existing output file if it exists
    touch {output.paralog_output}
    samtools mpileup {input.bam_ngsPara} -A -d 77000000 -q 0 -Q 0 --ff UNMAP,QCFAIL,DUP \
    -l {input.bed_file} -r {params.hap2scaffold} -f {input.ref} 2>> {log} | \
    /home/socamero/ngsParalog/ngsParalog calcLR -infile - -outfile {output.paralog_output} -allow_overwrite 1 \
    -minQ 20 -minind 20 -mincov 1 \
    -runinfo 1 \
    2>> {log} || true
    """


# Combine all ngsParalog outputs together into one file per population
rule concatenate_paralog_hap1:
  input:
    scaffold_files=lambda wildcards: expand("results/ngs_paralog/hap1/by_popln/{population}_scaffolds/{{population}}-{hap1scaffold_prefix}.lr", population=wildcards.population, hap1scaffold_prefix=HAP1SCAFFOLD_PREFIXES)
  output:
    paralog_final="results/ngs_paralog/hap1/concat/{population}_paralog_realign_WGS.lr"
  log:
    "results/logs/ngs_paralog/hap1/concat/{population}_paralog_realign_WGS.log"
  shell:
    """
    cat {input.scaffold_files} >> {output.paralog_final}
    """


rule concatenate_paralog_hap2:
  input:
    scaffold_files=lambda wildcards: expand("results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{{population}}-{hap2scaffold_prefix}.lr", population=wildcards.population, hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES)
  output:
    paralog_final="results/ngs_paralog/hap2/concat/{population}_paralog_realign_WGS.lr"
  log:
    "results/logs/ngs_paralog/hap2/concat/{population}_paralog_realign_WGS.log"
  shell:
    """
    cat {input.scaffold_files} >> {output.paralog_final}
    """


# Print summary of calcLR quantile ranges
rule quantile_summary:
  input:
    "results/ngs_paralog/hap{hap}/concat/{population}_paralog_realign_WGS.lr"
  output:
    "results/ngs_paralog/hap{hap}/quantile_summary/{population}_calcLR_quantile_summary.txt"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/calcLR_quantile_summary.R {input} {output}"


# Identify false positives from ngsParalog (grab only true Paralogs) and filter out
# NOTE: R script indexes positions back to 1bp start
rule ngsParalog_false_pos:
  input:
    lr_file="results/ngs_paralog/hap{hap}/concat/{population}_paralog_realign_WGS.lr"
  output:
    deviant_snps="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected.BED",
    deviant_snps_bp1="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_realign_bp1_BH_corrected.lr"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/ngs_paralog_false_pos.R {input.lr_file} {output.deviant_snps} {output.deviant_snps_bp1}
    """


# Identify false positives for varying levels of Benjamini Hochberg critical values
rule ngsParalog_false_pos_BH:
  input:
    lr_file="results/ngs_paralog/hap{hap}/concat/{population}_paralog_realign_WGS.lr"
  output:
    deviant_snps_bp1_BH="results/bed/hap{hap}/deviant_SNPs/BH_correction/{population}_deviant_SNPs_bp1_realign_BH{BH_VAR}.lr"
  params:
    BH_VAR="{BH_VAR}"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/ngs_paralog_false_pos_BH.R {input.lr_file} {output.deviant_snps_bp1_BH} {params.BH_VAR}
    """


# Create Manhattan Plots [NOTE: Specific Scaffold # otherwise it'll take forever for whole genome!] 
rule ngsParalog_manhattan:
  input:
    lr_file="results/ngs_paralog/hap{hap}/by_popln/{population}_scaffolds/{population}-Scaffold_1.lr"
  output:
    plot="results/plots/hap{hap}/ngsParalog/by_popln/{population}_manhattan_plot_realign_scaffold_1.png"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    "Rscript scripts/ngs_paralog_graphs.R {input.lr_file} {output.plot}"


## dupHMM - Discover paralogous regions rather than individual sites ##
# Requires calcLR .lr files and depth of coverage .tsv files 

# Convert ngsParalog calcLR outputs (.lr files) to BED format
# This is NOT BH corrected. This is just here for future use
rule convert_calcLR_output_to_BED:
  input:
    lr_file="results/ngs_paralog/hap{hap}/concat/{population}_paralog_realign_WGS.lr",
  output:
    bed="results/bed/hap{hap}/deviant_SNPs_realign/{population}_deviant_SNPs_calcLR.BED"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    "Rscript scripts/convert_calcLR_to_bed.R {input.lr_file} {output.bed}"


# Generate average depth of coverage per site for dupHMM
# Only estimate depth at sites deemed 'paralogous' from raw ngsParalog (pre BH filter)
# Use pre BH filter because take all data as input!
rule calculate_average_coverage_per_population:
  input:
    bam_files=lambda wildcards: expand("results/bam_realign/hap{hap}/{sample}_hap{hap}_realign.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)], hap=haps),
    raw_lr="results/ngs_paralog/hap{hap}/concat/{population}_paralog_realign_WGS.lr"
  output:
    bed="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR.BED",
    avg_cov_file="results/coverage/hap{hap}/{population}_average_deviant_SNP_coverage.tsv"
  threads: 2
  log:
    "results/logs/coverage/hap{hap}/{population}_average_coverage.log"
  envmodules:
    "samtools/1.20"
  shell:
    """
    # Convert .lr file to BED format (0-based start)
    awk '{{print $1 "\\t" ($2-1) "\\t" $2}}' {input.raw_lr} > {output.bed}

    # Calculate the average coverage per position for all BAM files together using the calcLR BED file
    samtools depth -q 0 -Q 0 -J -a -b {output.bed} {input.bam_files} | \
    awk '{{cov[$1"\\t"$2]+=$3; count[$1"\\t"$2]++}} END {{for (pos in cov) print pos, cov[pos] / count[pos]}}' | \
    sort -k1,1V -k2,2n > {output.avg_cov_file}
    """

# Print summary statistics
rule calculate_coverage_statistics:
  input:
    avg_cov_file="results/coverage/hap{hap}/{population}_average_deviant_SNP_coverage.tsv"
  output:
    stats="results/coverage_stats/hap{hap}/{population}_coverage_stats.txt"
  shell:
    """
    python scripts/calculate_coverage_stats.py {input.avg_cov_file} > {output.stats}
    """


# Reduce deviant SNP calls (aka calcLR LRs) to the number of known sites with depth data.
# ANGSD estimates depth globally and per individual, but NOT per site (see -doDepth 1)
rule filter_paralog_by_coverage:
  input:
    lr_file="results/ngs_paralog/hap{hap}/concat/{population}_paralog_realign_WGS.lr",
    cov_file="results/coverage/hap{hap}/{population}_average_deviant_SNP_coverage.tsv"
  output:
    filtered_lr="results/bed/hap{hap}/deviant_SNPs/{population}_coverage_only_paralog_realign_WGS.lr"
  log:
    "results/logs/filter_lr/hap{hap}/{population}_filtered_deviant_SNPs_realign.log"
  shell:
    """
    # Filter the .lr file based on the coverage file, ensuring no empty lines
    awk 'NR==FNR && $0!="" {{cov[$1" "$2]=1; next}} $0!="" && cov[$1" "$2] {{print}}' {input.cov_file} {input.lr_file} > {output.filtered_lr}
    """

# dupHMM_setup is run across the entire genome Whereas dupHMM_run is across each scaffold. 
rule dupHMM_setup:
  input:
    lr_file="results/bed/hap{hap}/deviant_SNPs/{population}_coverage_only_paralog_realign_WGS.lr",
    cov_file="results/coverage/hap{hap}/{population}_average_deviant_SNP_coverage.tsv"
  output:
    param_output="results/ngs_paralog/hap{hap}/dupHMM/{population}_dupHMM_realign.par"
  params:
    r_script_path="/home/socamero/ngsParalog/dupHMM.R",
    name_output="results/ngs_paralog/hap{hap}/dupHMM/{population}_dupHMM_realign"
  log:
    "results/logs/ngs_paralog/hap{hap}/dupHMM/{population}_dupHMM_realign_setup.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript {params.r_script_path}
    --lrfile {input.lr_file} \
    --outfile {params.name_output} \
    --covfile {input.cov_file} \
    --lrquantile 0.975 \
    --paramOnly 1 \
    2> {log}
    """


# Now generate coverage per scaffold and per population
# We need to input coverage into dupHMM!
rule calculate_coverage_per_population_and_scaffold_hap1:
  input:
    bam_files=lambda wildcards: expand("results/bam_realign/hap1/{sample}_hap1_realign.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)]),
    raw_lr="results/ngs_paralog/hap1/by_popln/{population}_scaffolds/{population}-{hap1scaffold_prefix}.lr"
  output:
    bed="results/bed/hap1/deviant_SNPs/{population}_scaffolds/{population}_{hap1scaffold_prefix}.BED",
    avg_cov_file="results/coverage/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_average_coverage.tsv"
  log:
    "results/logs/coverage/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_average_coverage.log"
  threads: 2
  envmodules:
    "samtools/1.20"
  shell:
    """
    # Check if the .lr file is empty and create a BED file accordingly
    if [ -s {input.raw_lr} ]; then
      # Convert non-empty .lr file to BED format (0-based start)
      awk '{{print $1 "\\t" ($2-1) "\\t" $2}}' {input.raw_lr} > {output.bed}
    else
      # Create an empty BED file if .lr file is empty
      touch {output.bed}
    fi

    # If the BED file is not empty, calculate the average coverage
    if [ -s {output.bed} ]; then
      samtools depth -@ 2 -q 0 -Q 0 -J -a -b {output.bed} {input.bam_files} | \
      awk '{{cov[$1"\\t"$2]+=$3; count[$1"\\t"$2]++}} END {{for (pos in cov) print pos, cov[pos] / count[pos]}}' | \
      sort -k1,1V -k2,2n > {output.avg_cov_file}
    else
      # Create an empty average coverage file if BED file is empty
      touch {output.avg_cov_file}
    fi
    """


rule calculate_coverage_per_population_and_scaffold_hap2:
  input:
    bam_files=lambda wildcards: expand("results/bam_realign/hap2/{sample}_hap2_realign.bam", sample=[s for s in sample_prefixes if s.startswith(wildcards.population)]),
    raw_lr="results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.lr"
  output:
    bed="results/bed/hap2/deviant_SNPs/{population}_scaffolds/{population}_{hap2scaffold_prefix}.BED",
    avg_cov_file="results/coverage/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_average_coverage.tsv"
  log:
    "results/logs/coverage/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_average_coverage.log"
  threads: 2
  envmodules:
    "samtools/1.20"
  shell:
    """
    # Check if the .lr file is empty and create a BED file accordingly
    if [ -s {input.raw_lr} ]; then
      # Convert non-empty .lr file to BED format (0-based start)
      awk '{{print $1 "\\t" ($2-1) "\\t" $2}}' {input.raw_lr} > {output.bed}
    else
      # Create an empty BED file if .lr file is empty
      touch {output.bed}
    fi

    # If the BED file is not empty, calculate the average coverage
    if [ -s {output.bed} ]; then
      samtools depth -@ 2 -q 0 -Q 0 -J -a -b {output.bed} {input.bam_files} | \
      awk '{{cov[$1"\\t"$2]+=$3; count[$1"\\t"$2]++}} END {{for (pos in cov) print pos, cov[pos] / count[pos]}}' | \
      sort -k1,1V -k2,2n > {output.avg_cov_file}
    else
      # Create an empty average coverage file if BED file is empty
      touch {output.avg_cov_file}
    fi
    """

# Filter .lr files by which there is depth data
# This is similar as before, but we need it by scaffold since...
# dupHMM setup is run across the entire genome ! Where as dupHMM run is run across each scaffold instead. 
rule filter_scaffold_lr_by_coverage_hap1:
  input:
    lr_file="results/ngs_paralog/hap1/by_popln/{population}_scaffolds/{population}-{hap1scaffold_prefix}.lr",
    cov_file="results/coverage/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_average_coverage.tsv"
  output:
    filtered_lr="results/bed/hap1/deviant_SNPs/{population}_scaffolds/{population}_{hap1scaffold_prefix}_coverage_only_paralog_realign.lr"
  log:
    "results/logs/filter_lr/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_filtered_deviant_SNPs_realign.log"
  shell:
    """
    # Filter the .lr file based on the coverage file, ensuring no empty lines
    awk 'NR==FNR && $0!="" {{cov[$1" "$2]=1; next}} $0!="" && cov[$1" "$2] {{print}}' {input.cov_file} {input.lr_file} > {output.filtered_lr}

    # If the filtered file is empty, create an empty file
    if [ ! -s {output.filtered_lr} ]; then
      touch {output.filtered_lr}
    fi
    """


rule filter_scaffold_lr_by_coverage_hap2:
  input:
    lr_file="results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.lr",
    cov_file="results/coverage/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_average_coverage.tsv"
  output:
    filtered_lr="results/bed/hap2/deviant_SNPs/{population}_scaffolds/{population}_{hap2scaffold_prefix}_coverage_only_paralog_realign.lr"
  log:
    "results/logs/filter_lr/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_filtered_deviant_SNPs_realign.log"
  shell:
    """
    # Filter the .lr file based on the coverage file, ensuring no empty lines
    awk 'NR==FNR && $0!="" {{cov[$1" "$2]=1; next}} $0!="" && cov[$1" "$2] {{print}}' {input.cov_file} {input.lr_file} > {output.filtered_lr}

    # If the filtered file is empty, create an empty file
    if [ ! -s {output.filtered_lr} ]; then
      touch {output.filtered_lr}
    fi
    """

# Run dupHMM with beginning estimated parameters
# adjust --lrquantile based on manhattan plot! 
# adjust --maxcoverage based on .tsv files! Seems like some coverages are HIGH! ~50.. To be slightly conservative, choosing 25 (see coverage stats)
# NOTE: Some scaffolds only have 1 SNP with coverage and/or paralog data and can't be used in dupHMM so we skip these. 
rule dupHMM_run_hap1:
  input:
    lr_file = "results/ngs_paralog/hap1/by_popln/{population}_scaffolds/{population}-{hap1scaffold_prefix}.lr",
    cov_file = "results/coverage/hap1/{population}_scaffolds/{population}_{hap1scaffold_prefix}_average_coverage.tsv",
    parameters = "results/ngs_paralog/hap1/dupHMM/{population}_dupHMM_realign.par"
  output:
    paralog_region = "results/ngs_paralog/hap1/dupHMM/{population}_scaffolds/{population}-{hap1scaffold_prefix}_dupHMM_realign_run.rf"
  params:
    r_script_path = "/home/socamero/ngsParalog/dupHMM.R",
    name_output = "results/ngs_paralog/hap1/dupHMM/{population}_scaffolds/{population}-{hap1scaffold_prefix}_dupHMM_realign_run"
  log:
    "results/logs/ngs_paralog/hap1/dupHMM/{population}_scaffolds/{population}-{hap1scaffold_prefix}_dupHMM_realign_run.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    if [ -s {input.lr_file} ] && [ $(cat {input.lr_file} | wc -l) -gt 1 ]; then
        Rscript {params.r_script_path} --lrfile {input.lr_file} \
        --outfile {params.name_output} \
        --covfile {input.cov_file} \
        --lrquantile 0.97 \
        --maxcoverage 25 \
        --paramfile {input.parameters} \
        2> {log} || touch {output.paralog_region}
    else
        echo "Skipping .lr file due to insufficient data (empty or only one row): {input.lr_file}" >> {log}
        touch {output.paralog_region}
    fi
    """

# NOTE: coverage is lower using hap2
rule dupHMM_run_hap2:
  input:
    lr_file = "results/ngs_paralog/hap2/by_popln/{population}_scaffolds/{population}-{hap2scaffold_prefix}.lr",
    cov_file = "results/coverage/hap2/{population}_scaffolds/{population}_{hap2scaffold_prefix}_average_coverage.tsv",
    parameters = "results/ngs_paralog/hap2/dupHMM/{population}_dupHMM_realign.par"
  output:
    paralog_region = "results/ngs_paralog/hap2/dupHMM/{population}_scaffolds/{population}-{hap2scaffold_prefix}_dupHMM_realign_run.rf"
  params:
    r_script_path = "/home/socamero/ngsParalog/dupHMM.R",
    name_output = "results/ngs_paralog/hap2/dupHMM/{population}_scaffolds/{population}-{hap2scaffold_prefix}_dupHMM_realign_run"
  log:
    "results/logs/ngs_paralog/hap2/dupHMM/{population}_scaffolds/{population}-{hap2scaffold_prefix}_dupHMM_realign_run.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    if [ -s {input.lr_file} ] && [ $(cat {input.lr_file} | wc -l) -gt 1 ]; then
        Rscript {params.r_script_path} --lrfile {input.lr_file} \
        --outfile {params.name_output} \
        --covfile {input.cov_file} \
        --lrquantile 0.97 \
        --maxcoverage 20 \
        --paramfile {input.parameters} \
        2> {log} || touch {output.paralog_region}
    else
        echo "Skipping .lr file due to insufficient data (empty or only one row): {input.lr_file}" >> {log}
        touch {output.paralog_region}
    fi
    """


# Combine all dupHMM outputs together into one file per population
rule concatenate_dupHMM_hap1:
  input:
    dupHMM_files=lambda wildcards: expand("results/ngs_paralog/hap1/dupHMM/{population}_scaffolds/{{population}}-{hap1scaffold_prefix}_dupHMM_realign_run.rf", population=wildcards.population, hap1scaffold_prefix=HAP1SCAFFOLD_PREFIXES)
  output:
    dupHMM_final="results/ngs_paralog/hap1/dupHMM/{population}_dupHMM_realign_WGS.lr"
  log:
    "results/logs/ngs_paralog/hap1/dupHMM/{population}_dupHMM_realign_WGS.log"
  shell:
    """
    cat {input.dupHMM_files} >> {output.dupHMM_final}
    """


rule concatenate_dupHMM_hap2:
  input:
    dupHMM_files=lambda wildcards: expand("results/ngs_paralog/hap2/dupHMM/{population}_scaffolds/{{population}}-{hap2scaffold_prefix}_dupHMM_realign_run.rf", population=wildcards.population, hap2scaffold_prefix=HAP2SCAFFOLD_PREFIXES)
  output:
    dupHMM_final="results/ngs_paralog/hap2/dupHMM/{population}_dupHMM_realign_WGS.lr"
  log:
    "results/logs/ngs_paralog/hap2/dupHMM/{population}_dupHMM_realign_WGS.log"
  shell:
    """
    cat {input.dupHMM_files} >> {output.dupHMM_final}
    """


# Convert dupHMM outputs (.lr files) to BED format
rule convert_dupHMM_output_to_BED:
  input:
    lr_file="results/ngs_paralog/hap{hap}/dupHMM/{population}_dupHMM_realign_WGS.lr"
  output:
    bed="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM.BED"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    "Rscript scripts/convert_dupHMM_to_bed.R {input.lr_file} {output.bed}"


## STEP 4: VERIFY PARALOGS REDUCED -> look at SFS with RE-RUN ANGSD

# Print out all sites (monomorphic or not) to later filter out paralogous regions
rule angsd_raw_sites_by_popln:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt"
  output:
    all_sites_gz="results/angsd/hap{hap}/raw/by_popln/{population}_all_sites.pos.gz",
    all_sites_arg="results/angsd/hap{hap}/raw/by_popln/{population}_all_sites.arg",
    all_sites_bed="results/bed/hap{hap}/raw_sites/{population}_all_sites.BED"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    angsd_out="results/angsd/hap{hap}/raw/by_popln/{population}_all_sites"
  log:
    angsd_log="results/logs/angsd/hap{hap}/raw/by_popln/{population}_all_sites.log"
  threads: 8
  envmodules:
    "angsd/0.940"
  shell:
    """
    # Extract all sites across the genome
    angsd -bam {input.bam_list} -ref {params.ref} -out {params.angsd_out} \
    -doCounts 1 -dumpCounts 1 -P 8 &> {log.angsd_log}

    # Convert the ANGSD output to a BED format file
    gunzip -c {params.angsd_out}.pos.gz | awk 'NR > 1 {{print $1, $2-1, $2}}' > {output.all_sites_bed}
    """


# Format all BED files accordingly and make sure they're in UNIX format + tab delimited
rule preprocess_bed_files_all_populations:
  input: 
    all_sites_bed="results/bed/hap{hap}/raw_sites/{population}_all_sites.BED",
    calcLR_deviant_snps_bed="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected.BED",
    dupHMM_deviant_regs_bed="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM.BED"
  output:
    processed_all_sites_bed="results/bed/hap{hap}/raw_sites/by_popln/{population}_all_sites_processed.BED",
    processed_calcLR_snps="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected_processed.BED",
    processed_dupHMM_regs="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM_processed.BED"
  shell:
    """
    # Convert all_sites_bed BED file from DOS to UNIX line endings and process
    dos2unix {input.all_sites_bed}
    awk 'NF >= 3 {{print $1"\t"$2"\t"$3}}' {input.all_sites_bed} > {output.processed_all_sites_bed}

    # Convert calcLR_snps BED file from DOS to UNIX line endings and process
    dos2unix {input.calcLR_deviant_snps_bed}
    awk 'NF >= 3 {{print $1"\t"$2"\t"$3}}' {input.calcLR_deviant_snps_bed} > {output.processed_calcLR_snps}

    # Convert dupHMM_regions BED file from DOS to UNIX line endings and process
    dos2unix {input.dupHMM_deviant_regs_bed}
    awk 'NF >= 3 {{print $1"\t"$2"\t"$3}}' {input.dupHMM_deviant_regs_bed} > {output.processed_dupHMM_regs}
    

    # Optional: Print lines before and after removal for auditing
    echo "Lines before processing in {input.all_sites_bed}:"
    wc -l {input.all_sites_bed}
    echo "Lines after processing in {output.processed_all_sites_bed}:"
    wc -l {output.processed_all_sites_bed}

    echo "Lines before processing in {input.calcLR_deviant_snps_bed}:"
    wc -l {input.calcLR_deviant_snps_bed}
    echo "Lines after processing in {output.processed_calcLR_snps}:"
    wc -l {output.processed_calcLR_snps}

    echo "Lines before processing in {input.dupHMM_deviant_regs_bed}:"
    wc -l {input.dupHMM_deviant_regs_bed}
    echo "Lines after processing in {output.processed_dupHMM_regs}:"
    wc -l {output.processed_dupHMM_regs}
    """


# Filter out deviant SNPs from all known sites
rule filter_all_sites_by_popln_calcLR:
  input:
    processed_all_sites_bed="results/bed/hap{hap}/raw_sites/by_popln/{population}_all_sites_processed.BED",
    processed_calcLR_snps="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected_processed.BED",
  output:
    filtered_calcLR_bed="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.BED",
    filtered_calcLR_txt="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.txt"
  log:
    calcLR_log="results/logs/bedtools/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_sites.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Filter out deviant calcLR SNPs using bedtools and save to BED
    bedtools subtract -a {input.processed_all_sites_bed} -b {input.processed_calcLR_snps} > {output.filtered_calcLR_bed} 2> {log.calcLR_log}

    # Convert the filtered BED files to a .txt file formatted for -sites in ANGSD
    awk '{{print $1, $3}}' {output.filtered_calcLR_bed} > {output.filtered_calcLR_txt}
    """


rule filter_all_sites_by_popln_dupHMM:
  input:
    processed_dupHMM_regs="results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM_processed.BED",
    filtered_calcLR_bed="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.BED"
  output:
    filtered_dupHMM_bed="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.BED",
    filtered_dupHMM_txt="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt"
  log:
    dupHMM_log="results/logs/bedtools/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Filter out dupHMM regions from the calcLR BED file
    bedtools subtract -a {input.filtered_calcLR_bed} -b {input.processed_dupHMM_regs} > {output.filtered_dupHMM_bed} 2> {log.dupHMM_log}

    # Convert the filtered BED files to a .txt file formatted for -sites in ANGSD
    awk '{{print $1, $3}}' {output.filtered_dupHMM_bed} > {output.filtered_dupHMM_txt}
    """


# Index filtered all sites BED file
rule index_all_sites_by_popln_calcLR:
  input: 
    canonical_calcLR_sites="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.txt"
  output: 
    calcLR_index="results/bed/hap{hap}/canonical_sites/filtered_calcLR/{population}_filtered_calcLR_sites.txt.bin"
  envmodules:
    "angsd/0.940"
  shell: 
    "angsd sites index {input.canonical_calcLR_sites}"


rule index_all_sites_by_popln_dupHMM:
  input:
    canonical_dupHMM_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt"
  output:
    dupHMM_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin"
  envmodules:
    "angsd/0.940"
  shell:
    "angsd sites index {input.canonical_dupHMM_sites}"


# ANGSD by population: To calculate SFS (check if heterozygote excess reduced) with filtered sites by ngsParalog
# Previously attempted just calcLR outputs but this did not remove paralogs. We now incorporate dupHMM and calcLR
rule angsd_SFS_by_population_on_all_sites:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin"
  output:
    arg_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.arg",
    mafs_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.mafs.gz",
    depth_sample="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.depthSample",
    depth_global="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.depthGlobal",
    saf_1="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.saf.idx",
    saf_2="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.saf.pos.gz",
    saf_3="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.saf.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt"
  log:
    "results/logs/angsd/hap{hap}/canonical/by_popln/angsd_canonical_sites_dupHMM_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 2\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -minMapQ 30\
    -minQ 20\
    -minInd 20\
    -minMaf 0.01\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doSaf 1\
    -anc {params.ref}\
    &> {log}
    """


# FILTERED SFS: Optimize and calculate SFS with folded spectra (--fold) since ancestral state unknown
# Can bootstrap to get confidence intervals
# We use 'winsfs' rather than 'realSFS' because realSFS is computational heavy (>100GB RAM required) plus less accurate
# Based on winsfs github, we allocate ~150GB RAM for ~1B sites (extracted from gzip and wc -l)
# We set seed at 1 for reproducibility
rule global_SFS_by_population:
  input:
    saf_idx="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_sites_dupHMM.saf.idx"
  output:
    sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM.sfs"
  log:
    sfs1="results/logs/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM.log"
  threads: 40
  shell:
    """
    winsfs {input.saf_idx} -t {threads} --seed 1 -v > {output.sfs} 2> {log.sfs1}
    """


rule fold_global_SFS_by_population:
  input:
    sfs="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM.sfs"
  output:
    sfs_folded="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM_folded.sfs"
  log:
    sfs2="results/logs/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM_folded.log"
  threads: 40
  shell:
    """
    winsfs view --fold {input.sfs} -v > {output.sfs_folded} 2> {log.sfs2}
    """


# Create SFS Plots on filtered data
rule global_SFS_by_population_plots:
  input:
    sfs_file="results/winSFS/hap{hap}/globalSFS/{population}_hap{hap}_globalSFS_dupHMM_folded.sfs"
  output:
    plot="results/plots/hap{hap}/SFS/{population}_globalSFS_dupHMM_folded.png"
  envmodules:
    "r/4.4.0"
  threads: 2
  shell:
    "Rscript scripts/SFS_1D_graph.R {input.sfs_file} {output.plot}"


# To check F distribution on filtered SNPs using dupHMM and calcLR
rule angsd_HWE_by_population_on_dupHMM_SNPs:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.arg",
    mafs_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.mafs.gz",
    hwe_file="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.hwe.gz",
    depth_sample="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.depthSample",
    depth_global="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.depthGlobal",
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt"
  log:
    "results/logs/angsd/hap{hap}/canonical/by_popln/angsd_canonical_SNPs_dupHMM_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 2\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -minMapQ 30\
    -minQ 20\
    -minInd 20\
    -minHWEpval 0.01\
    -minMaf 0.01\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -doHWE 1\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    &> {log}
    """


# Unzip hwe to extract desired variables (i.e F values)
rule unzip_hwe_dupHMM:
  input:
    zipped_hwe="results/angsd/hap{hap}/canonical/by_popln/{population}_canonical_SNPs_dupHMM.hwe.gz"
  output:
    unzip_hwe="results/angsd/hap{hap}/canonical/by_popln/{population}_hwe_canonical_SNPs_dupHMM.lr"
  shell:
    """
    zcat {input.zipped_hwe} > {output.unzip_hwe}
    """


# Create histogram of F values for canonical SNPs (excludes header)
rule hwe_histogram_dupHMM:
  input:
    lr_file="results/angsd/hap{hap}/canonical/by_popln/{population}_hwe_canonical_SNPs_dupHMM.lr"
  output:
    plot="results/plots/hap{hap}/HWE/{population}_hwe_canonical_SNPs_dupHMM.png"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/hetero_excess_header.R {input.lr_file} {output.plot}"


#####################################
#### POPULATION GENOMIC ANALYSES ####
#####################################


## ANALYSIS 1 : Identify potential clones using ngsRelate
# NOTE: set -minMaf to 0.001 which is quite liberal compared to 0.01.
# Previously, did not notice much difference in SFS plots between 0.01 and 0.0001
# ANGSD with -doGlf 3 to prepare input for ngsRelate
# Estimate allele frequencies, genotype likelihoods, and call SNPs
rule angsd_for_ngsRelate:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.arg",
    mafs_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.mafs.gz",
    hwe_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.hwe.gz",
    glf_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.glf.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt"
  log:
    "results/logs/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/angsd_canonical_SNPs_dupHMM_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 2\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -minMapQ 30\
    -minQ 20\
    -minInd 20\
    -minHWEpval 0.01\
    -minMaf 0.0001\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -doHWE 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doCounts 1\
    -doGlf 3\
    &> {log}
    """
  

rule ngsRelate_prep:
  input:
    mafs_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.mafs.gz",
  output:
    freq_file="results/ngsRelate/hap{hap}/{population}_freq"
  envmodules:
    "gcc",
    "htslib"
  log:
    "results/logs/ngsRelate/hap{hap}/{population}_freq_extraction.log"
  shell:
    """
    mkdir -p $(dirname {output.freq_file})
    zcat {input.mafs_file} | cut -f6 | sed 1d > {output.freq_file} 2>{log}
    """


# Identify potential clones using ngsRelate
# Will need to change -n 27 according to the number of samples per population (create some list)
rule ngsRelate_analysis:
  input:
    glf_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.glf.gz",
    freq_file="results/ngsRelate/hap{hap}/{population}_freq"
  output:
    ngsrelate_output="results/ngsRelate/hap{hap}/{population}_ngsrelate.out"
  log:
    "results/logs/ngsRelate/hap{hap}/{population}_ngsrelate_analysis.log"
  threads: 4
  shell:
    """
    ngsRelate -g {input.glf_file} \
              -n 27 \
              -f {input.freq_file} \
              -O {output.ngsrelate_output} \
              &> {log}
    """


rule plot_ngsRelate_output:
  input:
    ngsrelate_out=expand("results/ngsRelate/hap2/{population}_ngsrelate.out", population=POPULATIONS)
  output:
    boxplot="results/plots/hap2/ngsRelate/rab_boxplot.png"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_ngsRelate.R
    """


rule detect_clones:
  input:
    ngsrelate_out="results/ngsRelate/hap{hap}/{population}_ngsrelate.out"
  output:
    clone_summary="results/ngsRelate/hap{hap}/{population}_clone_summary.txt"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/detect_clones.R {input.ngsrelate_out} {output.clone_summary}"


rule ngsRelate_analysis_inbreeding:
  input:
    glf_file="results/angsd/hap{hap}/canonical/ngsRelate_input/by_popln/{population}_canonical_SNPs_dupHMM.glf.gz",
    freq_file="results/ngsRelate/hap{hap}/{population}_freq"
  output:
    ngsrelate_output="results/ngsRelate/hap{hap}/{population}_inbreeding.out"
  log:
    "results/logs/ngsRelate/hap{hap}/{population}_ngsrelate_inbreeding_analysis.log"
  threads: 4
  shell:
    """
    ngsRelate -g {input.glf_file} \
              -n 27 \
              -F 1\
              -f {input.freq_file} \
              -O {output.ngsrelate_output} \
              &> {log}
    """

rule plot_ngsRelate_inbreeding:
  input:
    ngsrelate_out=expand("results/ngsRelate/hap2/{population}_inbreeding.out", population=POPULATIONS)
  output:
    boxplot="results/plots/hap2/ngsRelate/F_boxplot.png"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_ngsRelate_inbreeding.R
    """


# ANALYSIS 2: individual-level inbreeding

# ANGSD with -doGlf 3 to prepare input for ngsF
# Estimate allele frequencies, genotype likelihoods, call SNPs, and get site frequency spectra
# NOTE: ngsF requires variable sites only, so call SNPs
rule angsd_for_ngsF:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin"
  output:
    arg_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.arg",
    mafs_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.mafs.gz",
    hwe_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.hwe.gz",
    saf_1="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.saf.idx",
    saf_2="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.saf.pos.gz",
    saf_3="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.saf.gz",
    glf_gz="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.gz",
    glf_pos="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.pos.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt"
  log:
    "results/logs/angsd/hap{hap}/canonical/ngsF_input/by_popln/angsd_canonical_SNPs_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 2\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -minMapQ 30\
    -minQ 20\
    -minInd 20\
    -minMaf 0.0001\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -doHWE 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doSaf 1\
    -doCounts 1\
    -anc {params.ref}\
    -doGlf 3\
    &> {log}
    """


# Grab number of variable sites (SNPs) per population
rule count_variable_sites:
  input:
    mafs_file="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.mafs.gz"
  output:
    site_count="results/ngsF/hap{hap}/by_popln/{population}_SNPs_count.txt"
  shell:
    """
    zcat {input.mafs_file} | wc -l | awk '{{print $1-1}}' > {output.site_count}
    """


# Estimate inbreeding coefficients
# For low-coverage data, set min_epsilon for lower threshold b/w 1e-5 and 1e-9 so algorithm keeps exploring before stopping
rule ngsF_analysis:
  input:
    GL3="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.gz",
    SNP_count="results/ngsF/hap{hap}/by_popln/{population}_SNPs_count.txt"
  output:
    ngsF_est="results/ngsF/hap{hap}/by_popln/{population}_ngsF_inbreeding_final.lrt"
  params:
    ngsF_output_base="results/ngsF/hap{hap}/by_popln/{population}_ngsF_inbreeding_iter",
    pop_size=27,
    n_sites=lambda wildcards: int(open(f"results/ngsF/hap{wildcards.hap}/by_popln/{wildcards.population}_SNPs_count.txt").read().strip()),
    iterations=5  # Number of iterations
  log:
    "results/logs/ngsF/hap{hap}/by_popln/{population}_inbreeding_estimate.log"
  threads: 4
  shell:
    """
    # Run the first iteration without the --init_values parameter
    echo "Running ngsF iteration 1"
    zcat {input.GL3} |\
    ngsF --glf -\
         --n_threads {threads}\
         --calc_LRT 1\
         --out {params.ngsF_output_base}_1\
         --n_ind {params.pop_size}\
         --n_sites {params.n_sites}\
         --init_values r\
         --min_epsilon 1e-7\
         &>> {log}

    # Loop over the remaining iterations
    for iter in $(seq 2 {params.iterations}); do
        echo "Running ngsF iteration $iter"

        # Use the .pars file from the previous iteration as the initial values for the current iteration
        zcat {input.GL3} |\
        ngsF --glf -\
             --n_threads {threads}\
             --calc_LRT 1\
             --out {params.ngsF_output_base}_$iter\
             --n_ind {params.pop_size}\
             --n_sites {params.n_sites}\
             --init_values {params.ngsF_output_base}_$((iter - 1)).pars\
             --min_epsilon 1e-7\
             &>> {log}
    done

    # Copy the final iteration output to the expected result
    cp {params.ngsF_output_base}_{params.iterations}.lrt {output.ngsF_est}
    """


# Use ngsF-HMM which uses a 2-step hidden-markov model
rule ngsF_HMM_analysis:
  input:
    GL3="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.gz",
    GL3_pos="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.pos.gz",
    SNP_count="results/ngsF/hap{hap}/by_popln/{population}_SNPs_count.txt"
  output:
    ngsFHMM_est="results/ngsF/hap{hap}/by_popln/{population}_ngsF-HMM_inbreeding.indF",
    GL3_unzipped="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf",
    GL3_pos_unzipped="results/angsd/hap{hap}/canonical/ngsF_input/by_popln/{population}_canonical_SNPs.glf.pos"
  params:
    ngsFHMM_output_base="results/ngsF/hap{hap}/by_popln/{population}_ngsF-HMM_inbreeding",
    pop_size=27,
    n_sites=lambda wildcards: int(open(f"results/ngsF/hap{wildcards.hap}/by_popln/{wildcards.population}_SNPs_count.txt").read().strip())
  log:
    "results/logs/ngsF/hap{hap}/by_popln/{population}_inbreeding_HMM_estimate.log"
  threads: 4
  shell:
    """
    zcat {input.GL3} > {output.GL3_unzipped}
    zcat {input.GL3_pos} > {output.GL3_pos_unzipped}
    ngsF-HMM --geno {output.GL3_unzipped}\
         --n_threads {threads}\
         --pos {output.GL3_pos_unzipped}\
         --out {params.ngsFHMM_output_base}\
         --n_ind {params.pop_size}\
         --n_sites {params.n_sites}\
         --freq r\
         --indF r\
         --loglkl\
         --min_epsilon 1e-7\
         --seed 12345\
         --log 1\
         &>> {log}
    """

rule plot_ngsF_HMM:
  input:
    ngsFHMM_est=expand("results/ngsF/hap2/by_popln/{population}_ngsF-HMM_inbreeding.indF", population=POPULATIONS)
  output:
    plot="results/plots/hap2/ngsF/ngsF-HMM_inbreeding_coeff.tiff"
  log:
    "results/logs/ngsF/hap2/plot_ngsF_HMM_metrics.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_ngsF_HMM.R 
    """




# ANALYSIS 3: PCA for population structure 
# Population structure uses all data from all populations
# Therefore, this requires consolidating all .BAM files as well the .BED files from ngsParalog for filtering

rule generate_bam_list_all_populations:
  input:
    expand("results/bam_realign/hap{hap}/{sample}_hap{hap}_realign.bam", sample=sample_prefixes, hap=(1,2)),
  output:
    "data/lists/hap{hap}/all_populations_realign_hap{hap}.txt"
  run:
    bam_files = input
    output_file = output[0]
    hap = wildcards.hap

    with open(output_file, "w") as output:
        for bam_file in bam_files:
            if f"_hap{hap}_" in bam_file:
                output.write(f"{bam_file}\n")


rule combine_population_calcLR_bed_files:
  input:
    lambda wildcards: expand("results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_calcLR_BH_corrected.BED", 
                             hap=wildcards.hap, 
                             population=POPULATIONS),
  output:
    "results/bed/hap{hap}/deviant_sites/hap{hap}_combined_deviant_SNPs_realign_BH_correction.BED"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Combine all BED files for a specific haplotype and remove duplicates
    bedtools merge -i <(sort -k1,1 -k2,2n $(echo {{' '.join(input)}})) > {output}
    """


rule combine_population_dupHMM_bed_files:
  input:
    lambda wildcards: expand("results/bed/hap{hap}/deviant_SNPs/{population}_deviant_SNPs_dupHMM.BED", 
                             hap=wildcards.hap, 
                             population=POPULATIONS),
  output:
    "results/bed/hap{hap}/deviant_sites/hap{hap}_combined_dupHMM_regions.BED"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Combine all BED files for a specific haplotype and remove duplicates
    bedtools merge -i <(sort -k1,1 -k2,2n $(echo {{' '.join(input)}})) > {output}
    """


# Extract all known sites from ALL populations. This is to create list sites and later filter from paralogs
rule angsd_raw_sites_all_poplns:
  input:
    bam_list="data/lists/hap{hap}/all_populations_realign_hap{hap}.txt"
  output:
    all_sites_gz="results/angsd/hap{hap}/raw/all_poplns/all_sites.pos.gz",
    all_sites_arg="results/angsd/hap{hap}/raw/all_poplns/all_sites.arg"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    angsd_out="results/angsd/hap{hap}/raw/all_poplns/all_sites"
  log:
    angsd_log="results/logs/angsd/hap{hap}/raw/all_poplns/all_sites.log"
  threads: 8
  envmodules:
    "angsd/0.940"
  shell:
    """
    # Extract all sites across the genome for all populations
    angsd -bam {input.bam_list} -ref {params.ref} -out {params.angsd_out} \
    -doCounts 1 -dumpCounts 1 -P {threads} &> {log.angsd_log}
    """


rule convert_raw_sites_all_poplns:
  input:
    all_sites_gz="results/angsd/hap{hap}/raw/all_poplns/all_sites.pos.gz",
    all_sites_arg="results/angsd/hap{hap}/raw/all_poplns/all_sites.arg"
  output:
    all_sites_bed="results/bed/hap{hap}/raw_sites/all_poplns/all_sites.BED"
  threads: 4
  shell:
    """    
    # Convert the ANGSD output to a BED format file and check for completeness of data
    gzip -cd {input.all_sites_gz} | awk 'NR > 1 && NF >= 3 {{print $1"\t"$2-1"\t"$2}}' > {output.all_sites_bed}
    dos2unix {output.all_sites_bed}
    """


rule filter_all_sites_all_populations_calcLR:
  input:
    all_sites_bed="results/bed/hap{hap}/raw_sites/all_poplns/all_sites.BED",
    deviant_snps="results/bed/hap{hap}/deviant_sites/hap{hap}_combined_deviant_SNPs_realign_BH_correction.BED"
  output:
    filtered_sites_bed="results/bed/hap{hap}/canonical_sites/filtered_calcLR/calcLR_filtered_sites.BED",
    filtered_sites_txt="results/bed/hap{hap}/canonical_sites/filtered_calcLR/calcLR_filtered_sites.txt"
  log:
    calcLR_log="results/logs/bedtools/hap{hap}/canonical_sites/filtered_calcLR/calcLR_filtered_sites.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Filter out deviant sites for all populations using bedtools
    bedtools subtract -a {input.all_sites_bed} -b {input.deviant_snps} > {output.filtered_sites_bed} \
    2> {log.calcLR_log}

    # Convert the filtered BED file to a .txt file formatted for -sites in ANGSD
    awk '{{print $1, $3}}' {output.filtered_sites_bed} > {output.filtered_sites_txt}
    """


rule filter_all_sites_all_populations_dupHMM:
  input:
    filtered_calcLR_bed="results/bed/hap{hap}/canonical_sites/filtered_calcLR/calcLR_filtered_sites.BED",
    dupHMM_sites="results/bed/hap{hap}/deviant_sites/hap{hap}_combined_dupHMM_regions.BED"
  output:
    filtered_dupHMM_bed="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.BED",
    filtered_dupHMM_txt="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt"
  log:
    dupHMM_log="results/logs/bedtools/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_filtered_sites.log"
  envmodules:
    "bedtools/2.31.0"
  shell:
    """
    # Filter out deviant sites for all populations using bedtools
    bedtools subtract -a {input.filtered_calcLR_bed} -b {input.dupHMM_sites} > {output.filtered_dupHMM_bed} \
    2> {log.dupHMM_log}

    # Convert the filtered BED file to a .txt file formatted for -sites in ANGSD
    awk '{{print $1, $3, $3 + 1}}' {output.filtered_dupHMM_bed} > {output.filtered_dupHMM_txt}
    """


rule index_all_sites_all_popln_calcLR:
  input: 
    calcLR_sites="results/bed/hap{hap}/canonical_sites/filtered_calcLR/calcLR_filtered_sites.txt"
  output: 
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_calcLR/calcLR_filtered_sites.txt.bin",
    idx_index="results/bed/hap{hap}/canonical_sites/filtered_calcLR/calcLR_filtered_sites.txt.idx"
  envmodules:
    "angsd/0.940"
  shell: 
    """
    angsd sites index {input.calcLR_sites}
    """


rule index_all_sites_all_popln_dupHMM:
  input: 
    dupHMM_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt"
  output: 
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.bin",
    idx_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.idx"
  envmodules:
    "angsd/0.940"
  shell: 
    """
    angsd sites index {input.dupHMM_sites}
    """


# Estimate SAF, HWE, GL with SNPs on entire population
rule angsd_SNP_on_all_populations:
  input:
    bam_list="data/lists/hap{hap}/all_populations_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.bin",
    idx_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/dupHMM_calcLR_filtered_sites.txt.idx"
  output:
    arg_file="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.arg",
    mafs_file="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.mafs.gz",
    hwe_file="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.hwe.gz",
    depth_sample="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.depthSample",
    depth_global="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.depthGlobal",
    saf_1="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.saf.idx",
    saf_2="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.saf.pos.gz",
    saf_3="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.saf.gz",
    beagle="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.beagle.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt"
  log:
    "results/logs/angsd/hap{hap}/canonical/all_poplns/angsd_canonical_SNPs_all_poplns.log"
  envmodules:
    "angsd/0.940"
  threads: 16
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 2\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -minMapQ 30\
    -minQ 20\
    -minInd 100\
    -minHWEpval 0.01\
    -minMaf 0.0001\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -SNP_pval 1e-6\
    -doHWE 1\
    -doCounts 1\
    -doDepth 1\
    -doMajorMinor 1\
    -doMaf 1\
    -doSaf 1\
    -anc {params.ref}\
    -doGlf 2\
    &> {log}
    """



# Calculate PCA on genotype likelihoods using PCAngsd
rule PCAngsd_all_populations:
  input:
    beagle="results/angsd/hap{hap}/canonical/all_poplns/all_poplns_canonical_SNPs.beagle.gz"
  output:
    cov_matrix="results/pcangsd/hap{hap}/canonical/all_poplns/all_popln_canonical_SNP_PCAngsd.cov"
  params:
    file_name="results/pcangsd/hap{hap}/canonical/all_poplns/all_popln_canonical_SNP_PCAngsd"
  log:
    "results/logs/pcangsd/hap{hap}/canonical/all_poplns/all_popln_canonical_SNP_PCAngsd.log"
  threads: 8
  envmodules:
    "python/3.10.2"
  shell:
    """
    pcangsd -b {input.beagle}\
    -o {params.file_name}\
    -t {threads}\
    --iter 1000\
    2> {log}
    """


# Plot PCAs
rule PCAngsd_all_populations_plots:
  input:
    cov_matrix="results/pcangsd/hap{hap}/canonical/all_poplns/all_popln_canonical_SNP_PCAngsd.cov",
    pop_info="data/lists/Batch1_pop_names.info"
  output:
    plot="results/plots/hap{hap}/PCAngsd/all_popln_canonical_SNP_PCAngsd.png"
  threads: 2
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/plot_PCA.R {input.cov_matrix} {input.pop_info} {output.plot}"


# Calculate geographical distances between pairwise individuals between populations
rule PCA_calc_geo_distances:
  input:
    csv="data/lists/bamlist_hap2_geo_coord.csv"
  output:
    dist="results/pcangsd/hap2/canonical_realign/all_poplns/pairwise_individ_geodist.csv"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/calc_popln_geodist.R {input.csv} {output.dist}"


# ANALYSIS 4: Thetas (nucleotide diversity, etc)
# We do not filter for MAFs
rule angsd_for_thetas:
  input:
    bam_list="data/lists/hap{hap}/{population}_realign_hap{hap}.txt",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    bin_index="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt.bin",
    fasta_fai="data/reference/hap{hap}/lupinehap{hap}.fasta.fai"
  output:
    arg_file="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.arg",
    mafs_file="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.mafs.gz",
    saf_1="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.idx",
    saf_2="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.pos.gz",
    saf_3="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.gz",
    glf_file="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.glf.gz"
  params:
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta",
    file_name="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites",
    scaffolds="results/scaffolds/hap{hap}_scaffolds.txt"
  log:
    "results/logs/angsd/hap{hap}/canonical/thetas/by_popln/angsd_canonical_sites_hap{hap}_{population}.log"
  envmodules:
    "angsd/0.940"
  threads: 8
  shell:
    """
    angsd -bam {input.bam_list}\
    -ref {params.ref}\
    -out {params.file_name}\
    -remove_bads 1\
    -rf {params.scaffolds}\
    -GL 2\
    -C 50\
    -sites {input.canonical_sites}\
    -setMinDepth 25\
    -setMaxDepth 3500\
    -minMapQ 30\
    -minQ 20\
    -minInd 20\
    -baq 2\
    -only_proper_pairs 1\
    -nThreads {threads}\
    -doMajorMinor 1\
    -doMaf 1\
    -doGlf 1\
    -doSaf 1\
    -doCounts 1\
    -anc {params.ref}\
    &> {log}
    """
# NOTE: Did not filter with -minMAF 0.001

rule global_SFS_theta_by_population:
  input:
    saf_idx="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.idx",
    canonical_sites="results/bed/hap{hap}/canonical_sites/filtered_dupHMM/{population}_filtered_sites_dupHMM_calcLR.txt",
    ref="data/reference/hap{hap}/lupinehap{hap}.fasta"
  output:
    sfs="results/realSFS/hap{hap}/globalSFS/{population}_globalSFS_folded_theta.sfs"
  log:
    "results/logs/realSFS/hap{hap}/globalSFS/{population}_globalSFS_folded_theta.log"
  envmodules:
    "angsd/0.940"
  threads: 40
  shell:
    """
    realSFS {input.saf_idx}\
    -P {threads}\
    -seed 1\
    -fold 1\
    -anc {input.ref}\
    > {output.sfs}\
    2> {log}
    """


rule theta_prep_by_population:
  input:
    sfs="results/realSFS/hap{hap}/globalSFS/{population}_globalSFS_folded_theta.sfs",
    saf_idx="results/angsd/hap{hap}/canonical/thetas/by_popln/{population}_canonical_sites.saf.idx"
  output:
    theta="results/theta/hap{hap}/{population}_out.thetas.gz",
    index="results/theta/hap{hap}/{population}_out.thetas.idx"
  log:
    "results/logs/theta/hap{hap}/{population}_estimate_theta.log"
  params:
    file="results/theta/hap{hap}/{population}_out"
  envmodules:
    "angsd/0.940"
  threads: 4
  shell:
    """
    realSFS saf2theta {input.saf_idx}\
    -sfs {input.sfs}\
    -outname {params.file}\
    -fold 1\
    -P {threads}\
    &> {log}
    """


# Extract log scaled estimates of 
rule estimate_theta_by_sites:
  input:
    theta_index="results/theta/hap{hap}/{population}_out.thetas.idx"
  output:
    file_name="results/theta/hap{hap}/{population}_log_scale.out"
  envmodules:
    "angsd/0.940"
  threads: 4
  shell:
    """
    thetaStat print {input.theta_index}\
    > {output.file_name}
    """

rule estimate_theta_sliding_window:
  input:
    theta_index="results/theta/hap{hap}/{population}_out.thetas.idx"
  output:
    window="results/theta/hap{hap}/{population}_out.thetasWindow.gz.pestPG"
  params:
    file_name="results/theta/hap{hap}/{population}_out.thetasWindow.gz"
  log:
    "results/logs/theta/hap{hap}/{population}_estimate_theta_sliding_window.log"
  envmodules:
    "angsd/0.940"
  threads: 4
  shell:
    """
    thetaStat do_stat {input.theta_index}\
    -win 10000\
    -step 1000\
    -outnames {params.file_name}\
    &> {log}
    """

rule plot_thetas_output:
  input:
    theta_files=expand("results/theta/hap2/{population}_out.thetasWindow.gz.pestPG", population=POPULATIONS)
  output:
    tajimasD_plot="results/plots/hap2/tajimasD_boxplot.png",
    nucleotide_div_plot="results/plots/hap2/nucleotide_diversity_boxplot.png"
  log:
    "results/logs/theta/hap2/plot_theta_metrics.log"
  envmodules:
    "r/4.4.0"
  shell:
    """
    Rscript scripts/plot_theta.R
    Rscript scripts/plot_theta_by_scaffold.R
    """


# ANALYSIS 5: Pairwise Fst Between Populations
# We re-use estimated .saf files from the theta analysis

rule fst_prep_by_population:
  input:
    pop1_saf_idx="results/angsd/hap2/canonical/thetas/by_popln/{pop1}_canonical_sites.saf.idx",
    pop2_saf_idx="results/angsd/hap2/canonical/thetas/by_popln/{pop2}_canonical_sites.saf.idx"
  output:
    sfs_prior="results/realSFS/hap2/fst/{pop1}_{pop2}_prior.ml"
  log:
    "results/logs/realSFS/hap2/globalSFS_Fst/{pop1}_{pop2}_prior.log"
  envmodules:
    "angsd/0.940"
  threads: 40
  shell:
    """
    realSFS {input.pop1_saf_idx} {input.pop2_saf_idx} \
      -P {threads} \
      -fold 1 \
      > {output.sfs_prior} \
      2> {log}
    """

rule fst_analysis:
  input:
    pop1_saf_idx="results/angsd/hap2/canonical/thetas/by_popln/{pop1}_canonical_sites.saf.idx",
    pop2_saf_idx="results/angsd/hap2/canonical/thetas/by_popln/{pop2}_canonical_sites.saf.idx",
    sfs_prior="results/realSFS/hap2/fst/{pop1}_{pop2}_prior.ml"
  output:
    fst_idx="results/realSFS/hap2/fst/{pop1}_{pop2}.fst.idx"
  params:
    fst_out_prefix="results/realSFS/hap2/fst/{pop1}_{pop2}"
  log:
    "results/logs/realSFS/hap2/globalSFS_Fst/{pop1}_{pop2}_Fst_estimation.log"
  envmodules:
    "angsd/0.940"
  threads: 40
  shell:
    """
    realSFS fst index {input.pop1_saf_idx} {input.pop2_saf_idx} \
      -sfs {input.sfs_prior} \
      -fold 1 \
      -fstout {params.fst_out_prefix} \
      -P {threads} \
      2> {log}
    """

rule estimate_fst_stats:
  input:
    fst_idx="results/realSFS/hap2/fst/{pop1}_{pop2}.fst.idx"
  output:
    global_fst="results/realSFS/hap2/fst/{pop1}_{pop2}_fst_global.txt"
  params:
    window_size=50000,  # Window size for sliding window Fst
    step_size=10000  # Step size for sliding window Fst
  log:
    "results/logs/realSFS/hap2/globalSFS_Fst/{pop1}_{pop2}_Fst_extract.log"
  envmodules:
    "angsd/0.940"
  shell:
    """
    # Global Fst estimate
    realSFS fst stats {input.fst_idx} > {output.global_fst}
    """

    #window_fst="results/realSFS/hap2/fst/{pop1}_{pop2}_fst_windows.txt"
    # Fst in sliding windows
    # realSFS fst stats2 {input.fst_idx} -win {params.window_size} -step {params.step_size} > {output.window_fst}

rule plot_fst_by_distance:
  input:
    global_fst=expand("results/realSFS/hap2/fst/{pop1}_{pop2}_fst_global.txt")
  output:
    fst_plot="results/plots/hap2/Fst/{pop1}_{pop2}_fst_by_distance.png"
  log:
    "results/logs/Fst/{pop1}_{pop2}_fst_by_distance.log"
  envmodules:
    "r/4.4.0"
  shell:
    "Rscript scripts/
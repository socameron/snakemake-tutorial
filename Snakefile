###########################
##    WILDCARDS SET-UP   ##
###########################

# Python based 
import os
import glob
import subprocess
import pandas as pd
import yaml
import re
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


###########################
## BEGINNING OF WORKFLOW ##
###########################


# Expand the final files using the updated configuration
rule all:
  input:
    expand("results/bam_merge/hap{hap}/merged_hap{hap}.bam.bai", hap=haps)
    

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

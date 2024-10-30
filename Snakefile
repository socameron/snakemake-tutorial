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

# Rule to update the configuration file - you can use the config.yaml file to store sample names. 
#I incorporated a custom python script to extract this information
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
  


# snakemake-tutorial for DRAC Users (or Canadian Researchers)
Hi there! This tutorial is tailored towards Canadian researchers who have access to the [Digital Research Alliance of Canada](alliancecan.ca/en) advanced computing clusters. All credit is owed to [Eric Anderson](https://github.com/eriqande/mega-lcwgs-pw-fst-snakeflow) (Research Geneticist, NOAA) whom I met at [ConGen2023](https://www.umt.edu/ces/conferences/congen/). For more information about Snakemake's amazing abilities, see Eric Anderson's lecture slides [here](https://eriqande.github.io/con-gen-2023/slides/snake-slides.html#/section). 

**What is Snakemake?**
[Snakemake](https://snakemake.readthedocs.io/en/stable/) is a workflow manager that allows you to parallelize multiple jobs at once, while submittin new jobs dynamically as older ones are completed.  [NextFlow](https://www.nextflow.io/) is another workflow manager but is Java based. Each manager has varying levels of capabilities and documentation, so choose your pick. Either way, both are better than working with bash scripts!

Setting up Snakemake on DRAC clusters is simpler than you think:

1. First load Python and create and a virtual environment.

```
# We create a virtual environment in the login node.
module load python/3.12.4
virtualenv --no-download ENV

# Just labelled the virtual environment folder as 'ENV'
# Activate virtual environment

source ENV/bin/activate
```

2. Install and upgrade pip (python installer) plus Snakemake

```
pip install --no-index --upgrade pip
pip install --no-index snakemake
pip install snakemake-executor-plugin-cluster-generic
pip install snakemake-executor-plugin-slurm

# Double check version of Snakemake  - this is compatible with v8.16.0
snakemake --version

```
** NOTE ** : The newest version of Snakemake has released 'plug-ins'. These are now required to interact with SLURM and config.yaml files. 

3. Download Snakemake files from my repository

Using Snakemake on DRAC clusters requires 4 items. (1) The `Snakemake` file; (2) a snakeprofile folder containing (3) a `config.yaml` and (4) `status.sacct-robust.sh` file. The Snakemake file is where all rules are written for your workflow. The `config.yaml` specifies the resources allocated for each rule (although this can also be written into your Snakemake file) plus other system settings. The `status-sacct-robust.sh` file interacts Snakemake with the SLURM scheduler to send jobs via `sbatch`. 

```
git clone https://github.com/socameron/snakemake-tutorial.git
```

4. Check that data files have readable permissions

Wherever your storing your data, say in `scratch/data`, then use `chmod ug+x *.fastq.gz` within that folder (logged onto your cluster, of course) to change permission rules. 

5. Call Snakemake using the profile set up

```
snakemake --profile snakeprofile --executor cluster-generic all
# Calls snakemake with--profile settings from snakeprofile. Tells Snakemake to specifically run the rule 'all' rather than calling for a specific file.
# Add -np to run a practice run
# --executor cluster-generic is required on snakemake version greater than 8.0
```

6. Call Snakemake in a SLURM job ** **RECOMMENDED**

Calling Snakemake (i.e building DAG of jobs) on a login-node may take extremely long depending on your workflow complexity. To decrease timing, DRAC suggests building a virtual environment within a job and activating Snakemake in the job. See https://docs.alliancecan.ca/wiki/Python#Creating_virtual_environments_inside_of_your_jobs for details. This should reduce wait times for 'Building DAG of jobs' and require no screen (if you wish to log off your cluster).

Here's an example of an SBATCH script that you could use for running snakemake from inside a job. Be sure to set the total time to equal or greater than the biggest snakemake rule required.

```
#!/bin/bash
#SBATCH --account= XXX
#SBATCH --mem-per-cpu=10G      # decrease/increase as needed
#SBATCH --time=8:00:00         # equal or greater than longest snakemake rule
#SBATCH --mail-user= yourusername@email.com
#SBATCH --output=/path/to/some/folder/where/snakemake/outputs/are/held/snakemake_%j.out # %j is the job number

module load python/3.10
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index snakemake

snakemake --profile snakeprofile all
```

Then you get a constant update of the snakemake in your terminal, for last 50 lines:

```
cd /path/to/some/folder/where/snakemake/outputs/are/held/snakemake_%j.out
tail -f SNAKEMAKE_%j.out file -n 50 
```

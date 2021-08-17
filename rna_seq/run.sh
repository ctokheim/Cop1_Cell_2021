#!/bin/bash
#SBATCH --job-name=COP1_RIMA
#SBATCH --mem=200G       # total memory need
#SBATCH --time=96:00:00
#SBATCH -c 32 #number of cores

export CONDA_ROOT=/liulab/linyang/rnaseq_env/miniconda3
export PATH=/liulab/linyang/rnaseq_env/miniconda3/bin:$PATH
source activate /liulab/linyang/rnaseq_env/miniconda3/envs/rnaseq

snakemake -s rnaseq.snakefile -j 4 --latency-wait 30

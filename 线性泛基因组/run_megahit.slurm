#!/bin/bash


#SBATCH --job-name="megahit"
#SBATCH -n 2 #threads
#SBATCH -N 1 #node number
#SBATCH --mem=2000
#SBATCH --partition=cuPartition
#SBATCH --mail-user=mingyan24@126.com
#SBATCH --mail-type=BEGIN,END,FAIL

source activate map2pan
snakemake --cluster "sbatch --output=/data/myan/slurm.out/%j.out \
--error=/data/myan/slurm.out/%j.out --cpus-per-task={threads} \
--mail-type=END,FAIL --mail-user=mingyan24@126.com --mem={resources.mem} --partition=tcum256c128Partition" \
--jobs 8 -s megahit.smk

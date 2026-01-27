#!/usr/bin/env bash

#SBATCH --job-name=cNMF_factorize
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=36G
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --partition=cpu2023,cpu2022-bf24,cpu2022,cpu2021-bf24,cpu2021,cpu2019-bf05,cpu2019,cpu2017-bf05,single,sherlock
#SBATCH --time=05:00:00
#SBATCH -o %A_%a.stdout.txt
#SBATCH -e %A_%a.stderr.txt
#SBATCH --export=ALL
#SBATCH --array=0-499

####### Set environment variables ###############
source ~/software/init-conda
conda activate mosaic

####### Run your script #########################
echo $1 "/" $2

mosaicmpi factorize --output_dir $1 --name $2 --total_workers 500 --worker_index ${SLURM_ARRAY_TASK_ID}

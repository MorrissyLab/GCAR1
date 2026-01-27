#!/usr/bin/env bash

#SBATCH --job-name=cNMF_GCAR1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --partition=cpu2023,cpu2022-bf24,cpu2022,cpu2021-bf24,cpu2021,cpu2019-bf05,cpu2019,cpu2017-bf05,single,sherlock,bigmem
#SBATCH --time=24:00:00
#SBATCH -o %j.stdout.txt
#SBATCH -e %j.stderr.txt
#SBATCH --export=ALL

echo $PWD

####### Set environment variables ###############
# Source and activate the mosaicMPI conda enviroment
source ~/software/init-conda
conda activate mosaic

RUN_NAME_VISHD="cNMF_5_30_5_ST"
COUNT_FILE_VISHD="primary_biopsy_DS_RA_M_count_filtered_024um.tsv"
METADATA_FILE_VISHD="primary_biopsy_DS_RA_M_metadata_24um.tsv"
OD_GENES="primary_biopsy_DS_RA_M_ODgenes_list_024um.txt"

# To have a single rank, set start == end and by = 1
RANK_START="5"
RANK_END="30"
RANK_BY="5"
RUN_RANGE=($(seq ${RANK_START} ${RANK_BY} ${RANK_END}))

BETA_LOSS="kullback-leibler"

SAMPLES="samples.csv"
SR_OUTPUT="/work/morrissy_lab/gurveer/gpnmb_car/datasets/"

####### Run the CNMF script #########################
# Prepare the data for MosaicMPI sfactorization.
mosaicmpi txt-to-h5ad --data_file ${COUNT_FILE_VISHD} --metadata ${METADATA_FILE_VISHD} -o dataset_vishd.h5ad
mosaicmpi check-h5ad -i dataset_vishd.h5ad -o dataset_vishd_filtered.h5ad
mosaicmpi model-odg --name ${RUN_NAME_VISHD} --input dataset_vishd_filtered.h5ad
mosaicmpi set-parameters --name ${RUN_NAME_VISHD} --k_range ${RANK_START} ${RANK_END} ${RANK_BY} --beta_loss ${BETA_LOSS} -m genes_file -p ${OD_GENES}

# Submit the factorization script using its own sh script, so that we can utilize a Slurm array to speed it up.
sbatch cnmf_factorize.sh $PWD ${RUN_NAME_VISHD}

# Postprocess the factorization step.
mosaicmpi postprocess --output_dir $PWD --name ${RUN_NAME_VISHD} --cpus ${SLURM_JOB_CPUS_PER_NODE}

# Create the usage heatmap plots
mosaicmpi annotated-heatmap --output_dir ${RUN_NAME_VISHD} -i ${RUN_NAME_VISHD}/${RUN_NAME_VISHD}.h5ad

####### Create GEP plots #########################
# Create gene expression spatial plots from the data in a plots folder in the RUN_NAME folder
echo ${#RUN_RANGE[@]}

METADATA_FILE_CSV=$(echo $METADATA_FILE_VISHD | cut -d'.' -f1)'.csv'
awk 'BEGIN { FS="\t"; OFS="," } {$1=$1; print}' $METADATA_FILE_VISHD > $METADATA_FILE_CSV
sed -i 's/\"//g' $METADATA_FILE_CSV

# Make the neccessary directories, but ideally this should be apart of the plots.py script
PLOT_DIRECTORY="plots_"${RUN_NAME_VISHD}

HIRES_DIRECTORY=${PLOT_DIRECTORY}"/hires_vectors/gene_expression"
LOWRES_DIRECTORY=${PLOT_DIRECTORY}"/lowres_jpeg/gene_expression"

COMBINED_DIRECTORY=${PLOT_DIRECTORY}"/combined"
COMBINED_HIRES_DIRECTORY=${COMBINED_DIRECTORY}"/hires_vectors/gene_expression/binarized"
COMBINED_LOWRES_DIRECTORY=${COMBINED_DIRECTORY}"/lowres_jpeg/gene_expression/binarized"

mkdir -p ${PLOT_DIRECTORY}
mkdir -p ${HIRES_DIRECTORY}
mkdir -p ${LOWRES_DIRECTORY}
mkdir -p ${COMBINED_DIRECTORY}
mkdir -p ${COMBINED_HIRES_DIRECTORY}
mkdir -p ${COMBINED_LOWRES_DIRECTORY}

sbatch --array=0-$(( ${#RUN_RANGE[@]} - 1)) cnmf_usage_plots.sh ${RUN_NAME_VISHD} ${PLOT_DIRECTORY} ${METADATA_FILE_CSV} ${SAMPLES} ${SR_OUTPUT} ${RUN_RANGE[@]}
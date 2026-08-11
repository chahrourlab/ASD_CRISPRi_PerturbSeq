#!/bin/bash
#SBATCH --job-name factorize cNMF
#SBATCH -N 1
#SBATCH -p 256GBv2
#SBATCH -t 7-0:0:0
#SBATCH -o ./cnmf_factor_%j.out
#SBATCH -e ./cnmf_factor_%j.err

module load load_rhel7_apps
module load python/latest-3.12.x-anaconda  
eval "$(conda shell.bash hook)"
conda activate cnmf_env
date

echo "factorize. Total worker: $1"
echo "factorize. worker index: $2"
cnmf factorize --output-dir ./ --name Neuron_Prod1_cNMF_run --total-workers $1 --worker-index $2

echo "finish"
date

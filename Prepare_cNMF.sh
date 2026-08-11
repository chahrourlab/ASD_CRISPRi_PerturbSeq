#!/bin/bash
#SBATCH --job-name prepare cNMF
#SBATCH -N 1
#SBATCH -p 512GB 
#SBATCH -t 7-0:0:0
#SBATCH -o cnmf_%j.out
#SBATCH -e cnmf_%j.err

TOTAL_WORKERS=15
K_LIST=(50 100 200 250 500)
NUM_ITER=20
K_TEXT=""
for k in ${K_LIST[@]}
do
    K_TEXT="$K_TEXT$k "
done
echo "Total workers: $TOTAL_WORKERS"
echo "k list: $K_TEXT"
echo "number of iteration: $NUM_ITER"

module load load_rhel7_apps
module load python/latest-3.12.x-anaconda  
eval "$(conda shell.bash hook)"
conda activate cnmf_env

echo "prepare"
date``
cnmf prepare --output-dir ./ --name Neuron_Prod1_cNMF_run -c raw_count.h5ad -k ${K_TEXT} --n-iter ${NUM_ITER} --total-workers ${TOTAL_WORKERS} --seed 14 --numgenes 2000 --beta-loss frobenius


for num in `seq 0 $((TOTAL_WORKERS-1))`
do
    sbatch factorize_cNMF.sh ${TOTAL_WORKERS} ${num}
done

echo "finish"

date

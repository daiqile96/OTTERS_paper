#!/bin/bash

# GET JOB INFO FROM LIST OF TRAIN JOBS 
i=${SGE_TASK_ID}

module load R/4.0.4

train_jobs_path=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/2_simulate_GWAS_Z_allN.txt

line=`head -n $i ${train_jobs_path} | tail -n1`
he2=`echo $line | cut -d' ' -f 1`
cp=`echo $line | cut -d' ' -f 2`
N=`echo $line | cut -d' ' -f 3`

Rscript /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/4_evaluate_Power.R ${he2} ${cp} ${N}

# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N TWAS_power -pe smp 2 -t 1-32 -tc 16 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/4_evaluate_Power.sh

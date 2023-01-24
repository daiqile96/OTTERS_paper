#!/bin/bash

# GET JOB INFO FROM LIST OF TRAIN JOBS 
i=${SGE_TASK_ID}

module load R/4.0.4

wk_dir=/home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu
train_jobs_path=${wk_dir}/4_evaluate_R2.txt

# for he2 in 0.01 0.05 0.1 ; do  
# 	for cp in 0.001 0.01 ; do
#         echo -e "${he2}\t${cp}" >> ${train_jobs_path}
#     done 
# done

n_row=`wc -l ${train_jobs_path} | cut -d' ' -f 1`
line=`head -n $i ${train_jobs_path} | tail -n 1`
he2=`echo ${line} | cut -d' ' -f 1`
cp=`echo ${line} | cut -d' ' -f 2`

Rscript ${wk_dir}/4_evaluate_R2.R ${he2} ${cp}

# qsub -q b.q -cwd -j y -wd /home/qdai/YangFSSdata2/qdai/new_logs -N Simu_R2 -pe smp 2 -t 1-6 -tc 6 /home/qdai/YangFSSdata2/qdai/BS_TWAS/Scripts/4_revise_simu/4_evaluate_R2.sh

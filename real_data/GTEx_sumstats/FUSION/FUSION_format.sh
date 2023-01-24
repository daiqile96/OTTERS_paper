#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=10:00:00
#SBATCH --partition=yanglab,preemptable,week-long-cpu,month-long-cpu,largemem
#SBATCH --job-name=FUSION_format
#SBATCH --error=/home/qdai8/projects/logs/job.%A_%a.error
#SBATCH --output=/home/qdai8/projects/logs/job.%A_%a.out

set=${SLURM_ARRAY_TASK_ID}
train=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/train_group.txt
line=`head -n $set ${train} | tail -n 1`
chr=`echo ${line} | cut -d' ' -f 1`
batch=`echo ${line} | cut -d' ' -f 2`

module load R 
Rscript ~/projects/OTTERS/Scripts/3_ReRDA/FUSION_format.R ${chr} ${batch}
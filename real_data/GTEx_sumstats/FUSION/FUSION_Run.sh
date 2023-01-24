#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=30:00:00
#SBATCH --partition=yanglab,preemptable,largemem
#SBATCH --job-name=FUSION
#SBATCH --error=/home/qdai8/projects/logs/job.%A_%a.error
#SBATCH --output=/home/qdai8/projects/logs/job.%A_%a.out

module load R

chr=${SLURM_ARRAY_TASK_ID}
exp_anno=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/GTEx_CHR${chr}_GeneAnno.txt
n_row=`wc -l ${exp_anno} | cut -d' ' -f 1`

### Run FUSION
for gene in $(seq 2 $n_row)
do

line=`head -n $gene ${exp_anno} | tail -n 1`
ID=`echo ${line} | cut -d' ' -f 4`

cd /home/qdai8/projects/OTTERS/Data/ReRDA/ReFUSION/CHR${chr}/${ID}

b_dir=${ID}
tmp_dir=${ID}.tmp
o_dir=${ID}

Rscript /home/qdai8/projects/OTTERS/Tool/fusion_twas-master/FUSION.compute_weights.R \
--bfile $b_dir \
--tmp $tmp_dir \
--out $o_dir \
--hsq_p 1 \
--PATH_plink /home/qdai8/projects/OTTERS/Tool/plink \
--PATH_gcta /home/qdai8/projects/OTTERS/Tool/fusion_twas-master/gcta_nr_robust \
--PATH_gemma /home/qdai8/projects/OTTERS/Tool/gemma-linux \
--models top1,lasso,enet,blup

done 
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-23:00:00
#SBATCH --partition=yanglab,preemptable,week-long-cpu,month-long-cpu,largemem
#SBATCH --job-name=RDA_II
#SBATCH --error=/home/qdai8/projects/logs/job.%A_%a.error
#SBATCH --output=/home/qdai8/projects/logs/job.%A_%a.out

set=${SLURM_ARRAY_TASK_ID}
train=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/train_group.txt
line=`head -n $set ${train} | tail -n 1`
chr=`echo ${line} | cut -d' ' -f 1`
batch=`echo ${line} | cut -d' ' -f 2`

clump=0.99
source ~/.bashrc
conda activate otters
module load R
N_THREADS=1
N_TASK=4

# prevent automatically using  all available cores on a compute node
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

######### Run OTTERS to prepare inputs for imputation models ############
### Set Tool Directory
OTTERS_DIR=/home/qdai8/projects/tool/OTTERS

### Set Required Inputs
exp_anno=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/Batch${batch}_GTEx_CHR${chr}_GeneAnno.txt
geno_dir=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Binary/CHR${chr}/CHR${chr}
eQTL_dir=/home/qdai8/projects/OTTERS/Data/ReRDA/Results/CHR${chr}/Batch${batch}
twas_dir=/home/qdai8/projects/OTTERS/Data/ReRDA/Results/CHR${chr}/Batch${batch}/TWAS
gwas_file=/home/qdai8/projects/Data/UKBB_GWAS/CARDIO_Hg38_Zscore.txt

python3 ${OTTERS_DIR}/testing.py \
--OTTERS_dir=${OTTERS_DIR} \
--weight_dir=${eQTL_dir} \
--models=FUSION \
--anno_dir=${exp_anno} \
--geno_dir=${geno_dir} \
--out_dir=${twas_dir} \
--gwas_file=${gwas_file} \
--chrom=${chr} \
--thread=${N_TASK}
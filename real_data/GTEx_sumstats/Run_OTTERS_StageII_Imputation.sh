#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition=yanglab,preemptable,week-long-cpu,month-long-cpu,largemem
#SBATCH --job-name=ReRDA_II
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
grex_dir=/home/qdai8/projects/OTTERS/Data/ReRDA/Results/CHR${chr}/Batch${batch}/GReX
#sample_dir=/home/qdai8/projects/OTTERS/Data/RDA/sampleID.txt
sample_dir=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/SampleID.txt

python3 ${OTTERS_DIR}/imputing.py \
--OTTERS_dir=${OTTERS_DIR} \
--weight_dir=${eQTL_dir} \
--models=SDPR,PRScs,lassosum,P0.001,P0.05,FUSION \
--anno_dir=${exp_anno} \
--geno_dir=${geno_dir} \
--out_dir=${grex_dir} \
--chrom=${chr} \
--samples=${sample_dir} \
--thread=${N_TASK}



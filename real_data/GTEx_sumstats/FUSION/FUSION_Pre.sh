#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition=yanglab,preemptable,week-long-cpu,month-long-cpu
#SBATCH --job-name=FUSION_pre
#SBATCH --error=/home/qdai8/projects/logs/job.%A_%a.error
#SBATCH --output=/home/qdai8/projects/logs/job.%A_%a.out

chr=${SLURM_ARRAY_TASK_ID}
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
exp_anno=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/GTEx_CHR${chr}_GeneAnno.txt
geno_dir=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Binary/CHR${chr}/CHR${chr}
out_dir=/home/qdai8/projects/OTTERS/Data/ReRDA/ReFUSION/CHR${chr}
sampleID=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/SampleID.txt
mkdir -p ${out_dir}

python3 ${OTTERS_DIR}/FUSION_Pre.py \
--anno_dir=${exp_anno} \
--geno_dir=${geno_dir} \
--out_dir=${out_dir} \
--sampleID_dir=${sampleID} \
--chrom=${chr} \
--thread=${N_TASK}

### Edit .fam to add gene expression data
module load R
Rscript /home/qdai8/projects/OTTERS/Scripts/3_ReRDA/FUSION_Add_GeneExp.R ${chr}





#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=5-00:00:00
#SBATCH --partition=yanglab,preemptable,week-long-cpu,month-long-cpu,largemem
#SBATCH --job-name=RDA_I
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

######### Run OTTERS to prepare inputs for imputation models ############
### Set Tool Directory
OTTERS_DIR=/home/qdai8/projects/tool/OTTERS
SDPR_DIR=/home/qdai8/projects/bin/SDPR
# make sure the dynamic libraries are not changed
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib
# prevent automatically using  all available cores on a compute node
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

# ######### Convert VCF genotype files to PLINK 1 binary files ############
# VCF_dir=/projects/YangLabData/SharedData/GTExV8/GenotypeFiles/
# BIN_dir=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Binary

# VCF_file=${VCF_dir}/CHR${chr}_GTEx_WGS.vcf.gz
# BIN_prefix=${BIN_dir}/CHR${chr}/CHR${chr}
# mkdir -p ${BIN_dir}/CHR${chr}
# plink --vcf ${VCF_file} --keep-allele-order --make-bed --maf 0.01 --out ${BIN_prefix}

### Set Required Inputs
exp_anno=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/Batch${batch}_GTEx_CHR${chr}_GeneAnno.txt
geno_dir=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Binary/CHR${chr}/CHR${chr}
sst_file=/home/qdai8/projects/Data/eQTLGen/eQTLGen_Hg38.txt
output_dir=/home/qdai8/projects/OTTERS/Data/RDA/Results/CHR${chr}/Batch${batch}

######### Run Imputation Models ################
python3 ${OTTERS_DIR}/training.py \
--OTTERS_dir=${OTTERS_DIR} \
--SDPR_dir=${SDPR_DIR} \
--anno_dir=${exp_anno} \
--geno_dir=${geno_dir} \
--sst_file=${sst_file} \
--out_dir=${output_dir} \
--chrom=${chr} \
--r2=${clump} \
--models=lassosum,PT,PRScs,SDPR \
--thread=${N_TASK}



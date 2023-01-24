
# calculate Test R2 in 315 GTEx v8 samples
module load R
script_dir=/home/qdai8/projects/OTTERS/Scripts/1_RDA/
Rscript ${script_dir}/Summary_OTTERS_GReX.R

# merge TWAS results
res_dir=/home/qdai8/projects/OTTERS/Data/RDA/Results
twas_res_dir=${res_dir}/TWAS
mkdir -p ${twas_res_dir}

# create file for each chromosome 
for chr in {1..22}
do
   mkdir ${twas_res_dir}/CHR${chr}
   for method in lassosum P0.001 P0.05 PRScs SDPR 
   do 
      out_twas=${twas_res_dir}/CHR${chr}/${method}.txt
      echo -e 'CHROM\tGeneStart\tGeneEnd\tTargetID\tn_snps\tFUSION_Z\tFUSION_PVAL' > ${out_twas}
   done
done

# merge results from different batches
train=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno/train_group.txt
for set in {1..207}
do
   line=`head -n $set ${train} | tail -n 1`
   chr=`echo ${line} | cut -d' ' -f 1`
   batch=`echo ${line} | cut -d' ' -f 2`
   
   batch_twas_dir=${res_dir}/CHR${chr}/Batch${batch}/TWAS

   for method in lassosum P0.001 P0.05 PRScs SDPR 
   do 
      out_twas=${twas_res_dir}/CHR${chr}/${method}.txt
      tail -n +2 ${batch_twas_dir}/${method}.txt >> ${out_twas}
   done
done

# Perform ACAT
GReX_r2_dir=${res_dir}/GReX_R2.txt
anno_dir=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Anno
methods=SDPR,PRScs,lassosum,P0.001,P0.05
filter=T
Rscript ${script_dir}/Summary_OTTERS_TWAS.R ${twas_res_dir} ${GReX_r2_dir} ${anno_dir} ${methods} ${filter}

filter=F
Rscript ${script_dir}/Summary_OTTERS_TWAS.R ${twas_res_dir} ${GReX_r2_dir} ${anno_dir} ${methods} ${filter}
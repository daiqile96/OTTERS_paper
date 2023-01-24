
######### Convert VCF genotype files to PLINK 1 binary files ############
VCF_dir=/projects/YangLabData/SharedData/GTExV8/GenotypeFiles/
BIN_dir=/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V8/Binary

VCF_file=${VCF_dir}/CHR${chr}_GTEx_WGS.vcf.gz
BIN_prefix=${BIN_dir}/CHR${chr}/CHR${chr}
mkdir -p ${BIN_dir}/CHR${chr}
plink --vcf ${VCF_file} --keep-allele-order --make-bed --maf 0.01 --out ${BIN_prefix}

######### Liftover and format eQTLGen summary statistics ###############

### download eQTLGen summary statistics ###
eQTL_dir=/home/qdai8/projects/Data/eQTLGen
wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz -P ${eQTL_dir}

### Format the summary statistics to perform liftover ###
eQTLraw_file=${eQTL_dir}/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
unlift_file=${eQTL_dir}/eQTLGen_unlift.txt
lifted_file=${eQTL_dir}/eQTLGen_Hg38.txt

# create the header of the new file (unlifted summary statistics)
echo -e "#CHROM\tPOS\tA1\tA2\tZ\tTargetID\tN" > ${unlift_file}
# select chrom, pos, A1, A2, Zscore, TargetID, N from the raw eQTL summary statistics
# also filter SNPs with sample size N < 3000
zcat ${eQTLraw_file} | \
tail -n +2 | \
awk -F "\t" '$13 > 3000' | \
awk 'BEGIN { OFS = "\t" } { print $3, $4, $5, $6, $7, $8, $13}' >> ${unlift_file}

# use liftover tool to do liftover
lift_over_tool=/home/qdai8/projects/bin/lift-over-vcf/
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -P ${lift_over_tool}
chmod a+x ${lift_over_tool}/liftOver

${lift_over_tool}/awk_lift.sh -chain hg19ToHg38 -in ${unlift_file} -out ${lifted_file} -liftOver_dir ${lift_over_tool}

### BGZIP the lifted eQTLgen summary statistics ###
echo bgziping.
bgzip -f ${lifted_file}
### TABIX the lifted eQTLGen summary statistics ###
echo tabixing.
tabix -f -b2 -e2 -S1 ${lifted_file}.gz

######### Prepare GWAS summary statistics ###############

### download GWAS summary statistics ###
gwas_dir=/home/qdai8/projects/Data/UKBB_GWAS
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/disease_CARDIOVASCULAR.sumstats.gz -P ${gwas_dir} 

### Format the GWAS summary statistics to perform liftover ###
gwasraw_file=${gwas_dir}/disease_CARDIOVASCULAR.sumstats.gz
unlift_file=${gwas_dir}/CARDIO_unlift.txt
lifted_file=${gwas_dir}/CARDIO_Hg38.txt

# create the header of the new file (unlifted summary statistics)
echo -e "#CHROM\tPOS\tA1\tA2\tBeta\tP" > ${unlift_file}
# select chrom, pos, A1, A2, Beta, P from the raw eQTL summary statistics
# NOTE here, in this GWAS summary statistics, the A1 is the reference allele (which is A2 in our definition)
zcat ${gwasraw_file} | \
tail -n +2 | \
awk 'BEGIN { OFS = "\t" } { print $2, $3, $5, $4, $8, $10}' >> ${unlift_file}

${lift_over_tool}/awk_lift.sh -chain hg19ToHg38 -in ${unlift_file} -out ${lifted_file} -liftOver_dir ${lift_over_tool}

### Sorting ###
echo bgziping.
bgzip -f ${lifted_file}

### Convert p-values in the summary statistics to Z-scores ###
convert_tool=/home/qdai8/projects/OTTERS/tools/convert_pvalue_zscore.py
conda activate otters
python3 ${convert_tool} --p_dir ${lifted_file}.gz --out_dir ${gwas_dir}/CARDIO_Hg38_Zscore.txt


### BGZIP the converted file ###
echo bgziping.
bgzip -f ${gwas_dir}/CARDIO_Hg38_Zscore.txt
echo tabixing.
tabix -f -b2 -e2 -S1 ${gwas_dir}/CARDIO_Hg38_Zscore.txt.gz

######### Prepare GWAS summary statistics ###############

### download GWAS summary statistics ###
gwas_dir=/home/qdai8/projects/Data/UKBB_GWAS
wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB/body_BMIz.sumstats.gz -P ${gwas_dir} 

### Format the GWAS summary statistics to perform liftover ###
gwasraw_file=${gwas_dir}/body_BMIz.sumstats.gz
unlift_file=${gwas_dir}/BMI_unlift.txt
lifted_file=${gwas_dir}/BMI_Hg38.txt

# create the header of the new file (unlifted summary statistics)
echo -e "#CHROM\tPOS\tA1\tA2\tBeta\tP" > ${unlift_file}
# select chrom, pos, A1, A2, Beta, P from the raw eQTL summary statistics
# NOTE here, in this GWAS summary statistics, the A1 is the reference allele (which is A2 in our definition)
zcat ${gwasraw_file} | \
tail -n +2 | \
awk 'BEGIN { OFS = "\t" } { print $2, $3, $5, $4, $8, $10}' >> ${unlift_file}

${lift_over_tool}/awk_lift.sh -chain hg19ToHg38 -in ${unlift_file} -out ${lifted_file} -liftOver_dir ${lift_over_tool}

### Sorting ###
echo bgziping.
bgzip -f ${lifted_file}

### Convert p-values in the summary statistics to Z-scores ###
convert_tool=/home/qdai8/projects/OTTERS/tools/convert_pvalue_zscore.py
conda activate otters
python3 ${convert_tool} --p_dir ${lifted_file}.gz --out_dir ${gwas_dir}/BMI_Hg38_Zscore.txt


# ### BGZIP the converted file ###
# echo bgziping.
# bgzip -f ${gwas_dir}/CARDIO_Hg38_Zscore.txt
# echo tabixing.
# tabix -f -b2 -e2 -S1 ${gwas_dir}/CARDIO_Hg38_Zscore.txt.gz


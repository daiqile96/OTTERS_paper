# This script prepare the genotype data and SNP information for the simulation study
# snp information include chrom, pos, A1, A2, maf
# Reference: http://alimanfoo.github.io/2017/06/14/read-vcf.html

import sys
import allel
import numpy as np
import pandas as pd
import getopt


def parse_param():
    long_opts_list = ['vcf_in=', 'geno_out=', 'help']

    param_dict = {'vcf_in': None, 'geno_out': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)
        except:
            print('Option not recognized.')
            print('Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--vcf_in":
                param_dict['vcf_in'] = arg
            elif opt == "--geno_out":
                param_dict['geno_out'] = arg
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['vcf_in'] is None:
        print('* Please specify the path of vcf file using --vcf_in\n')
        sys.exit(2)
    elif param_dict['geno_out'] is None:
        print('* Please specify the output path using --out_dir\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict

############################################################

param_dict = parse_param()

# read in vcf files
myvcf = param_dict['vcf_in']
out_geno_path = param_dict['geno_out']

# Some fields like ‘ALT’ can have a variable number of values.
# I.e., each variant may have a different number of data values for these fields.
# One trade-off you have to make when loading data into NumPy arrays is that
# you cannot have arrays with a variable number of items per row.
# Here we only consider the first value for ALT.
callset = allel.read_vcf(myvcf, numbers={'ALT': 1})

# The callset object returned by read_vcf() is a Python dictionary (dict). 
# It contains several NumPy arrays, each of which can be accessed via a key. 
# Here are the available keys:
sorted(callset.keys())
# Results: 
# ['calldata/GT', 'samples', 'variants/ALT', 
# 'variants/CHROM', 'variants/FILTER_PASS', 
# 'variants/ID', 'variants/POS', 'variants/QUAL', 
# 'variants/REF']

# extract information for SNPs (chrom, position, A1, A2 )
chrom = callset['variants/CHROM']
pos = callset['variants/POS']
alt = callset['variants/ALT']
ref = callset['variants/REF']
SNP = ["_".join([str(chrom[i]),str(pos[i]),str(alt[i]),str(ref[i])]) for i in range(0,len(chrom))]

sample_id = callset['samples']
sample_id_clean = [x.split("_")[0] for x in sample_id]

# extract genotype arrays 
gt = allel.GenotypeArray(callset['calldata/GT'])
genotype = np.sum(gt, -1)

# caculate maf for each SNP
nsample = len(sample_id)


def calc_maf(x):
    maf = np.sum(x)/(2 * nsample)
    return maf
maf = np.apply_along_axis(calc_maf, 1, genotype)

# extract the SNP information and genotype data for the remaining SNPs
snpinfo_frame = pd.DataFrame({'Chrom': chrom, 'SNP': SNP, 'SNPPos': pos, 'A1': alt, 'A2': ref, 'MAF': maf},
                             columns=['Chrom', 'SNP', 'SNPPos', 'A1', 'A2', 'MAF'])
genotype_frame = pd.DataFrame(genotype, columns=sample_id_clean)

print(genotype_frame.shape)
# merge snp_info and genotype adta across columns
result = pd.concat([snpinfo_frame, genotype_frame], axis=1)

# save the results
result.to_csv(out_geno_path,
              sep='\t',
              index=None,
              header=True,
              mode='w')


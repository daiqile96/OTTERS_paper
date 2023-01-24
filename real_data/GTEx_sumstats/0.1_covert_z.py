
from scipy.stats import norm

import numpy as np
import pandas as pd

import getopt
import sys

from time import time
from scipy.stats import norm

import numpy as np
import pandas as pd

############################################################
# time calculation
start_time = time()


############################################################
def parse_param():
    long_opts_list = ['stat_dir=', 'out_dir=', 'N=', 'help']

    param_dict = {'stat_dir': None, 'out_dir': None, 'N': None}

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
            elif opt == "--stat_dir":
                param_dict['stat_dir'] = arg
            elif opt == "--out_dir":
                param_dict['out_dir'] = arg
            elif opt == "--N":
                param_dict['N'] = arg
    else:
        print(__doc__)
        sys.exit(0)

    if param_dict['stat_dir'] is None:
        print('* Please specify the directory to the summary statistics with p-values --stat_dir\n')
        sys.exit(2)
    elif param_dict['out_dir'] is None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict

############################################################

param_dict = parse_param()

# read data
stat_chunks = pd.read_csv(
    param_dict['stat_dir'],
    sep='\t',
    chunksize=50000,
    iterator=True,
    header=0,
    dtype={'pval_nominal': np.float64, 'slope': np.float64},
    compression='gzip')

stat = pd.concat([x for x in stat_chunks]).reset_index(drop=True)

print('Done Reading')

# get ID
stat['TargetID'] = [x.split('.')[0] for x in stat.gene_id]

# get position, chrom, A1, A2
variant_id = stat.variant_id
stat['CHROM'] = [x.split('chr')[1] for x in [x.split('_')[0] for x in variant_id]]
stat['POS'] = [x.split('_')[1] for x in variant_id]
stat['A1'] = [x.split('_')[3] for x in variant_id]
stat['A2'] = [x.split('_')[2] for x in variant_id]

# only select eQTL data for CHROM 1 - 22
stat = stat[stat.CHROM.isin([str(x) for x in range(1, 23)])]

# calculate Z
stat['Z'] = np.sign(stat.slope) * abs(norm.ppf(stat.pval_nominal/2.0))
stat['N'] = param_dict['N']

# change data types
stat['CHROM'] = stat.CHROM.astype('int32')
stat['POS'] = stat.POS.astype('int32')
stat['Z'] = stat.Z.astype('float64')
stat['N'] = stat.N.astype('int32')

# output
sort_stat = stat.sort_values(by=['CHROM', 'POS'], ascending=[True, True])
sort_stat[['CHROM', 'POS', 'A1', 'A2', 'Z', 'TargetID', 'N']].to_csv(
    param_dict['out_dir'],
    sep='\t',
    index=None,
    header=True,
    mode='w')




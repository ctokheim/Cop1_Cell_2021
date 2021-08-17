"""
File: merge_counts.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Merge counts for peaks
"""
import os
import glob
import argparse
import pandas as pd

def parse_arguments():
    info = 'Merge counts for peaks'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input-dir',
                        type=str, required=True,
                        help='Directory containing count data')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Merged output file')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # merge MUTs
    myfiles = glob.glob(os.path.join(opts['input_dir'], 'COP1_KO*_v2.txt'))
    suffix_len = len('_unique.sorted.bam_v2.txt')
    for ix, f in enumerate(myfiles):
        mybase = os.path.basename(f)[:-suffix_len]
        if ix==0:
            df = pd.read_csv(f, header=None, names=['chrom', 'start', 'end', mybase], sep='\t')
        else:
            tmp = pd.read_csv(f, header=None, names=['chrom', 'start', 'end', mybase], sep='\t')
            df = pd.merge(df, tmp, on=['chrom', 'start', 'end'], how='left')

    # add ROSA
    for f in glob.glob(os.path.join(opts['input_dir'], 'ROSA_*_v2.txt')):
        mybase = os.path.basename(f)[:-suffix_len]
        tmp = pd.read_csv(f, header=None, names=['chrom', 'start', 'end', mybase], sep='\t')
        df = pd.merge(df, tmp, on=['chrom', 'start', 'end'], how='left')

    # add WT data
    for f in glob.glob(os.path.join(opts['input_dir'], 'WT_*_v2.txt')):
        mybase = os.path.basename(f)[:-suffix_len]
        tmp = pd.read_csv(f, header=None, names=['chrom', 'start', 'end', mybase], sep='\t')
        df = pd.merge(df, tmp, on=['chrom', 'start', 'end'], how='left')

    # save output
    df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)



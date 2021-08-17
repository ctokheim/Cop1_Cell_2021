"""
File: split_diff_bed.py
Author: Collin Tokheim
Email: ctokheim@mail.dfci.harvard.edu
Github: ctokheim
Description: Split diff peaks into separated bed files
"""
import pandas as pd
import argparse
import os


def parse_arguments():
    info = 'Split diff peaks into separated bed files'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Diff peak result')
    parser.add_argument('-o', '--output-dir',
                        type=str, required=True,
                        help='Output directory')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    df = pd.read_csv(opts['input'], sep='\t')
    df = df[df['chrom'].str.startswith('chr')]

    # format col names
    cop1_col = 'COP1_KO'
    #rosa_col = 'ROSA'
    wt_col = 'WT'

    # check whether peak is up / down regulated
    is_up_cop1 = df[cop1_col]==1
    #is_up_rosa = df[rosa_col]==1
    is_up_wt = df[wt_col]==1
    is_down_cop1 = df[cop1_col]==-1
    #is_down_rosa = df[rosa_col]==-1
    is_down_wt = df[wt_col]==-1

    # save bed files
    mycols = ['chrom', 'start', 'end']
    mypath = os.path.join(opts['output_dir'], 'up_cop1_ko_vs_rosa.bed')
    df.loc[is_up_cop1, mycols].to_csv(mypath, header=None, index=None, sep='\t')
    mypath = os.path.join(opts['output_dir'], 'down_cop1_ko_vs_rosa.bed')
    df.loc[is_down_cop1, mycols].to_csv(mypath, header=None, index=None, sep='\t')
    mypath = os.path.join(opts['output_dir'], 'neutral_cop1_ko_vs_rosa.bed')
    df.loc[~(is_down_cop1) & ~(is_up_cop1), mycols].to_csv(mypath, header=None, sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)


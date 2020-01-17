# -*- coding:utf-8 -*-
from functools import reduce
import pandas as pd
import numpy as np
import optparse
import glob
import sys
import os


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    python merge_cnvnator_cnv_freq.py --all_sample all_sample_number --chrY_sample male_sample_number \
    --cnv_del_frq "chr*.deletion.cnv.freq.txt" \
    --cnv_dup_frq "chr*.duplication.cnv.freq.txt" \
    --chrY_del_frq "chrY.deletion.cnv.freq.txt" \
    --chrY_dup_frq "chrY.duplication.cnv.freq.txt" \
    --cnv_del "chr*/deletion/*.cnvnator.bed" \
    --cnv_dup "chr*/duplication/*.cnvnator.bed" \
    --chrY_del "chrY/deletion/*.cnvnator.bed" \
    --chrY_dup "chrY/duplication/*.cnvnator.bed"
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def cnv_freq(frq, sample_count):
    files = glob.glob(frq)
    stat = []
    for i in files:
        df = pd.read_csv(i, sep='\t')
        stat.append(df)
    stat_df = reduce(lambda x, y: x.append(y), stat)
    stat_df['Frequency'] = stat_df['Count'] / sample_count
    stat_df.sort_values(by=['#Chr', 'Start', 'End'], ascending=[True, True, True], inplace=True)
    return stat_df


def cnv_real_count_add_sample(cnvs):
    files = glob.glob(cnvs)
    stat = []
    for i in files:
        sample = os.path.basename(i).split('.')[0]
        df = pd.read_csv(i, sep='\t', header=None)
        df.columns = ['#Chr', 'Start', 'End', 'Size', 'Type']
        df['Samples'] = sample
        stat.append(df)
    stat_df = reduce(lambda x, y: x.append(y), stat)
    stat_df['real count'] = 0
    stat_df = stat_df.groupby(['#Chr', 'Start', 'End', 'Type']).agg({'real count': np.size, 'Samples': ';'.join}).reset_index()
    stat_df.sort_values(by=['#Chr', 'Start', 'End'], ascending=[True, True, True], inplace=True)
    return stat_df


def cnv_real_count(cnvs):
    files = glob.glob(cnvs)
    stat = []
    for i in files:
        df = pd.read_csv(i, sep='\t', header=None)
        df.columns = ['#Chr', 'Start', 'End', 'Size', 'Type']
        stat.append(df)
    stat_df = reduce(lambda x, y: x.append(y), stat)
    stat_df['real count'] = 0
    stat_df = stat_df.groupby(['#Chr', 'Start', 'End', 'Type'], as_index=False).agg({'real count': np.size})
    stat_df.sort_values(by=['#Chr', 'Start', 'End'], ascending=[True, True, True], inplace=True)
    return stat_df


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('--cnv_del', dest='cnv_del', default=None, type='string')
    parser.add_option('--cnv_dup', dest='cnv_dup', default=None, type='string')
    parser.add_option('--cnv_del_frq', dest='cnv_del_frq', default=None, type='string')
    parser.add_option('--cnv_dup_frq', dest='cnv_dup_frq', default=None, type='string')
    parser.add_option('--chrY_del', dest='chrY_del', default=None, type='string')
    parser.add_option('--chrY_dup', dest='chrY_dup', default=None, type='string')
    parser.add_option('--chrY_del_frq', dest='chrY_del_frq', default=None, type='string')
    parser.add_option('--chrY_dup_frq', dest='chrY_dup_frq', default=None, type='string')
    parser.add_option('--chrY_sample', dest='chrY_sample', type=int)
    parser.add_option('--all_sample', dest='all_sample', type=int)
    # parser.add_option('-p', '--pwd', dest='pwd', default=None, type='string')
    # parser.add_option('-o', '--out', dest='out_file', default='none', type='string')
    (opts, args) = parser.parse_args()
    # pwd = opts.pwd
    # out_file = opts.out_file
    cnv_del = opts.cnv_del
    cnv_dup = opts.cnv_dup
    cnv_del_frq = opts.cnv_del_frq
    cnv_dup_frq = opts.cnv_dup_frq
    chrY_del = opts.chrY_del
    chrY_dup = opts.chrY_dup
    chrY_del_frq = opts.chrY_del_frq
    chrY_dup_frq = opts.chrY_dup_frq
    chrY_sample = opts.chrY_sample
    all_sample = opts.all_sample
    df_del_cnv = cnv_freq(cnv_del_frq, all_sample)
    df_dup_cnv = cnv_freq(cnv_dup_frq, all_sample)
    df_del_cnv_chrY = cnv_freq(chrY_del_frq, chrY_sample)
    df_dup_cnv_chrY = cnv_freq(chrY_dup_frq, chrY_sample)
    all_cnv_del_frq = df_del_cnv.append(df_del_cnv_chrY)
    all_cnv_dup_frq = df_dup_cnv.append(df_dup_cnv_chrY)
    all_cnv_del_frq.to_csv("cnv.freq.del.with_chrY.txt", sep='\t', index=False)
    all_cnv_dup_frq.to_csv("cnv.freq.dup.with_chrY.txt", sep='\t', index=False)
    real_count_del = cnv_real_count_add_sample(cnv_del)
    real_count_dup = cnv_real_count_add_sample(cnv_dup)
    real_count_del_chrY = cnv_real_count_add_sample(chrY_del)
    real_count_dup_chrY = cnv_real_count_add_sample(chrY_dup)
    all_cnv_real_count_del = real_count_del.append(real_count_del_chrY)
    all_cnv_real_count_dup = real_count_dup.append(real_count_dup_chrY)
    all_cnv_real_count_del.to_csv("cnv.real_count.add_sample.del.with_chrY.txt", sep='\t', index=False)
    all_cnv_real_count_dup.to_csv("cnv.real_count.add_sample.dup.with_chrY.txt", sep='\t', index=False)
    cols = ['#Chr', 'Start', 'End', 'Type']
    merged_del = pd.merge(all_cnv_del_frq, all_cnv_real_count_del, on=cols, how='left')
    merged_dup = pd.merge(all_cnv_dup_frq, all_cnv_real_count_dup, on=cols, how='left')
    merged_del.to_csv("cnv.freq.real_count.add_sample.del.with_chrY.txt", sep='\t', index=False)
    merged_dup.to_csv("cnv.freq.real_count.add_sample.dup.with_chrY.txt", sep='\t', index=False)

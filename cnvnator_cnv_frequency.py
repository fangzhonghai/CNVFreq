# -*- coding:utf-8 -*-
import pandas as pd
import sys
import optparse
import numpy as np
import os


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    python cnvnator_cnv_frequency.py -p /path/to/work -o cnv.freq.txt --overlap 0.7 --cnv_type deletion \
    --cnv cnv.info.list
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def cnvnator_cnv2bed(cnv_file, sample, chrom, wkdir, del_dup):
    cnv_df = pd.read_csv(cnv_file, sep='\t', header=None)
    cnv_df['chr'] = cnv_df[1].str.split(":").str[0]
    cnv_df['interval'] = cnv_df[1].str.split(':').str[1]
    cnv_df['start'] = cnv_df['interval'].str.split('-').str[0]
    cnv_df['end'] = cnv_df['interval'].str.split('-').str[1]
    cnv_df = cnv_df[cnv_df[0] == del_dup].copy()
    cnv_df.to_csv(wkdir + "/" + sample + "." + chrom + ".cnvnator.bed", sep='\t', index=False, header=False, columns=['chr', 'start', 'end', 2, 0])


def cnv_intersect_tril(cnv_infos, wkdir):
    bedtools = "bedtools intersect -wo"
    for i in range(cnv_infos.shape[0]-1):
        for j in range(i+1, cnv_infos.shape[0]):
            sample_a = cnv_infos.iloc[i, 1]
            sample_b = cnv_infos.iloc[j, 1]
            chrom = cnv_infos.iloc[i, 2]
            cnv_a = wkdir + "/" + sample_a + "." + chrom + ".cnvnator.bed"
            cnv_b = wkdir + "/" + sample_b + "." + chrom + ".cnvnator.bed"
            command = bedtools + " -a " + cnv_a + " -b " + cnv_b + " > " + wkdir + "/" + sample_a + "." + sample_b + "." + chrom + ".intersect.bed"
            os.system(command)


def cnv_intersect(cnv_infos, wkdir):
    bedtools = "bedtools intersect -wo"
    for i in range(cnv_infos.shape[0]):
        for j in range(cnv_infos.shape[0]):
            if i != j:
                sample_a = cnv_infos.iloc[i, 1]
                sample_b = cnv_infos.iloc[j, 1]
                chrom = cnv_infos.iloc[i, 2]
                cnv_a = wkdir + "/" + sample_a + "." + chrom + ".cnvnator.bed"
                cnv_b = wkdir + "/" + sample_b + "." + chrom + ".cnvnator.bed"
                command = bedtools + " -a " + cnv_a + " -b " + cnv_b + " > " + wkdir + "/" + sample_a + "." + sample_b + "." + chrom + ".intersect.bed"
                os.system(command)


def cnv_count(cnv_infos, wkdir, per, out):
    sta_count = pd.DataFrame()
    sta_sample = {}
    cols = ['#Chr', 'Start', 'End', 'Size', 'Type', 'hit']
    for i in range(cnv_infos.shape[0]):
        tmp_df = pd.DataFrame()
        sample_a = cnv_infos.iloc[i, 1]
        chrom = cnv_infos.iloc[i, 2]
        for j in range(cnv_infos.shape[0]):
            if i != j:
                sample_b = cnv_infos.iloc[j, 1]
                cnv_inter = pd.read_csv(wkdir + "/" + sample_a + "." + sample_b + "." + chrom + ".intersect.bed", sep='\t', header=None)
                cnv_inter['hit'] = 0
                cnv_inter['over_rate_a'] = cnv_inter[10] / cnv_inter[3]
                cnv_inter['over_rate_b'] = cnv_inter[10] / cnv_inter[8]
                cnv_inter.loc[(cnv_inter['over_rate_a'] >= per) & (cnv_inter['over_rate_b'] >= per), 'hit'] = cnv_inter.loc[(cnv_inter['over_rate_a'] >= per) & (cnv_inter['over_rate_b'] >= per), 'hit'] + 1
                cnv_inter.sort_values(by='hit', ascending=False, inplace=True)
                cnv_inter_a = cnv_inter[[0, 1, 2, 3, 4, 'hit']].copy()
                cnv_inter_a.columns = cols
                cnv_inter_a.drop_duplicates(subset=cols[0:5], keep='first', inplace=True)
                tmp_df = tmp_df.append(cnv_inter_a)
        cnv_inter_a_all = pd.read_csv(wkdir + "/" + sample_a + "." + chrom + ".cnvnator.bed", sep='\t', header=None)
        cnv_inter_a_all['Count'] = 1
        cnv_inter_a_all.columns = cols[0:5] + ['Count']
        cnv_inter_a_all['Size'] = cnv_inter_a_all['Size'].astype(int)
        tmp_df['Count'] = tmp_df['hit']
        tmp_df['Size'] = tmp_df['Size'].astype(int)
        tmp_df_freq = tmp_df.groupby(cols[0:5], as_index=False).agg({'Count': np.sum})
        tmp_df_freq['Count'] = tmp_df_freq['Count'] + 1
        tmp_df_freq = tmp_df_freq.append(cnv_inter_a_all)
        tmp_df_freq.drop_duplicates(subset=cols[0:5], keep='first', inplace=True)
        sta_sample[sample_a] = tmp_df_freq
    for k in sta_sample.keys():
        sta_count = sta_count.append(sta_sample[k])
    sta_count.sort_values(by=['Start', 'End'], ascending=True, inplace=True)
    sta_count.drop_duplicates(inplace=True)
    sta_count.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('--cnv', dest='cnv', default=None, type='string')
    parser.add_option('-p', '--pwd', dest='pwd', default=None, type='string')
    parser.add_option('-o', '--out', dest='out_file', default=None, type='string')
    parser.add_option('--overlap', dest='overlap', type=float)
    parser.add_option('--cnv_type', dest='cnv_type', default='deletion', help='deletion or duplication', type='string')
    (opts, args) = parser.parse_args()
    cnv = opts.cnv
    pwd = opts.pwd
    out_file = opts.out_file
    overlap = opts.overlap
    cnv_type = opts.cnv_type
    cnv_info = pd.read_csv(cnv, header=None, sep='\t')
    for i in range(cnv_info.shape[0]):
        cnvnator_cnv2bed(cnv_info.iloc[i, 0], cnv_info.iloc[i, 1], cnv_info.iloc[i, 2], pwd, cnv_type)
    cnv_intersect(cnv_info, pwd)
    cnv_count(cnv_info, pwd, overlap, out_file)

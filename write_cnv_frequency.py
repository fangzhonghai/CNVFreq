# -*- coding:utf-8 -*-
import sys
import optparse
import os


def print_usage(option, opt, value, parser):
    usage_message = """
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    python write_cnv_frequency.py -p /path/to/work --overlap 0.7 --cnv /cnv/info/path
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def write_cnv_freq_sh(wkdir, per, cnv_info_path):
    chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX"]
    cnv_type = ["deletion", "duplication"]
    python = "/zfssz4/B2C_RD_P2/PMO/fangzhonghai/software/anaconda3/envs/python27/bin/python"
    py = "/zfssz4/B2C_RD_P2/PMO/fangzhonghai/py_scripts/cnvnator_cnv_frequency.py"
    for chrom in chroms:
        path = os.path.join(wkdir, chrom)
        if not os.path.exists(path):
            os.makedirs(path + "/deletion/scripts")
            os.makedirs(path + "/duplication/scripts")
        for ct in cnv_type:
            fp = open(path + "/" + ct + "/scripts/" + chrom + ".sh", "w")
            shell = '''#!/bin/bash
{python} {py} -p {path}/{ct} -o {wkdir}/{chrom}.{ct}.cnv.freq.txt --overlap {per} --cnv_type {ct} --cnv {cnv_info_path}/cnv.info.{chrom}.list
'''.format(**locals())
            fp.write(shell)
            fp.close()


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('--cnv', dest='cnv', default=None, type='string')
    parser.add_option('-p', '--pwd', dest='pwd', default=None, type='string')
    parser.add_option('--overlap', dest='overlap', type=float)
    (opts, args) = parser.parse_args()
    cnv = opts.cnv
    pwd = opts.pwd
    overlap = opts.overlap
    write_cnv_freq_sh(pwd, overlap, cnv)

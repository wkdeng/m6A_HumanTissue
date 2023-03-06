'''
read in peaks and compute t-statistic
ZZJ
8.27.2018
'''

import os
import pandas as pd
import numpy as np
import scipy.stats as ss
import pysam
from collections import defaultdict
import math
from tqdm import tqdm
import yaml
import sys
from multiprocessing import Pool
from functools import partial
import json 
import re


try:
    config = yaml.safe_load(open(sys.argv[1]))
except:
    sys.exit('command insufficient {config files}')

cur_dir = os.path.dirname(os.path.realpath(__file__))

SNAKEMAKE_FILE_DIR = config['SCRIPTS_DIR']
PROJECT_DIR = config['PROJECT_DIR']

OUT_DIR = PROJECT_DIR + '/summary/'
TMP_DIR = PROJECT_DIR + '/tmp/'
if not os.path.isdir(TMP_DIR):
    os.mkdir(TMP_DIR)
if not os.path.isdir(OUT_DIR):
    os.mkdir(OUT_DIR)

SAMPLE_TYPE_DICT = config['sample_type_dict']
COMPARISON_LIST = config['sample_comparison']
MAX_TAGS = config['clam']['max_tags']
IP_SUFFIX = config['ip_suffix']
CON_SUFFIX = config['con_suffix']
GTF_PATH=config[config['genome']]['gtf']


def read_peak_file(fn, peak_to_gene=None):
    peak_dict = {}
    if peak_to_gene is None:
        peak_to_gene = {}
    with open(fn, 'r') as f:
        for line in f:
            ele = line.strip().split()
            peak = ':'.join(ele[0:3])+':'+ele[5]
            score = float(ele[-4])
            gene = ele[3].split('-')[0]
            peak_dict[peak] = score
            peak_to_gene[peak] = gene
    return peak_dict, peak_to_gene


def read_project_clam(_project_dir):
    project_dir = os.path.join(_project_dir, 'clam/peaks')
    project_dict = {}
    peak_fn_list = [os.path.join(project_dir, x, 'narrow_peak.combined.bed.all') for x in os.listdir(project_dir) if
                    x.startswith('peaks-')]
    for peak_fn in peak_fn_list:
        peak_name = peak_fn.split('/')[-2].split(IP_SUFFIX + '__')[0].replace('peaks-', '').rstrip(
            '-')  # .rstrip('_IP')
        print( peak_name)
        project_dict[peak_name] = read_peak_file(peak_fn)[0]

    df = pd.DataFrame.from_dict(project_dict)

    sig_peaks = set()
    peak_to_gene = dict()
    peak_dict = dict()
    peak_fn_list2 = [os.path.join(project_dir, x, 'narrow_peak.combined.bed') for x in os.listdir(project_dir) if
                     x.startswith('peaks-')]
    for peak_fn in peak_fn_list2:
        print( peak_fn)
        peak_name = peak_fn.split('/')[-2].split(IP_SUFFIX + '__')[0].replace('peaks-', '').rstrip(
            '-')  # .rstrip('_IP')
        tmp, peak_to_gene = read_peak_file(peak_fn, peak_to_gene)
        _ = list(map(sig_peaks.add, tmp.keys()))
        peak_dict[peak_name] = tmp

    df_with_sig = df.loc[sig_peaks]
    peak_df = pd.DataFrame.from_dict(peak_dict)
    return df_with_sig, peak_to_gene, peak_df


def get_total_reads(bam_filename,unique=True):
    tot = 0

    if unique:
        idxstats = pysam.idxstats(bam_filename).split('\n')
        for l in idxstats:
            if not l:
                continue
            ele = l.split('\t')
            tot += int(ele[-2])
    else:
        for x in pysam.Samfile(bam_filename,'rb').fetch():
            tot += x.opt('AS')
    return tot


def get_reads_window(bam, chrom, start, end,unique=True):
    reads = [1 for x in bam.fetch(chrom, start, end)]
    if unique:
        return len(reads)
    else:
        return sum([x.opt('AS') for x in bam.fetch(chrom, start, end)])


def make_file_handlers(_project_dir):
    project_dir = os.path.join(_project_dir, 'clam')
    bam_dict = {}
    bam_fn_list = [[os.path.join(project_dir, x, 'unique.sorted.bam'), os.path.join(project_dir, x, 'realigned.sorted.bam')] 
                   for x in os.listdir(project_dir) if not x.startswith('peak') and os.path.isdir(project_dir + '/' + x)]

    for bam_fn in bam_fn_list:
        sample_name = bam_fn[0].split('/')[-2]
        # bam_dict[sample_name] = [pysam.Samfile(
        #     bam_fn[0], 'rb'), pysam.Samfile(bam_fn[1], 'rb')]
        bam_dict[sample_name] = bam_fn
    return bam_dict


def count_sample(peak_list, bam_dict, input_total_counts, ip_total_counts,sample):
    # print('%s started' % sample)
    # peak_list=peak_list[:1000]
    printed=False
    input_readcounts = pd.DataFrame(
        0, index=peak_list, columns=[sample])
    ip_readcounts = pd.DataFrame(
        0, index=peak_list, columns=[sample])
    input_name = sample + CON_SUFFIX
    ip_name = sample + IP_SUFFIX
    input_bams = [pysam.Samfile(x, 'rb') for x in bam_dict[input_name]]
    ip_bams = [pysam.Samfile(x, 'rb') for x in bam_dict[ip_name]]
    for peak in peak_list:
        if not printed:
            print(peak)
            printed=True
        chrom, start, end = peak.split(':')[:-1]
        start = int(start)
        end = int(end)
        this_input_u = get_reads_window(input_bams[0], chrom, start, end)
        this_input_m = get_reads_window(
            input_bams[1], chrom, start, end, False)
        this_input = this_input_u + this_input_m
        input_readcounts.loc[peak, sample] = \
            this_input / \
            float(input_total_counts[input_name]) / ((end - start) / 1000.)
        
        this_ip_u = get_reads_window(ip_bams[0], chrom, start, end)
        this_ip_m = get_reads_window(ip_bams[1], chrom, start, end, False)
        this_ip = this_ip_u+this_ip_m
        ip_readcounts.loc[peak, sample] = \
            this_ip / \
            float(ip_total_counts[ip_name]) / ((end - start) / 1000.)
    print('%s finished' % sample)
    input_readcounts.to_csv(PROJECT_DIR+'/parsed_peaks/tmp/%s_input.csv'%sample,sep='\t')
    ip_readcounts.to_csv(PROJECT_DIR+'/parsed_peaks/tmp/%s_ip.csv'%sample,sep='\t')
    # return [input_readcounts, ip_readcounts]


def read_gtf(fn):
    gene_annot = {}
    with open(fn, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            ele = line.strip().split('\t')
            if ele[2] != 'gene':
                continue
            chr, start, end, strand = ele[0], int(ele[3]), int(ele[4]), ele[6]
            try:
                gene_id = re.findall(r"(\w+)", ele[-1])[1]
            except AttributeError:
                continue
            gene_annot[gene_id] = [chr, start, end, strand, gene_id]
    return gene_annot


def get_gene_bin(chr, start, end, strand, bin_size):
    ret=[]
    current_end=start
    count=0
    while current_end<end:
        ret.append([chr,start+count*bin_size,start+(count+1)*bin_size,strand])
        current_end=start+(count+1)*bin_size
        count=count+1
    ret.append([chr,start+count*bin_size, start+(count+1)*bin_size,strand])
    return ret


def count_gene(peak_list, peak_to_gene, bam_dict, target_genes, gene_bins, sample):
    input_name = sample + CON_SUFFIX
    ip_name = sample + IP_SUFFIX
    input_bams = [pysam.Samfile(x, 'rb') for x in bam_dict[input_name]]
    ip_bams = [pysam.Samfile(x, 'rb') for x in bam_dict[ip_name]]

    gene_bin_count_ip = defaultdict(lambda: defaultdict(float))
    gene_bin_count_inp = defaultdict(lambda: defaultdict(float))
    p_=0
    for gene in target_genes:
        print('processing %s, %s/%s'%(gene,p_,len(target_genes)))
        bins = gene_bins[gene]
        for bin_ in bins:
            bin_id = ':'.join([str(x) for x in bin_])
            gene_bin_count_ip[gene][bin_id] = get_reads_window(
                ip_bams[0], bin_[0], bin_[1], bin_[2], unique=True)
            gene_bin_count_ip[gene][bin_id] += get_reads_window(
                ip_bams[1], bin_[0], bin_[1], bin_[2], unique=False)
            gene_bin_count_inp[gene][bin_id] = get_reads_window(
                input_bams[0], bin_[0], bin_[1], bin_[2], unique=True)
            gene_bin_count_inp[gene][bin_id] += get_reads_window(
                input_bams[1], bin_[0], bin_[1], bin_[2], unique=False)
        p_+=1
    json.dump(gene_bin_count_ip, open(os.path.join(
        PROJECT_DIR, 'parsed_peaks/tmp2', 'gene_bin_count_ip_%s.json' % sample), 'w'))
    json.dump(gene_bin_count_inp, open(os.path.join(
        PROJECT_DIR, 'parsed_peaks/tmp2', 'gene_bin_count_inp_%s.json' % sample), 'w'))
    
    gene_median_ip = defaultdict(float)
    gene_median_inp = defaultdict(float)
    for gene in target_genes:
        gene_ip_rc = list(gene_bin_count_ip[gene].values())
        gene_ip_rc.sort()
        gene_median_ip[gene] = gene_ip_rc[int(len(gene_ip_rc)/2)]+1

        gene_inp_rc = list(gene_bin_count_inp[gene].values())
        gene_inp_rc.sort()
        gene_median_inp[gene] = gene_inp_rc[int(len(gene_inp_rc)/2)]+1

    winscore_ip = pd.DataFrame(
        0, index=peak_list, columns=[sample])
    winscore_inp = pd.DataFrame(
        0, index=peak_list, columns=[sample])

    for peak in peak_list:
        gene = peak_to_gene[peak]
        if gene not in gene_bin_count_ip.keys():
            continue
        ip_ = gene_bin_count_ip[gene][peak]+1
        winscore_ip.loc[peak, sample] = ip_/gene_median_ip[gene]

        inp_ = gene_bin_count_inp[gene][peak]+1
        winscore_inp.loc[peak, sample] = inp_/gene_median_inp[gene]

    winscore_inp.to_csv(
        PROJECT_DIR+'/parsed_peaks/tmp2/%s_input_winscore.csv' % sample, sep='\t')
    winscore_ip.to_csv(
        PROJECT_DIR+'/parsed_peaks/tmp2/%s_ip_winscore.csv' % sample, sep='\t')


def count_sep(df_with_sig, bam_dict, project_dir,peak_to_gene,sample,bin_size=100):

    target_genes = list(set(peak_to_gene.values()))
    gene_anno = read_gtf(GTF_PATH)

    gene_bins={}
    for gene in target_genes:
        gene_bins[gene] = get_gene_bin(
            gene_anno[gene][0], gene_anno[gene][1], gene_anno[gene][2], gene_anno[gene][3], bin_size)
    
    
    if not os.path.isdir(PROJECT_DIR+'/parsed_peaks/tmp2'):
        os.mkdir(os.path.join(PROJECT_DIR,'parsed_peaks/tmp2'))
    peak_list = df_with_sig.index.tolist()
    
    count_gene(peak_list, peak_to_gene,
                                  bam_dict, target_genes, gene_bins,sample)

def count_mp(df_with_sig, bam_dict, project_dir,peak_to_gene,bin_size=100):
    '''
    Given a list of peaks and samples, count the RPKM,
    then save to file.
    '''
    input_readcounts = pd.DataFrame(
        0, index=df_with_sig.index, columns=df_with_sig.columns)
    ip_readcounts = pd.DataFrame(
        0, index=df_with_sig.index, columns=df_with_sig.columns)
    
    
    
    # target_genes = list(set(peak_to_gene.values()))
    # gene_anno = read_gtf(GTF_PATH)

    # gene_bins={}
    # for gene in target_genes:
    #     gene_bins[gene] = get_gene_bin(
    #         gene_anno[gene][0], gene_anno[gene][1], gene_anno[gene][2], gene_anno[gene][3], bin_size)
    
    
    # if not os.path.isdir(PROJECT_DIR+'/parsed_peaks/tmp'):
    #     os.mkdir(os.path.join(PROJECT_DIR,'parsed_peaks/tmp'))
    # peak_list = df_with_sig.index.tolist()
    
    # counted=0
    # for sample in df_with_sig.columns.tolist():
    #     count_gene(peak_list, peak_to_gene,
    #                               bam_dict, target_genes, gene_bins,sample)
    #     counted+=1
    #     print('finished sample %s, %s/%s' %
    #           (sample, counted, len(df_with_sig.columns.tolist())))
    # pool=Pool(processes=20)
    # func = partial(count_gene, peak_list, peak_to_gene,
    #                bam_dict, target_genes, gene_bins)
    # results=pool.map(func,df_with_sig.columns.tolist())
    
    # pool.close()
    # pool.join()
    input_winscore=pd.concat(
        [pd.read_csv(os.path.join(PROJECT_DIR, 'parsed_peaks/tmp2', '%s_input_winscore.csv'%x),sep='\t',index_col=0) for x in df_with_sig.columns],axis=1)
    ip_winscore = pd.concat(
        [pd.read_csv(os.path.join(PROJECT_DIR, 'parsed_peaks/tmp2', '%s_ip_winscore.csv'%x),sep='\t',index_col=0) for x in df_with_sig.columns],axis=1) 

    return input_winscore, ip_winscore


def get_peak_intensity(ip_readcounts, input_readcounts):
    ratio = pd.DataFrame(0, index=ip_readcounts.index, columns=ip_readcounts.columns)
    for peak in ratio.index:
        input_rpkm = input_readcounts.loc[peak]
        ip_rpkm = ip_readcounts.loc[peak]
        ratio.loc[peak] = (ip_rpkm) / (input_rpkm)
    return ratio


def run_and_save():
    # read in clam peaks
    print('read peak')
    df_with_sig, peak_to_gene, peak_df = read_project_clam(PROJECT_DIR)

    # get peak intensities
    bam_dict = make_file_handlers(PROJECT_DIR)
    input_readcounts, ip_readcounts = count_mp(
        df_with_sig, bam_dict, PROJECT_DIR, peak_to_gene)

    # save the RPKM values
    if not os.path.isdir(PROJECT_DIR + '/parsed_peaks'):
        os.mkdir(PROJECT_DIR + '/parsed_peaks')
    input_readcounts.to_csv(PROJECT_DIR + '/parsed_peaks/input_peak.winscore.csv', sep='\t')
    ip_readcounts.to_csv(
        PROJECT_DIR + '/parsed_peaks/ip_peak.winscore.csv', sep='\t')

    # load previously counts
    input_readcounts = pd.read_table(PROJECT_DIR + '/parsed_peaks/input_peak.winscore.csv', index_col=0)
    ip_readcounts = pd.read_table(PROJECT_DIR + '/parsed_peaks/ip_peak.winscore.csv', index_col=0)

    # compute the peak intensity
    ratio = get_peak_intensity(ip_readcounts, input_readcounts)
    ratio.to_csv(
        PROJECT_DIR + '/parsed_peaks/peak_intensity_winscore.csv', sep='\t')

def run_sep():
    sample=sys.argv[2]
    df_with_sig, peak_to_gene, peak_df = read_project_clam(PROJECT_DIR)

    bam_dict = make_file_handlers(PROJECT_DIR)
    count_sep(
        df_with_sig, bam_dict, PROJECT_DIR, peak_to_gene,sample)

if __name__ == '__main__':
    run_and_save()
    # run_sep()

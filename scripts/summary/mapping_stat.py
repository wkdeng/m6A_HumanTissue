# read star output for mapping stats
# Zijun Zhang
# 8.11.2017
# revised 9.22.2017: incorporate into Snakemake
# revised 3.3.2018: update exon read criteria

import sys
import os
import subprocess
import shutil


def read_star_log(fn):
    assert os.path.isfile(fn)
    stats = {}
    stats_list = [
        'Number of input reads',
        #'Uniquely mapped reads number',
        'Uniquely mapped reads %',
        #'Number of reads mapped to multiple loci',
        '% of reads mapped to multiple loci',
        'Number of splices: Total',
        '% of reads unmapped: too short'
    ]
    with open(fn, 'r') as f:
        for line in f:
            ele = [x.strip() for x in line.strip().split('|')]
            if ele[0] in stats_list:
                stats[ele[0]] = ele[1]
    return stats


def read_junction_reads(fn):
    assert os.path.isfile(fn)
    cmd = '''samtools view -q 10 %s | awk -F"\t" '$6~"N"&&$6!~"D"&&$6!~"I"&&$0~"NH:i:1"' | wc -l''' % fn
    sj = subprocess.check_output(cmd, shell=True).decode('utf-8')
    return sj.strip()


def read_exon_reads(fn, read_len=None):
    assert os.path.isfile(fn)
    # assert read_len is None or type(read_len) == int
    # if read_len is None:
    #     cmd2 = "samtools view -q 10 %s | head -n 1 | awk '{print length($10)}'" % fn
    #     print >>sys.stderr, cmd2
        
    #     read_len = subprocess.check_output(cmd2, shell=True)
    #     read_len = read_len.strip()
    #cmd = '''samtools view -q 10 %s | awk '$6=="%sM"&&$0~"NH:i:1"' | wc -l'''%(fn,read_len)
    cmd = '''samtools view -q 10 %s | awk '$0~/NH:i:1\>/&&$6~"M"&&$6!~"N"' | wc -l''' % (
        fn)
    #print >>sys.stderr, cmd
    sj = subprocess.check_output(cmd, shell=True).decode('utf-8')
    return sj.strip()


def read_barcode(fn):
    bar_dict = {}
    with open(fn, 'r') as f:
        firstline = True
        for line in f:
            ele = line.strip().split()
            if firstline:
                header = {ele[x]: x for x in range(len(ele))}
                firstline = False
                continue
            bar_dict[ele[header['barcode']]] = ele[1:]
    return bar_dict


def read_fastqc(file_):
    if os.path.isfile(file_[:-4]):
        os.remove(file_[:-4])
    elif os.path.isdir(file_[:-4]):
        shutil.rmtree(file_[:-4])
    subprocess.call(
        'unzip -q -d {} {}'.format('/'.join(file_[:-4].split('/')[:-1]), file_), shell=True)
    qc_data = open(os.path.join(file_[:-4], 'fastqc_data.txt'))
    per_baseqc = 0.0
    record_num = 0
    module_flag = False
    gc_content = 0
    for line in qc_data:
        if line.startswith(r'%GC'):
            gc_content = line.strip().split('\t')[1]
        if module_flag and line.startswith('>>END_MODULE'):
            break
        if line.startswith('>>Per base sequence quality'):
            module_flag = True
            continue
        if module_flag and not line.startswith('#'):
            per_baseqc = per_baseqc+float(line.strip().split('\t')[1])
            record_num = record_num+1

    return per_baseqc/record_num, gc_content


def mapping_stat(bam_fn, qc_zip, verbose=True):
    stats_list = [
        'Number of input reads',
        #'Uniquely mapped reads number',
        'Uniquely mapped reads %',
        #'Number of reads mapped to multiple loci',
        '% of reads mapped to multiple loci',
        'Number of splices: Total',
        'Number of splice junction reads',
        'Number of exon reads',
        '% of reads unmapped: too short', 
        'Average quality score',
        '% of GC content'
    ]
    #if verbose: print '\t'.join(['Dataset']+stats_list)
    qc_data = read_fastqc(qc_zip)
    source_dir = os.path.dirname(bam_fn)
    stat_fn = os.path.join(source_dir, 'Log.final.out')
    stats = read_star_log(stat_fn)
    sj = read_junction_reads(bam_fn)
    #ec = str(int(stats['Uniquely mapped reads number'])*2 - int(sj))
    ec = read_exon_reads(bam_fn)
    stats['Number of splice junction reads'] = sj
    stats['Number of exon reads'] = ec
    stats['Average quality score'] = qc_data[0]
    stats['% of GC content'] = qc_data[1]
    if verbose:
        for x in stats_list:
            print( x + '\t' + str(stats[x]))
    return stats


if __name__ == '__main__':
    if len(sys.argv) > 1:
        bam_fn = sys.argv[1]
        qc_zip = sys.argv[2]
        mapping_stat(bam_fn, qc_zip)
    else:
        sys.exit(1)

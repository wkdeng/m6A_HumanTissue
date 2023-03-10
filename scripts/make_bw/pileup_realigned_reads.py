# @author [Wankun Deng]
# @email [dengwankun@hotmail.com]
# @create date 2020-03-09 17:25:41
# @modify date 2020-03-09 17:25:41
# @desc [Generate BigWig track of read signal using realigned multiple mapped reads and unique reads generated by CLAM realigner]

import pysam
import os
import re
import sys
import numpy as np
import yaml
import subprocess
import multiprocessing
from collections import defaultdict
from optparse import OptionParser, OptionGroup
import logging
options = None
logger = None

def read_gtf(fn):
    """read in the gene annotation from GTF file
    """
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


def parse_chunk_gene(args):
    bam_list = args[0]
    gene_anno = args[1]
    genes = args[2]
    unique_only = args[3]
    scale_factor=args[4]
    unstranded=args[5] if len(args)>=6 else False
    rt_method = args[6] if len(args) > 6 else 'start'
    crosslinked_only=args[7] if len(args)>=8 else False
    bam_list = [pysam.Samfile(x, 'rb') for x in bam_list]
    positive_strand = []
    negative_strand = []
    for gene in genes:
        chr, start, end = gene_anno[gene][:3]
        gene_intervals = [np.zeros(end-start+1),np.zeros(end-start+1)]
        for bam, unique in zip(bam_list, [False, True]):
            if unique_only and not unique:
                continue
            for x in bam.fetch(chr, start, end):
                ## if unstranded, signals all added to positive_strand
                idx=1 if x.is_reverse and not unstranded else 0 
                read_start = x.opt('RT')
                read_end = x.opt('RL')+read_start-1
                if rt_method == 'median':
                    read_start = x.opt('RT')-int(int(x.opt('RL')/2))
                    read_end = x.opt('RL')+read_start-1
                if crosslinked_only:
                    region_start=read_start-start
                    region_end=region_start+1
                else:
                    region_start = read_start-start
                    region_start = region_start if region_start >= 1 else 1
                    region_end = read_end-start
                    region_end = region_end if region_end < len(
                        gene_intervals[idx]) else len(gene_intervals[idx])-1
                to_add = 1 if unique else x.opt('AS')
                gene_intervals[idx][region_start-1:region_end] = np.add(
                    gene_intervals[idx][region_start-1:region_end], to_add*scale_factor)
        positive_strand.append([start, end, chr, gene_intervals[0]])
        negative_strand.append([start, end, chr, gene_intervals[1]])
        # if len(positive_strand)>100 and len(negative_strand)>100:
        #     break
    _ = map(lambda x: x.close(), [x for x in bam_list])

    return [positive_strand, negative_strand]


def chunkify(a, n):
    k, m = len(a) / n, len(a) % n
    return (a[int(i * k + min(i, m)):int((i + 1) * k + min(i + 1, m))] for i in range(n))


def process_one_chr(signal_list):
    covered_position = []
    ret = ''
    count = 0
    for item in signal_list:
        count = count+1
        current_position = item[0]
        current_chr = item[2]
        if len(ret) == 0:
            ret = ret+('variableStep   chrom=%s\n' % (current_chr))

        for read_num in item[3]:
            if current_position not in covered_position:
                if read_num > 0:
                    ret = ret+('%s %s\n' % (current_position, read_num))
                    covered_position.append(current_position)
            else:
                print( 'site %s:%s covered' % (current_chr, current_position))
            current_position = current_position+1
        if not count % 100:
            print( '%s progress: %s/%s' % (current_chr, count, len(signal_list)))
    return ret


def get_non_overlap_region(gene_anno):
    chr_dict = defaultdict(list)
    result_list = []
    for item in gene_anno.values():
        chr_dict[item[0]].append(item)
    for chr in chr_dict:
        temp_list = chr_dict[chr]
        temp_list.sort(key=lambda x: x[1])
        last_end = 0
        tmp_sorted = []
        for item in temp_list:
            chr, start, end, strand, gene_id = item
            if start <= last_end:
                if end <= last_end:
                    continue
                else:
                    start = last_end+1
            last_end = max(end, last_end)
            if not start > end:
                tmp_sorted.append([chr, start, end, strand, gene_id])
        result_list.extend(tmp_sorted)
    result_dict = {}
    for item in result_list:
        result_dict[item[-1]] = item
    return result_dict


def parse_one_sample(bam_list, chr_size, annotation, output_folder, sample, threadN, unique_only,unstranded,rt_method,scale,crosslinked_only):
    if output_folder.endswith('/'):
        output_folder=output_folder[:-1]
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    gene_anno = read_gtf(annotation)
    gene_anno = get_non_overlap_region(gene_anno)

    ## scale read signals by total read number.
    scale_factor=1
    if scale:
        scale_factor=float(subprocess.check_output('samtools view -F 260 %s| wc -l'%bam_list[1],shell=True).decode('utf-8'))
        if not unique_only:
            scale_factor += float(subprocess.check_output('samtools view -F 260 %s| wc -l' % bam_list[0], shell=True).decode('utf-8'))
        scale_factor=1E06/scale_factor
    
    for bam in bam_list:
        if not os.path.isfile(bam+'.bai'):
            subprocess.check_output('samtools index %s' % bam, shell=True)
    child_gene_list = [x for x in chunkify(list(gene_anno.keys()), threadN)]
    pool = multiprocessing.Pool(processes=threadN)
    list_list = pool.map(parse_chunk_gene, [
        (bam_list, gene_anno, child_gene_list[i], unique_only, scale_factor, unstranded, rt_method,crosslinked_only) for i in range(threadN)])
    pool.terminate()
    pool.join()

    positive_strand = []
    negative_strand = []
    for item in list_list:
        positive_strand.extend(item[0])
        negative_strand.extend(item[1])
    print( 'Strand extracted')
    
    pile_up_p = '%s/%s_pile_up_positive_%s.txt' % (
        output_folder,sample, 'unique' if unique_only else 'multi')

    wiggle_p = open(pile_up_p, 'w')
    for item in positive_strand:
        wiggle_p.write('fixedStep chrom=%s start=%s step=1\n' %
                       (item[2], item[0]))
        for _x in item[3]:
            wiggle_p.write('%s\n' % _x)
    wiggle_p.close()
    
    subprocess.call('wigToBigWig %s %s %s/%s_%s_plus.bw' %
                    (pile_up_p, chr_size, output_folder, sample, 'unique' if unique_only else 'multi'), shell=True)
    print( 'Positive track generated')
    os.remove(pile_up_p)
    if unstranded:
        subprocess.call('mv {0}/{1}_{2}_plus.bw {0}/{1}_{2}.bw'.format(output_folder,
                                                 sample, 'unique' if unique_only else 'multi'),shell=True)


    if not unstranded:
        pile_up_n = '%s/%s_pile_up_negative_%s.txt' % (
            output_folder, sample, 'unique' if unique_only else 'multi')
        wiggle_n = open(pile_up_n, 'w')
        for item in negative_strand:
            wiggle_n.write('fixedStep chrom=%s start=%s step=1\n' %
                        (item[2], item[0]))
            for _x in item[3]:
                wiggle_n.write('%s\n' % _x)
        wiggle_n.close()
        subprocess.call('wigToBigWig %s %s %s/%s_%s_minus.bw' %
                        (pile_up_n, chr_size, output_folder, sample, 'unique' if unique_only else 'multi'), shell=True)
        print( 'Negative track generated')
        os.remove(pile_up_n)

# if __name__ == '__main__':
    # sample = sys.argv[1]
    # unique_only = len(sys.argv) >= 3 and sys.argv[2] == 'True'
    # parse_one_sample(['/home/dengw1/workspace/eCLIP/clam/%s/realigned.sorted.bam' % sample, '/home/dengw1/workspace/eCLIP/clam/%s/unique.sorted.bam' % sample],
    #                  '/home/dengw1/workspace/genome/raw/hg19/chr_size.txt', '/home/dengw1/workspace/genome/raw/hg19/hg19.gtf',
    #                  '/home/dengw1/workspace/eCLIP/bigwig/clam_realigned/%s' % sample, threadN=20, unique_only=unique_only)


def parse():
    global logger, options

    parser = OptionParser()
    basicGroup = OptionGroup(parser, "Basic Options")

    basicGroup.add_option('-s', '--sample', dest='sample_name', #required=True,
                          help='Name of this sample')

    basicGroup.add_option('--ub',  dest='unique_bam', #required=True,
                          help='Uniquely mapped reads BAM')

    basicGroup.add_option('--rb',  dest='realigned_bam',#required=True,
                          help='Realigned reads BAM')

    basicGroup.add_option('--uo',  dest='unique_only', default=False, 
                          help='Using unique reads to build BigWig')

    basicGroup.add_option('-g', '--gtf', dest='gtf', #required=True,
                          help='Genome annotation GTF file')

    basicGroup.add_option('-c', '--chr', dest='chr_size',# required=True,
                          help='Chromosome size file')
    
    basicGroup.add_option('-u', '--unstranded', dest='unstranded', default=False, action="store_true",
                          help='merge forward and reverse strand and output only one bigwig file')

    basicGroup.add_option(
        '-t', '--thread', dest='threadN', default=1, type='int', help='Thread number')

    basicGroup.add_option('-o', '--output', dest='output_folder',# required=True,
                          help='Output folder')

    basicGroup.add_option('-m', '--read-tag', dest='rt_method', default='start',# required=False,
                          help='Read tag method (start or median)')
    basicGroup.add_option(
        '-q', '--quiet', dest='quiet', action="store_true",  default=False, help='Hide log information')

    basicGroup.add_option('--scale',  dest='scale', default=False, action="store_true",
                          help='Scale to RPM')
    basicGroup.add_option('--crosslinked_only',  dest='crosslinked_only', default=False, action="store_true",
                          help='Only add read to crosslinking site')

    parser.add_option_group(basicGroup)
    options, _ = parser.parse_args()

    if options.quiet:
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "ERROR"))
    else:
        logging.basicConfig(level=os.environ.get("LOGLEVEL", "DEBUG"))
    logger = logging.getLogger(__name__)

    parse_one_sample([options.realigned_bam, options.unique_bam], options.chr_size, options.gtf, options.output_folder, 
                    options.sample_name,  options.threadN, options.unique_only, options.unstranded, options.rt_method, 
                    options.scale,options.crosslinked_only)


if __name__ == '__main__':
    parse()

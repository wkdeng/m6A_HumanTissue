import yaml
from pybedtools import BedTool
from collections import defaultdict
import os

config = yaml.safe_load(
    open('/home/dengw1/workspace/snakerun/eCLIP/config.yaml'))
yrna = BedTool('/home/dengw1/workspace/genome/raw/hg19/YRNA.bed')
yrna_group = defaultdict(lambda: defaultdict(list))
out_file = open('yrna_group.txt', 'w')
# expression=defaultdict(list)
t2g = defaultdict(list)
for line in open('/home/dengw1/workspace/genome/raw/hg19/tx2gene.bed'):
    info = line.strip().split('\t')
    t2g[info[1]] = info[0]
ensg_gn = defaultdict(lambda: 'N/A')
expression = defaultdict(lambda: [0, 0])
for line in open('/home/dengw1/workspace/genome/raw/hg19/ensg_gn.txt'):
    info = line.strip().split('\t')
    ensg_gn[info[1]] = info[0]
for group in config['groups']:
    exp_file = '/home/dengw1/workspace/eCLIP/kallisto/{}_Inp_rep1/abundance.tsv'.format(
        group)
    gene_ = group.split('_')
    index = 0 if gene_[0] == 'K562' else 1
    if os.path.isfile(exp_file):
        for line in open(exp_file):
            info = line.strip().split('\t')
            tx = info[0].split('|')[0].split('.')[0]
            ensg = t2g[tx]
            if ensg == ensg_gn[gene_[1]]:
                expression[gene_[1]][index] = expression[gene_[
                    1]][index]+float(info[-1])
    print( group)
    peak = BedTool(
        '/home/dengw1/workspace/eCLIP/clam/peaks/{}/narrow_peak.combined.bed'.format(group))
    inter = peak.intersect(yrna, f=0.5, wa=True, wb=True)
    for item in inter:
        if group not in yrna_group[item[13]]:
            yrna_group[item[13]][group] = [item[6], item[7], item[8]]
# order: K562, HepG2
collapsed_list = defaultdict(lambda: defaultdict(
    lambda: [['N/A', 'N/A'], ['N/A', 'N/A'], ['N/A', 'N/A']]))
for yrna in yrna_group:
    for group in yrna_group[yrna]:
        gene_ = group.split('_')
        index = 0 if gene_[0] == 'K562' else 1
        collapsed_list[yrna][gene_[1]][0][index] = yrna_group[yrna][group][0]
        collapsed_list[yrna][gene_[1]][1][index] = yrna_group[yrna][group][1]
        collapsed_list[yrna][gene_[1]][2][index] = yrna_group[yrna][group][2]
for yrna in collapsed_list:
    for gene in collapsed_list[yrna]:
        out_array = [collapsed_list[yrna][gene][0][0],
                     collapsed_list[yrna][gene][1][0], collapsed_list[yrna][gene][2][0], str(
                         expression[gene][0]), collapsed_list[yrna][gene][0][1],
                     collapsed_list[yrna][gene][1][1], collapsed_list[yrna][gene][2][1], str(expression[gene][1])]
        if out_array[0] == 'N/A':
            out_array[:4] = ['', '', '', '']
        if out_array[4] == 'N/A':
            out_array[4:] = ['', '', '', '']
        out_file.write(yrna+'\t'+gene+'\t' + '\t'.join(out_array) + '\n')
out_file.close()

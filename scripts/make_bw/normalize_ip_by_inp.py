###
## @author [Wankun Deng]
## @email [dengwankun@gmail.com]
## @create date 2023-03-06 01:12:31
## @modify date 2023-03-06 01:12:31
## @desc [description]
###
import pyBigWig
import sys
import re
import pybedtools
import pandas as pd


def read_gtf(fn):
    gene_annot = []
    with open(fn, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            ele = line.strip().split('\t')
            if ele[2] != 'gene':
                continue
            chr, start, end = ele[0], int(ele[3]), int(ele[4])
            gene_annot.append([chr, str(start), str(end)])
    gene_annot = pd.DataFrame(gene_annot, columns=['chr', 'start', 'end'])
    bed_ = pybedtools.BedTool.from_dataframe(gene_annot)
    bed_ = bed_.sort()
    bed_ = bed_.merge()
    return bed_


ip_bw = pyBigWig.open(sys.argv[1])
inp_bw = pyBigWig.open(sys.argv[2])
gene_anno = read_gtf(sys.argv[3])


bw_header = inp_bw.header()
out_bw = pyBigWig.open(sys.argv[4], "wb")

header = []
for chr in inp_bw.chroms():
    header.append((chr, inp_bw.chroms()[chr]))
out_bw.addHeader(header)

count = 0
for gene in gene_anno:
    count += 1
    chr, start, end = gene.fields[0], int(gene.fields[1]), int(gene.fields[2])
    ip_val = ip_bw.values(chr, start, end)
    inp_val = inp_bw.values(chr, start, end)
    values = [(x+1)/(y+1) if x and y else None for x,
              y in zip(ip_val, inp_val)]
    try:
        out_bw.addEntries(chr, start, span=1, step=1, values=values)
    except RuntimeError:
        print("Error: ", chr, start, end)
        continue

out_bw.close()

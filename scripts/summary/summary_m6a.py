# @author [Zijun Zhang, modified by Wankun Deng]
# @email [dengwankun@hotmail.com]
# @create date 2019-12-10 15:05:06
# @modify date 2019-12-10 15:05:06
# @desc [description]
import sys
import re
import numpy as np
import pandas as pd
import os
from subprocess import *
import time
from collections import defaultdict
from scipy.stats import zscore
import yaml
from PyPDF2 import PdfFileMerger, PdfFileReader


try:
    config = yaml.safe_load(open(sys.argv[1]))
except:
    sys.exit('command insufficient {config files}')

cur_dir = os.path.dirname(os.path.realpath(__file__))
PROJECT_DIR = config['PROJECT_DIR']
SNAKE_SCRIPT_DIR=config['SCRIPTS_DIR']
OUT_DIR = PROJECT_DIR + '/summary/'
TMP_DIR = PROJECT_DIR + '/tmp/'
GENOME = config['genome']
SAMPLE_TYPE_DICT = config['sample_type_dict']
COMPARISON_LIST = config['sample_comparison']
CON_SUFFIX = config['con_suffix']
T2G = config['t2g']
REGULATORS = config['regulator']
ENSG2GN = config['ensg2gn']
DEG_GROUP = config['deg_group'] if 'deg_group' in config else ''
NAME_SEP = config['name_sep'] if 'name_sep' in config else '-'
IS_TRIPSEQ = config['is_tripseq']==True if 'is_tripseq' in config else False
sys.path.insert(0, config['SCRIPTS_DIR']+'/scripts')

from peakComposition.peakPieCharts import intersect_gtf_regions
from GO_enrichment import hyper_geometry

def read_file(filename):
    with open(filename) as fp:
        List = [x.strip() for x in fp if len(x.strip()) > 0]
        return List


def convert2log2pwm(pwm, scale=False):
    if scale:
        scale_size = 0.01
    else:
        scale_size = 0.00
    array = np.array([[float(_x) for _x in x.split('\t')]
                      for x in pwm.strip().split('\n')[1:]])
    array += scale_size
    array /= np.sum(array, axis=1).reshape((len(array), 1))
    array /= 0.25
    return array


def star_mapping_stats():
    fw = open(OUT_DIR + 'mapping_stats.counts.summary.txt', 'w')
    fw.write(
        'sample\ttype\ttotal reads\tsplice junction reads\texon reads\tuniquely mapped reads (%)\tmultiple mapped reads(%)\tsequencing quality (avg)\tGC content(%)\n')
    star_dir = PROJECT_DIR + '/star/'
    for sample in sorted(SAMPLE_TYPE_DICT.keys()):
        sample_type = SAMPLE_TYPE_DICT[sample]
        if sample_type == 'ip':
            sample_type = 'IP'
        if sample_type == 'con':
            sample_type = 'Inp'
        mapping_info = read_file(star_dir + sample + '/mapping_stats.txt')
        input_reads = "{:,}".format(int(mapping_info[0].split()[-1]))
        splice_junction_reads = "{:,}".format(int(mapping_info[4].split()[-1]))
        exon_reads = "{:,}".format(int(mapping_info[5].split()[-1]))
        unique_mapped = mapping_info[1].split()[-1]
        multip_mapped = mapping_info[2].split()[-1]
        qc_quality = '%.2f' % float(mapping_info[7].split()[-1])
        gc_content = mapping_info[8].split()[-1]
        fw.write('\t'.join(
            [sample, sample_type, input_reads, splice_junction_reads, exon_reads, unique_mapped, multip_mapped, qc_quality, gc_content]) + '\n')
    fw.close()


def table_generate(table_inp, out_file):
    table_row = len(read_file(table_inp))
    print( "library(gtable)")
    print( "library(gridExtra)")
    print( "library(grid)")
    print( "data <- read.table('{}', header = T, sep = '\\t', check.names=FALSE)".format(table_inp))
    print("g <- tableGrob(data, rows = NULL)")
    print("g <- gtable_add_grob(g,")
    print("                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),")
    print("                     t = 2, b = nrow(g), l = 1, r = ncol(g))")
    print("g <- gtable_add_grob(g,")
    print("                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),")
    print("                     t = 1, l = 1, r = ncol(g))")
    print("pdf('{}', height = {}, width = 15)".format(
        out_file, round(table_row * 0.3)))
    print("grid.draw(g)")
    print("dev.off()")


def peaks_num_stats():
    peaks_dir = PROJECT_DIR + '/clam/peaks/'
    fw = open(OUT_DIR + 'peaks.counts.summary.txt', 'w')
    fw.write('sample\tpeak_num\n')
    for samples in COMPARISON_LIST:
        peaks_num = len(read_file(
            peaks_dir + 'peaks-{}__{}/narrow_peak.unique.bed'.format(samples[0], samples[1])))
        fw.write('{}\t{}\n'.format(
            NAME_SEP.join(samples[0].split(NAME_SEP)[:-1]), peaks_num))
    fw.close()


def load_theme():
    print( 'library(ggplot2)')
    print( 'library(ggthemes)')
    for line in read_file(cur_dir + '/theme_pub.R'):
        print( line.strip())


def plot_peaks_num(out_file):
    load_theme()
    print("data <- read.table('{}', header = T , sep = '\\t')".format(
        OUT_DIR + 'peaks.counts.summary.txt'))
    print("p1 <- ggplot(data, aes(x = sample, y = peak_num)) + geom_bar(stat='identity', aes(fill = peak_num)) + scale_fill_Publication() +theme_Publication() + ")
    print("  theme(axis.text.x = element_text(angle = 45, hjust = 1),  legend.position = 'none',")
    print("        strip.text = element_text(size=15), axis.title.x = element_blank(), legend.title = element_blank()) + ")
    print("  ylab('m6A Peak number') + scale_y_continuous(expand = c(0, 0)) + scale_fill_gradient2()")
    print("ggsave('{}', height = 5, width = 20)".format(out_file))


def get_motif_rScript(motif_file, title, index):
    pwm = '\n'.join(read_file(motif_file))
    pssm = np.array([[float(_x) for _x in x.split('\t')]
                     for x in pwm.strip().split('\n')[1:]]).T
    pssm_flatten = pssm.flatten()
    seq_len = len(pssm[0])
    print('seq_profile = c(' + ','.join([str(x) for x in pssm_flatten]) + ')')
    print('seq_matrix = matrix(seq_profile, 4, {}, byrow = T)'.format(seq_len))
    print("rownames(seq_matrix) <- c('A', 'C', 'G', 'U')")
    print("library(ggseqlogo)")
    print("library(ggplot2)")
    print("p{} <- ggplot() + geom_logo(seq_matrix, method = 'prob') + theme_logo() + theme(axis.text.x = element_blank(), ".format(
        index))
    print("       panel.spacing = unit(0.5, 'lines'),")
    print("       axis.text.y = element_blank(), ")
    print("       axis.title.y = element_blank(),")
    print("       axis.title.x = element_blank(), ")
    print("       plot.title = element_text(hjust = 0.5, size=30, face='bold'),")
    print("       legend.position = 'none') + ggtitle('{}')".format(title))


def get_dist_topology(outfn, index):
    print("data_fn='{}'".format(outfn))
    print("data = read.table(data_fn)")
    print("tot = nrow(data)")
    print("sep1 = tot/3+1")
    print("sep2 = tot/3*2+1")
    print("library(ggplot2)")
    print("colnames(data) = c('Peak_Tx', 'Peak_Tx_WithPeak', 'NumPeak')")
    print("data$coord = seq(1, nrow(data))")
    print("p{} <- ggplot(data, aes(x=coord, y=Peak_Tx_WithPeak, fill=NumPeak, colour=NumPeak)) + geom_point(colour='#386cb0', size=2) +".format(
        index))
    print("  geom_line() +")
    print("  guides(fill=F, colour=F) +")
    print("  theme_bw() +")
    print("  geom_vline(xintercept=c(sep1, sep2), linetype='dashed', size = 1.2) +")
    print("  annotate('text', label=c(\"5'UTR\", \"CDS\", \"3'UTR\"), size = 7, ")
    print(
        "           x=c( sep1/2, (sep1+sep2)/2, sep2+sep1/2 ), y=max(data[,2])*0.9) +  theme(axis.text.x = element_blank(), ")
    print("                                                                                    panel.spacing = unit(0.5, 'lines'),")
    print("                                                                                    panel.grid.major = element_blank(), ")
    print("                                                                                    panel.grid.minor = element_blank(),")
    print("                                                                                    axis.title.x = element_blank(), ")
    print("                                                                                    text = element_text(size=15),")
    print("                                                                                    axis.ticks.x = element_blank(),")
    print("                                                                                    legend.position = 'none') +")
    print("  scale_x_continuous(expand=c(0,0)) + ylab('Number of m6A peaks per transcript')")


def get_intersection_pie(inpfn, index):
    inp_list = read_file(inpfn)[1:]
    labels = []
    names = []
    for _inp in inp_list:
        sp = _inp.strip().split('\t')
        labels.append(str(round(float(sp[-1]) * 100, 1)) + '%')
        names.append(sp[0])
    labels[-1] = names[-1] + ' ' + labels[-1]
    label_info = ', '.join(["'{}'".format(label) for label in labels[::-1]])
    print("library(ggrepel)")
    print("data <- read.table('{}', header = T, sep = '\\t', quote=\"\\\"\")".format(inpfn))
    print("data$category <- factor(data$category, levels = c('3\\'UTR', '5\\'UTR', 'CDS', 'Other exon', 'Intron'))")
    print("data <- data[c(5, 4, 3, 2, 1),]")
    print(
        "data$pos = (cumsum(c(0, data$norm_count)) + c(data$norm_count /2, .01))[1:nrow(data)]")
    print("data$x <- c(1.5, 1.2, 1, 1, 1)")
    print("data$label <- c({})".format(label_info))
    print("data$nudge_x <- c(0.3, 0, 0.0, 0.0, 0.0)")
    print("data$nudge_y <- c(-0.05, 0, 0, 0, 0)")
    print('data$label_color <- c("black","white","white","white","white")')
    print('data$color <- rev(c("#7fc97f","#e79c47","#ef3b2c","#386cb0","#662506"))')
    print("p{} <- ggplot(data, aes(\"\", norm_count, fill = category)) + ".format(index))
    print("  geom_bar(width = 1, size = 1, color = 'white', stat = 'identity') +")
    print('  coord_polar("y") + ')
    print(
        "  geom_text_repel(aes(x = x, y = pos, label = label), color = data$label_color,")
    print("                  nudge_x = data$nudge_x,")
    print("                  nudge_y = data$nudge_y,")
    print("                  segment.size = .7, ")
    print("                  fontface = 'bold',")
    print("                  size = 8,")
    print("                  show.legend = FALSE) +")
    print("  scale_fill_manual(values = data$color) +")
    print("  theme_classic() +")
    print("  theme(axis.line = element_blank(),")
    print("        axis.text = element_blank(),")
    print("        axis.ticks = element_blank(),")
    print("        legend.key = element_rect(colour = NA),")
    print("        legend.position = c(0.5, 0.03),")
    print("        legend.direction = 'horizontal',")
    print("        legend.key.size= unit(0.2, 'cm'),")
    print("        legend.title = element_blank(),")
    print("        legend.spacing.x = unit(0.2, 'cm'),")
    print("        legend.key.width =unit(1,'line'),")
    print("        legend.text=element_text(size=20),")
    print("        axis.title.x = element_blank(),")
    print("        axis.title.y = element_blank(),")
    print("        plot.title = element_text(hjust = 0.5, color = '#666666')) ")


def draw_peaks_dis_figures(samples):
    motif_file = PROJECT_DIR + \
        '/homer/{}__{}/clam_unique/homerResults/motif1.motif'.format(
            samples[0], samples[1])
    peak = PROJECT_DIR + \
        '/clam/peaks/peaks-{}__{}/narrow_peak.unique.bed'.format(
            samples[0], samples[1])
    dist_data = PROJECT_DIR + \
        '/topology/{}__{}/clam_unique/dist.data'.format(
            samples[0], samples[1])
    title = NAME_SEP.join(samples[0].split(NAME_SEP)[:-1])
    gtf_bed_dir = '{}/scripts/peakComposition/{}'.format(SNAKE_SCRIPT_DIR, GENOME)
    get_motif_rScript(motif_file, title, 1)
    get_dist_topology(dist_data, 2)
    outfn = '{}{}data.pie.txt'.format(TMP_DIR, title)
    intersect_gtf_regions(peak, outfn, gtf_bed_dir)
    get_intersection_pie(outfn, 3)
    pdf_output = OUT_DIR + '{}.pdf'.format(title)
    print('library(ggpubr)')
    print('pdf(file="{}", width=7, height=11)'.format(pdf_output))
    print('ggarrange(p1, ')
    print('          p2,')
    print('          p3,')
    print('          heights = c(1, 1.3, 2.5),')
    print('          ncol = 1, nrow = 3, align = \'v\')')
    print('dev.off()')
    return pdf_output


def get_regulator_exp():
    t2g = {}
    for line in open(T2G):
        info = line.strip().split('\t')
        t2g[info[0]] = info[1]
    ensg2gn = {}
    for line in open(ENSG2GN):
        info = line.strip().split('\t')
        ensg2gn[info[0]] = info[1]

    gene_exp = defaultdict(lambda: defaultdict(float))
    for sample in SAMPLE_TYPE_DICT:
        if SAMPLE_TYPE_DICT[sample] == 'con':
            for line in open(os.path.join(PROJECT_DIR, 'kallisto', sample, 'abundance.tsv')):
                info = line.strip().split('\t')
                info[0] = info[0].split('.')[0]
                if info[0] in t2g and t2g[info[0]] in REGULATORS:
                    gene_exp[ensg2gn[t2g[info[0]]]
                             ][sample] = gene_exp[ensg2gn[t2g[info[0]]]][sample]+float(info[-1])
    data_file = open(os.path.join(
        PROJECT_DIR, 'tmp', 'regulator_exp.txt'), 'w')
    first_line = True
    header = "Gene"
    for gene in gene_exp:
        to_write = [gene]
        for sample in SAMPLE_TYPE_DICT:
            if SAMPLE_TYPE_DICT[sample] == 'con':
                if first_line:
                    header = header+'\t'+sample
                to_write.append(str(gene_exp[gene][sample]))
        if first_line:
            first_line = False
            data_file.write(header+'\n')
        data_file.write('\t'.join(to_write)+'\n')
    data_file.flush()
    out_pdf = os.path.join(PROJECT_DIR, 'tmp', 'regulators_exp.pdf')
    print("library(ComplexHeatmap)")
    print("library(circlize)")
    print("library(RColorBrewer)")
    print("data_table<-read.csv('{}',sep = '\\t',header = T,row.names = 1)".format(
        os.path.join(PROJECT_DIR, 'tmp', 'regulator_exp.txt')))
    print("color_fun=colorRamp2(c(-2, 0, 2), c('cyan', 'white', 'violet'))")
    print("data_table = t(data_table)")
    print("data_table <- scale(data_table[, ])")
    print("data_table <- t(data_table)")
    print("pdf('{}',width=20)".format(out_pdf))
    print("Heatmap(data_table,cluster_columns = T,cluster_rows = T,col=color_fun, name='z-score',show_row_names = T,show_column_names = T)")
    print("dev.off()")
    return out_pdf


def sig_gene(file_in, from_deseq2, log2fc_cutoff=1, fdr_cutoff=0.01):
    df = pd.read_csv(file_in, sep='\t')
    df = df.dropna()
    if from_deseq2:
        sig_genes = df[
            (df.padj < fdr_cutoff) &
            ((df.log2FoldChange > log2fc_cutoff) |
             (df.log2FoldChange < -log2fc_cutoff))].index
    else:
        sig_genes = df[
            (df.log2FoldChange > log2fc_cutoff |
             (df.log2FoldChange < -log2fc_cutoff))].index
    return sig_genes


def enrichment_of_sig_gene(fdr_cutoff=0.01, log2fc_cutoff=1):
    data_list = []
    list_names = []
    for group in DEG_GROUP:
        from_deseq2 = True
        for samples in DEG_GROUP[group]:
            if len(samples) == 1:
                from_deseq2 = False
                break
        data_list.append(sig_gene(os.path.join(
            PROJECT_DIR, 'deg', group+'.deg.txt'), from_deseq2, log2fc_cutoff, fdr_cutoff).tolist())
        list_names.append(group)
    ret = hyper_geometry.count_numbers(
        data_list, list_names, id_type='ENSG', result_dir=os.path.join(PROJECT_DIR, 'tmp'), width=12, height=10, text_width=100)
    if ret == 0:
        return os.path.join(PROJECT_DIR, 'tmp', 'plot_result_GO.pdf')
    else:
        return ''


def plot_trip_seq():
    script="""\
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
data_table<-read.csv('{}/sig_peaks.csv_',sep = '\t',header = T,row.names = 1)
color_fun=colorRamp2(c(-2, 0, 2), c('cyan', 'white', 'violet'))
data_table = t(data_table)
data_table <- scale(data_table[, ])
data_table <- t(data_table)
pdf('{}',width=5)
Heatmap(data_table,cluster_columns = F,cluster_rows = T,col=color_fun, name='z-score',show_row_names = F,show_column_names = T)
dev.off()
""".format(PROJECT_DIR+'/parsed_peaks',PROJECT_DIR+'/summary/trip_sig.pdf')
    print( script)
    script = """\
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
data_table<-read.csv('{}/peak_intensity.csv',sep = '\t',header = T,row.names = 1)
color_fun=colorRamp2(c(-2, 0, 2), c('cyan', 'white', 'violet'))
data_table = t(data_table)
data_table <- scale(data_table[, ])
data_table <- t(data_table)
pdf('{}',width=5)
Heatmap(data_table,cluster_columns = F,cluster_rows = T,col=color_fun, name='z-score',show_row_names = F,show_column_names = T)
dev.off()
""".format(PROJECT_DIR+'/parsed_peaks', PROJECT_DIR+'/summary/peak_intensity.pdf')
    print( script)
    return PROJECT_DIR+'/summary/peak_intensity.pdf',PROJECT_DIR+'/summary/trip_sig.pdf'


def get_pdf_from_sample():
    if not os.path.isdir(TMP_DIR):
        os.mkdir(TMP_DIR)
    rscript_out = TMP_DIR + 'generate_figure.R'
    pdf_output_files = []
    sys.stdout = open(rscript_out, "w+")
    star_mapping_stats()
    table_generate(OUT_DIR + 'mapping_stats.counts.summary.txt',
                   OUT_DIR + 'mapping.stats.table.pdf')
    pdf_output_files.append(OUT_DIR + 'mapping.stats.table.pdf')
    peaks_num_stats()
    plot_peaks_num(OUT_DIR + 'peaks.stats.pdf')
    pdf_output_files.append(OUT_DIR + 'peaks.stats.pdf')

    regulator_exp = get_regulator_exp()
    pdf_output_files.append(regulator_exp)

    if len(DEG_GROUP) > 0:
        go_enrich = enrichment_of_sig_gene()
        if len(go_enrich) > 1:
            pdf_output_files.append(go_enrich)

    dict_comparison = {}
    for samples in COMPARISON_LIST:
        sample_name = NAME_SEP.join(samples[0].split(NAME_SEP)[:-1])
        dict_comparison[sample_name] = samples
    for sample_name in sorted(dict_comparison.keys()):
        samples = dict_comparison[sample_name]
        pdf_output = draw_peaks_dis_figures(samples)
        pdf_output_files.append(pdf_output)

    if IS_TRIPSEQ:
        if os.path.isfile(PROJECT_DIR+'/parsed_peaks/sig_peaks.csv'):
            pdf_output_files.extend(plot_trip_seq())
    sys.stdout.flush()
    sys.stdout.close()
    call('Rscript {}'.format(rscript_out), shell=True)
    sys.stdout = sys.__stdout__
    merger = PdfFileMerger()
    print( pdf_output_files)

    for filename in pdf_output_files:
        merger.append(PdfFileReader(open(filename, 'rb')))
    merger.write(OUT_DIR + "merged_summary.pdf")


# peaks_num_stats()
get_pdf_from_sample()

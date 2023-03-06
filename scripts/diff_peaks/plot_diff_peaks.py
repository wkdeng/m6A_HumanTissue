'''
read in peaks_intensity.csv and differential_sites_log1.5.csv, plot differential peaks
DWK
6.3.2019
'''

import os
import yaml
import sys
import pandas as pd
import fileinput

try:
    config = yaml.safe_load(open(sys.argv[1]))
except:
    sys.exit('command insufficient {config files}')
    # config = yaml.safe_load(open(r'H:\yi_lab\m6a\src\mouse_m6a\config.yaml'))

PROJECT_DIR = config['PROJECT_DIR']
IP_SUFFIX = config['ip_suffix']
CON_SUFFIX = config['con_suffix']
DEG_GROUP = dict(config['deg_group'])


def plot_(tissue_specific=False):
    if not os.path.exists(PROJECT_DIR + '/diff_peaks/differential_sites_log1.5' + (
    '-tissue_specific' if tissue_specific else '') + '.processed.csv'):
        processed_diff_file = open(PROJECT_DIR + '/diff_peaks/differential_sites_log1.5' + (
            '-tissue_specific' if tissue_specific else '') + '.processed.csv', 'w')
        first_line = True
        for line in fileinput.input(PROJECT_DIR + '/diff_peaks/differential_sites_log1.5' + (
        '-tissue_specific' if tissue_specific else '') + '.csv'):
            if first_line:
                processed_diff_file.write(line)
                first_line = False
                continue

            info = line.strip().split('\t')
            info[0] = info[0].replace('-', ':')
            processed_diff_file.write('\t'.join(info) + '\n')
        processed_diff_file.flush()
        processed_diff_file.close()

    peak_intensitys = pd.read_csv(PROJECT_DIR + '/parsed_peaks/peak_intensity.csv', sep='\t', index_col=0)
    diff_peaks = pd.read_csv(PROJECT_DIR + '/diff_peaks/differential_sites_log1.5' + (
        '-tissue_specific' if tissue_specific else '') + '.processed.csv', sep='\t', index_col=0)
    # peak_intensitys = pd.read_csv(r'H:\yi_lab\m6a\src\scripts\diff_peaks\peak_intensity.csv', sep='\t', index_col=0)
    # diff_peaks = pd.read_csv(r'H:\yi_lab\m6a\src\scripts\diff_peaks\differential_sites_log1.5.processed.csv', sep='\t',
    #                          index_col=0)

    for group in DEG_GROUP:
        if group == 'All_Sample':
            continue
        group_diff_peaks = diff_peaks.loc[diff_peaks[group] != 0].index

        all_sample = []
        [all_sample.append(x[0:-len(CON_SUFFIX)] if x.endswith(CON_SUFFIX) else x) for x in DEG_GROUP[group][0]]
        [all_sample.append(x[0:-len(CON_SUFFIX)] if x.endswith(CON_SUFFIX) else x) for x in DEG_GROUP[group][1]]
        df_group = peak_intensitys.loc[group_diff_peaks, all_sample]
        df_group.to_csv(PROJECT_DIR + '/diff_peaks/{group}.diff_peaks'.format(group=group) + (
            '-tissue_specific' if tissue_specific else '') + '.csv', sep='\t')

        script_file = PROJECT_DIR + '/diff_peaks/{group}.diff_peak'.format(group=group) + (
            '-tissue_specific' if tissue_specific else '') + '.R'
        r_script = open(script_file, 'w')
        reps = "rep('1', 10),rep('2', 10)"
        if group in DEG_GROUP:
            reps = "rep('1', " + str(len(DEG_GROUP[group][0])) + "),rep('2', " + str(
                len(DEG_GROUP[group][1])) + ")"
        r_script.write("""library(ComplexHeatmap)
    library(circlize)
    library(RColorBrewer)
    y <- read.table('{file_in}', row.names = 1, header = T, sep = '\\t', check.names = FALSE)
    
    y = t(scale(t(y)))
    df <- data.frame(sample = c({reps}))
    col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    col_fun(seq(-4, 4))
    top_annotation <- HeatmapAnnotation(df = df, col = list(sample = c('1' = "#386cb0", '2' = "#ef3b2c",'3'='#fdb462','4'='#7fc97f')), which = 'column', show_legend = FALSE)
    ha1 = Heatmap(y, col=col_fun, clustering_distance_rows = 'spearman',clustering_distance_columns = 'pearson', show_row_names =FALSE, show_row_dend=FALSE, top_annotation = top_annotation, name = 'Z-score')
    pdf('{file_out}')
    draw(ha1)
    dev.off()
    """.format(file_in=PROJECT_DIR + '/diff_peaks/{group}.diff_peaks'.format(group=group) + (
            '-tissue_specific' if tissue_specific else '') + '.csv',
               file_out=PROJECT_DIR + '/diff_peaks/{group}.diff_peaks'.format(group=group) + (
                   '-tissue_specific' if tissue_specific else '') + '.pdf', reps=reps))
        r_script.flush()
        os.system('Rscript ' + script_file)
        os.remove(script_file)


plot_()
plot_(True)

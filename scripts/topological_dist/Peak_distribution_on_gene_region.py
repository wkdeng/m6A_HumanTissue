# /**
#  * @author [Wankun Deng]
#  * @email [dengwankun@hotmail.com]
#  * @create date 2019-09-24 15:48:58
#  * @modify date 2019-09-24 15:48:58
#  * @desc [Modified to require >= 50% of overlap and 300nt up/down stream proximal]
#  */
import sys
import yaml
import os
from pybedtools import BedTool
from multiprocessing import Pool
from multiprocessing import Manager
from functools import partial
import subprocess

config = yaml.safe_load(open(sys.argv[1]))
TMP_DIR = config['genomic_region']['tmp_dir']
out_array = {'K562': [open(TMP_DIR + '/all_peak_postion_K562.txt', 'w'),
                      open(TMP_DIR + '/alu_peak_postion_K562.txt', 'w')], 'HepG2': [open(TMP_DIR + '/all_peak_postion_HepG2.txt', 'w'),
                                                                                       open(TMP_DIR + '/alu_peak_postion_HepG2.txt', 'w')]}
for pre in out_array:
    for fp in out_array[pre]:
        fp.write('Group\tPosition\tPercentage\n')
        fp.flush()

PROJECT_DIR = config['PROJECT_DIR']

GENOME = config['genome']
ALU_FILE = config[GENOME]['alu_repetitive']

for repeat_file in config[GENOME]['repetitive']:
    if repeat_file.endswith(ALU_FILE):
        ALU_FILE = repeat_file
        break

gene_region = BedTool(config['genomic_region']['gene'])
non_coding_region = BedTool(config['genomic_region']['non_coding'])
cds_region = BedTool(config['genomic_region']['cds'])
utr5_region = BedTool(config['genomic_region']['5UTR'])
utr3_region = BedTool(config['genomic_region']['3UTR'])
up200_region = BedTool(config['genomic_region']['up_proximal'])
down200_region = BedTool(config['genomic_region']['down_proximal'])

print_index = ['Non_coding', 'CDS', 'UTR5', 'UTR3', 'Up_stream_300nt',
               'Down_stream_300nt', 'Distal_intronic',  'Alu_non_coding',
               'Alu_CDS', 'Alu_5UTR', 'Alu_3UTR', 'Alu_up_stream_300nt', 'Alu_down_stream_300nt', 'Alu_distal_intronic']

def count_group(lock, group):
    # non_coding count, cds count, 5UTR count, 3UTR count, up stream 200 proximal count, down stream 200 proximal count, distal intronic count.

    group_count = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    k = 0

    tmp_fp = TMP_DIR + '/'+group+'.all.mid.bed'
    tmp_file = open(tmp_fp, 'w')
    for line in open(PROJECT_DIR+'/clam/peaks/'+group+'/narrow_peak.combined.bed'):
        info = line.strip().split('\t')
        start = int(info[1])
        end = int(info[2])
        mid = str(int((start + end) / 2))
        tmp_file.write(
            '\t'.join([info[0], mid, mid, info[3], info[4], info[5]]) + '\n')
    tmp_file.flush()
    tmp_file.close()
    all_peak_bed = BedTool(tmp_fp)

    alu_peak = BedTool(PROJECT_DIR + '/clam/peaks/' + group +
                       '/narrow_peak.combined.bed').intersect(BedTool(ALU_FILE), f=0.5, s=True, wa=True, u=True)
    alu_peak_bed = ''
    if len(alu_peak) > 0:
        alu_peak.saveas(TMP_DIR + '/'+group+'.alu.full.bed')
        tmp_fp2 = TMP_DIR + '/'+group+'.alu.mid.bed'
        tmp_file2 = open(tmp_fp2, 'w')
        for line in open(TMP_DIR + '/'+group+'.alu.full.bed'):
            info = line.strip().split('\t')
            start = int(info[1])
            end = int(info[2])
            mid = str(int((start + end) / 2))
            tmp_file2.write(
                '\t'.join([info[0], mid, mid, info[3], info[4], info[5]]) + '\n')
        tmp_file2.flush()
        tmp_file2.close()
        alu_peak_bed = BedTool(tmp_fp2)

    for peak_bed in [all_peak_bed, alu_peak_bed]:
        # all_peak_file =
        # alu_peak_file =

        # # count intergenic
        # n1 = len(peak_bed)
        # subtracted_intergenic = peak_bed.subtract(gene_region, A=True, s=True)
        # subtracted_intergenic.saveas(TMP_DIR + '/' + group + '.intergenic.bed')
        # group_count.append(len(subtracted_intergenic))
        # peak_bed = peak_bed.intersect(gene_region, wa=True, u=True, s=True)
        # peak_bed.saveas(TMP_DIR + '/' + group + '.intergenic.removed.bed')
        # n2 = len(subtracted_intergenic)
        # n3 = len(peak_bed)

        # count non_coding
        if len(peak_bed) <= 0:
            continue
        subtracted_non_coding = peak_bed.intersect(
            non_coding_region, f=0.5, wa=True, u=True, s=True)
        subtracted_non_coding.saveas(TMP_DIR + '/' + group + '.noncoding.bed')
        group_count[0+k*7] = (len(subtracted_non_coding))
        peak_bed = peak_bed.subtract(non_coding_region, f=0.5, A=True, s=True)

        # count cds
        subtracted_cds = peak_bed.intersect(
            cds_region, f=0.5, wa=True, u=True, s=True)
        subtracted_cds.saveas(TMP_DIR + '/' + group + '.cds.bed')
        group_count[1+k*7] = (len(subtracted_cds))
        peak_bed = peak_bed.subtract(cds_region, f=0.5, A=True, s=True)

        # count 5UTR
        subtracted_5utr = peak_bed.intersect(
            utr5_region, f=0.5, wa=True, u=True, s=True)
        subtracted_5utr.saveas(TMP_DIR + '/' + group + '.5utr.bed')
        group_count[2+k*7] = (len(subtracted_5utr))
        peak_bed = peak_bed.subtract(utr5_region, f=0.5, A=True, s=True)

        # count 3UTR
        subtracted_3utr = peak_bed.intersect(
            utr3_region, f=0.5, wa=True, u=True, s=True)
        subtracted_3utr.saveas(TMP_DIR + '/' + group + '.3utr.bed')
        group_count[3+k*7] = (len(subtracted_3utr))
        peak_bed = peak_bed.subtract(utr3_region, f=0.5, A=True, s=True)

        # count up stream proximal 200
        subtracted_up200 = peak_bed.intersect(
            up200_region,f=0.5,  wa=True, u=True, s=True)
        subtracted_up200.saveas(TMP_DIR + '/' + group + '.up200.bed')
        group_count[4+k*7] = (len(subtracted_up200))
        peak_bed = peak_bed.subtract(up200_region, f=0.5, A=True, s=True)

        # count down stream proximal 200
        subtracted_down200 = peak_bed.intersect(
            down200_region, f=0.5, wa=True, u=True, s=True)
        subtracted_down200.saveas(TMP_DIR + '/' + group + '.down200.bed')
        group_count[5+k*7] = (len(subtracted_down200))
        peak_bed = peak_bed.subtract(down200_region, f=0.5, A=True, s=True)

        # the rest was treated as intergenic region peaks
        group_count[6 + k * 7] = (len(peak_bed))
        k = k + 1
    sum_up_1 = float(sum(group_count[0:7]))
    sum_up_2 = float(sum(group_count[7:]))

    for i in range(7):
        group_count[i] = group_count[i] / sum_up_1
        if sum_up_2 > 0:
            group_count[i+7] = group_count[i+7]/sum_up_2

    lock.acquire()
    splitted_group=group.split('_')
    for j in range(2):
        for i in range(8):
            out_array[splitted_group[0]][j].write(splitted_group[1]+'\t'+print_index[i+j*7] +
                               '\t' + str(group_count[i + j * 7]) + '\n')
        out_array[splitted_group[0]][j].flush()
    # print group+'\t'+'\t'.join([str(x) for x in group_count])
    lock.release()


# print 'Group\t'+'\t'.join(['Non_coding', 'CDS', '5\'UTR', '3\'UTR', 'Up_stream_300nt',
#                            'Down_stream_300nt', 'Distal_intronic', 'Intergenic', 'Alu_non_coding', 'Alu_CDS', 'Alu_5UTR', 'Alu_3UTR', 'Alu_up_stream_200_proximal', 'Alu_down_stream_200_proximal', 'Alu_distal_intronic', 'Alu_intergenic'])

pool = Pool(processes=20)
manager = Manager()
lock = manager.Lock()
func = partial(count_group, lock)
pool.map(func, config['groups'].keys())
pool.close()
pool.join()
for pre in out_array:
    for fp in out_array[pre]:
        fp.flush()
        fp.close()

rscript = open(TMP_DIR + '/plot_region_dist_org_order.R', 'w')
rscript.write('''
library(ggplot2)
library(ggpubr)
c0 = c('Non_coding', 'CDS', 'UTR5', 'UTR3',
       'Distal_intronic', 'Down_stream_300nt', 'Up_stream_300nt')
c1 = c('Alu_non_coding',
       'Alu_CDS', 'Alu_5UTR', 'Alu_3UTR', 'Alu_distal_intronic', 'Alu_down_stream_300nt', 'Alu_up_stream_300nt')


all_k562_order=c('DDX21', 'DDX3X','DGCR8','ILF3','DROSHA','FMR1','FXR1','HNRNPU','DHX30','SUPV3L1','SF3B4','U2AF2','SMNDC1','U2AF1','GPKOW','EFTUD2','AQR','BUD13','XRN2','KHSRP','SF3B1','DDX42','PRPF8','EWSR1','RBM22','HNRNPA1',
                 'GTF2F1','SSB','DDX52','RBFOX2','ZRANB2','AATF','TAF15','DDX51','NONO','PHF6','CSTF2T','AKAP8L','ABCF1',
                 'FASTKD2','HNRNPC','HNRNPUL1','UTP3','FUS','NCBP2','WDR3','SAFB','LSM11','FTO','LARP7','GEMIN5','NPM1','WDR43','AARS',
                 'SRSF1','TIA1','PPIL4','DDX24','YWHAG','CPSF6','EIF4G2','PUS1','NSUN2','WRN','SLTM','QKI','HNRNPK','UTP18','XRCC6','ZC3H8',
                 'SERBP1','SDAD1','HNRNPM','KHDRBS1','CPEB4','TRA2A','GNL3','SLBP','ZNF622','MTPAP','TROVE2','UCHL5','TBRG4','PCBP1','HNRNPL','RPS11',
                 'GRWD1','RBM15','TARDBP','SRSF7','LIN28B','EIF3G','PUM1','AGGF1','RPS3','NIPBL','ZNF800','METAP2','SND1','HLTF','LARP4',
                 'PTBP1','SBDS','YBX3','DDX6','APOBEC3C','SAFB2','FXR2','IGF2BP2','NOLC1','DDX55','ZC3H11A','MATR3','IGF2BP1','FAM120A','AKAP1',
                 'EXOSC5','PABPC4','PUM2','UPF1')
alu_k562_order=c('DDX21','DDX3X','DGCR8','ILF3','DROSHA','FMR1','FXR1','HNRNPU','DHX30','SUPV3L1','SERBP1','SF3B4','GRWD1','U2AF2','PCBP1','ABCF1','KHSRP','GEMIN5','U2AF1','BUD13','GPKOW','NOLC1','XRN2','TRA2A',
                 'HNRNPK','SMNDC1','SRSF1','RBM15','YWHAG','MTPAP','EWSR1','SLBP','FXR2','DDX24','RBFOX2','DDX55','AQR','AATF','METAP2',
                 'HNRNPA1','TROVE2','LARP7','SAFB','HNRNPC','ZNF622','FUS','SSB','GTF2F1','TARDBP','LSM11','ZC3H11A','WDR43','KHDRBS1','UCHL5',
                 'EFTUD2','TAF15','DDX52','HNRNPM','LARP4','PRPF8','LIN28B','XRCC6','AARS','HNRNPL','IGF2BP1','RPS3','CPEB4','NONO','AGGF1',
                 'ZRANB2','RBM22','WRN','NSUN2','UTP3','TIA1','SLTM','SAFB2','CSTF2T','PUM2','IGF2BP2','NIPBL','EIF4G2','PABPC4','SF3B1','UPF1','YBX3',
                 'AKAP1','FAM120A','PUS1','QKI','RPS11','EIF3G','SRSF7','DDX42','PTBP1','UTP18','FASTKD2','DDX6','MATR3','SBDS','APOBEC3C','PUM1',
                 'SND1','NPM1','HLTF','PPIL4','HNRNPUL1','PHF6','CPSF6','WDR3','TBRG4','ZNF800','DDX51','ZC3H8','NCBP2','GNL3','EXOSC5','AKAP8L','FTO','SDAD1')
all_hepg2_order=c('DDX3X','DGCR8','ILF3','DROSHA','STAU2','DHX30','SUPV3L1','SF3B4','SF3A3','U2AF2','AQR','SMNDC1','BUD13','U2AF1','EFTUD2','RBM5','XRN2','KHSRP','CDC40','RBM22','RBFOX2','HNRNPU',
                  'PRPF8','HNRNPA1','FKBP4','NKRF','GTF2F1','TAF15','CSTF2T','NIP7','FAM120A','PRPF4','TIAL1','HNRNPUL1','FUS','SSB','BCCIP',
                  'TIA1','XPO5','HNRNPK','DDX59','XRCC6','SLTM','DDX52','LARP7','UCHL5','TBRG4','FTO','GRSF1','TROVE2','SRSF1','FASTKD2','AGGF1',
                  'NCBP2','UTP18','SRSF7','RBM15','PPIG','LSM11','BCLAF1','QKI','PCBP2','POLR2G','PCBP1','TRA2A','EXOSC5','NOLC1','GRWD1','DKC1','YBX3',
                  'WDR43','SND1','SRSF9','SDAD1','EIF3D','HNRNPC','SAFB','ZNF800','CSTF2','SFPQ','LIN28B','NOL12','PABPN1','FXR2',
                  'DDX6','PTBP1','ZC3H11A','LARP4','RPS3','HLTF','HNRNPL','G3BP1','IGF2BP3','HNRNPM','EIF3H','DDX55','IGF2BP1','AKAP1','MATR3','SUGP2','SUB1','FUBP3','UPF1')
alu_hepg2_order=c('DDX3X','DGCR8','ILF3','DROSHA','HNRNPU','STAU2','DHX30','SUPV3L1','NIP7','SF3B4','SF3A3','BUD13','RBM15','G3BP1','SMNDC1','AQR','U2AF2','SSB','KHSRP','SFPQ','RBM5','HNRNPUL1','EFTUD2',
                  'EIF3D','PPIG','GRSF1','BCLAF1','FASTKD2','XRN2','RPS3','FUS','TROVE2','NKRF','GRWD1','U2AF1','HNRNPA1','XPO5','PRPF4',
                  'DDX52','HNRNPM','PRPF8','SLTM','SUB1','TAF15','LIN28B','NCBP2','AGGF1','CSTF2','RBFOX2','NOLC1','CDC40','XRCC6','HNRNPK',
                  'RBM22','QKI','HNRNPC','SAFB','TIA1','TBRG4','FAM120A','HNRNPL','ZC3H11A','DKC1','CSTF2T','GTF2F1','DDX59','YBX3','SUGP2',
                  'PCBP2','HLTF','PTBP1','UPF1','IGF2BP1','PABPN1','LARP7','TRA2A','MATR3','NOL12','LARP4','TIAL1','AKAP1','IGF2BP3','WDR43','SRSF7',
                  'FTO','FKBP4','UCHL5','SRSF1','FUBP3','EIF3H','SND1','SDAD1','POLR2G','SRSF9','DDX55','ZNF800','PCBP1','DDX6','EXOSC5',
                  'BCCIP','LSM11','FXR2','UTP18')
                  
data <- read.csv('{}', sep='\t', header=T)
data$Group<-factor(data$Group,levels =all_k562_order )
p1=ggplot(data, aes(y=Percentage, x=Group, fill=factor(Position, levels=c0)))+geom_bar(position='stack', stat='identity') + theme(axis.text.x=element_text(angle=90)) +\
ggtitle('{}')+ theme(plot.title = element_text(face = "bold",hjust=0.5))+ guides(fill=guide_legend(title="Position"))

data <- read.csv('{}', sep='\t', header=T)
data$Group<-factor(data$Group,levels =alu_k562_order )
p2=ggplot(data, aes(y=Percentage, x=Group, fill=factor(Position, levels=c1)))+geom_bar(position='stack', stat='identity') + theme(axis.text.x=element_text(angle=90)) +\
ggtitle('{}')+ theme(plot.title = element_text(face = "bold",hjust=0.5))+ guides(fill=guide_legend(title="Position"))
pdf('{}',width=17)
ggarrange(p1,p2,heights = c(1, 1),ncol = 1, nrow = 2, align = 'v')
dev.off()

data <- read.csv('{}', sep='\t', header=T)
data$Group<-factor(data$Group,levels =all_hepg2_order )
p3=ggplot(data, aes(y=Percentage, x=Group, fill=factor(Position, levels=c0)))+geom_bar(position='stack', stat='identity') + theme(axis.text.x=element_text(angle=90)) +\
ggtitle('{}')+ theme(plot.title = element_text(face = "bold",hjust=0.5))+ guides(fill=guide_legend(title="Position"))

data <- read.csv('{}', sep='\t', header=T)
data$Group<-factor(data$Group,levels =alu_hepg2_order )
p4=ggplot(data, aes(y=Percentage, x=Group, fill=factor(Position, levels=c1)))+geom_bar(position='stack', stat='identity') + theme(axis.text.x=element_text(angle=90)) +\
ggtitle('{}')+ theme(plot.title = element_text(face = "bold",hjust=0.5))+ guides(fill=guide_legend(title="Position"))
pdf('{}',width=17)
ggarrange(p3,p4,heights = c(1, 1),ncol = 1, nrow = 2, align = 'v')
dev.off()
'''.format(TMP_DIR+'/all_peak_postion_K562.txt', 'K562_All',
           TMP_DIR + '/alu_peak_postion_K562.txt', 'K562_Alu',
           TMP_DIR + '/K562_region_dist_org_order.pdf',
           TMP_DIR + '/all_peak_postion_HepG2.txt',  'HepG2_All',
           TMP_DIR + '/alu_peak_postion_HepG2.txt', 'HepG2_Alu',
           TMP_DIR + '/HepG2_region_dist_org_order.pdf',))
rscript.flush()
rscript.close()

subprocess.call('Rscript '+TMP_DIR +
                '/plot_region_dist_org_order.R', shell=True)


rscript = open(TMP_DIR + '/plot_region_dist_sort_by_order.R', 'w')
rscript.write('''
library(ggplot2)
library(ggpubr)
c0 = c('Non_coding', 'CDS', 'UTR5', 'UTR3',
       'Distal_intronic', 'Down_stream_300nt', 'Up_stream_300nt')
c1 = c('Alu_non_coding',
       'Alu_CDS', 'Alu_5UTR', 'Alu_3UTR', 'Alu_distal_intronic',  'Alu_down_stream_300nt', 'Alu_up_stream_300nt')


all_k562_order=c('SF3B4','U2AF2','SMNDC1','U2AF1','GPKOW','EFTUD2','AQR','BUD13','XRN2','KHSRP','SF3B1','DDX42','PRPF8','EWSR1','RBM22','HNRNPA1',
                 'GTF2F1','SSB','DROSHA','DGCR8','DDX52','RBFOX2','ZRANB2','AATF','TAF15','DDX51','HNRNPU','NONO','PHF6','CSTF2T','AKAP8L','ABCF1',
                 'FASTKD2','HNRNPC','HNRNPUL1','UTP3','FUS','DHX30','NCBP2','WDR3','SAFB','LSM11','FTO','LARP7','GEMIN5','NPM1','DDX21','WDR43','AARS',
                 'SRSF1','TIA1','PPIL4','DDX24','YWHAG','CPSF6','EIF4G2','PUS1','SUPV3L1','NSUN2','WRN','SLTM','QKI','HNRNPK','UTP18','XRCC6','ZC3H8',
                 'SERBP1','SDAD1','HNRNPM','KHDRBS1','CPEB4','TRA2A','GNL3','SLBP','ZNF622','MTPAP','TROVE2','UCHL5','TBRG4','PCBP1','HNRNPL','RPS11',
                 'GRWD1','RBM15','TARDBP','SRSF7','LIN28B','EIF3G','PUM1','AGGF1','RPS3','NIPBL','ZNF800','METAP2','ILF3','DDX3X','SND1','HLTF','LARP4',
                 'FMR1','PTBP1','SBDS','YBX3','DDX6','APOBEC3C','SAFB2','FXR2','IGF2BP2','NOLC1','DDX55','ZC3H11A','MATR3','IGF2BP1','FAM120A','AKAP1',
                 'EXOSC5','FXR1','PABPC4','PUM2','UPF1')
alu_k562_order=c('SERBP1','SF3B4','GRWD1','U2AF2','PCBP1','ABCF1','KHSRP','DDX3X','GEMIN5','U2AF1','BUD13','GPKOW','NOLC1','DHX30','XRN2','TRA2A',
                 'HNRNPK','SMNDC1','SRSF1','RBM15','YWHAG','MTPAP','EWSR1','SLBP','FXR2','DDX24','RBFOX2','DDX55','AQR','SUPV3L1','AATF','METAP2',
                 'HNRNPA1','TROVE2','LARP7','SAFB','HNRNPU','HNRNPC','ZNF622','FUS','SSB','GTF2F1','TARDBP','LSM11','ZC3H11A','WDR43','KHDRBS1','UCHL5',
                 'EFTUD2','TAF15','DDX52','HNRNPM','LARP4','PRPF8','LIN28B','XRCC6','AARS','HNRNPL','FMR1','IGF2BP1','RPS3','CPEB4','NONO','DGCR8','AGGF1',
                 'ZRANB2','RBM22','WRN','NSUN2','UTP3','TIA1','SLTM','SAFB2','CSTF2T','PUM2','IGF2BP2','NIPBL','EIF4G2','PABPC4','SF3B1','ILF3','UPF1','YBX3',
                 'AKAP1','FAM120A','PUS1','QKI','RPS11','EIF3G','SRSF7','DDX42','PTBP1','UTP18','FASTKD2','DDX6','MATR3','DDX21','SBDS','APOBEC3C','PUM1',
                 'SND1','NPM1','FXR1','DROSHA','HLTF','PPIL4','HNRNPUL1','PHF6','CPSF6','WDR3','TBRG4','ZNF800','DDX51','ZC3H8','NCBP2','GNL3','EXOSC5','AKAP8L','FTO','SDAD1')

all_hepg2_order=c('SF3B4','SF3A3','U2AF2','AQR','SMNDC1','BUD13','U2AF1','EFTUD2','RBM5','SUPV3L1','XRN2','KHSRP','CDC40','RBM22','DGCR8','RBFOX2',
                  'PRPF8','HNRNPA1','FKBP4','NKRF','GTF2F1','DROSHA','TAF15','CSTF2T','NIP7','FAM120A','PRPF4','TIAL1','HNRNPUL1','FUS','SSB','BCCIP',
                  'TIA1','XPO5','HNRNPK','DDX59','XRCC6','SLTM','DDX52','LARP7','UCHL5','TBRG4','FTO','GRSF1','HNRNPU','TROVE2','SRSF1','FASTKD2','AGGF1',
                  'NCBP2','UTP18','SRSF7','RBM15','PPIG','LSM11','BCLAF1','QKI','PCBP2','POLR2G','PCBP1','TRA2A','EXOSC5','NOLC1','GRWD1','DKC1','YBX3',
                  'DHX30','WDR43','SND1','SRSF9','SDAD1','EIF3D','HNRNPC','SAFB','ZNF800','CSTF2','DDX3X','SFPQ','LIN28B','NOL12','PABPN1','STAU2','FXR2',
                  'DDX6','PTBP1','ZC3H11A','ILF3','LARP4','RPS3','HLTF','HNRNPL','G3BP1','IGF2BP3','HNRNPM','EIF3H','DDX55','IGF2BP1','AKAP1','MATR3','SUGP2','SUB1','FUBP3','UPF1')
alu_hepg2_order=c('NIP7','SF3B4','SF3A3','BUD13','RBM15','G3BP1','SMNDC1','AQR','U2AF2','SSB','SUPV3L1','KHSRP','SFPQ','RBM5','HNRNPUL1','EFTUD2',
                  'EIF3D','PPIG','GRSF1','BCLAF1','FASTKD2','XRN2','RPS3','FUS','TROVE2','NKRF','GRWD1','U2AF1','HNRNPA1','DDX3X','XPO5','PRPF4',
                  'DDX52','HNRNPM','PRPF8','SLTM','SUB1','TAF15','LIN28B','NCBP2','AGGF1','CSTF2','RBFOX2','NOLC1','HNRNPU','CDC40','XRCC6','HNRNPK',
                  'RBM22','QKI','HNRNPC','SAFB','TIA1','TBRG4','FAM120A','HNRNPL','ZC3H11A','DKC1','ILF3','CSTF2T','GTF2F1','DDX59','YBX3','SUGP2',
                  'PCBP2','HLTF','PTBP1','UPF1','IGF2BP1','PABPN1','LARP7','TRA2A','MATR3','NOL12','LARP4','TIAL1','AKAP1','IGF2BP3','WDR43','SRSF7',
                  'FTO','FKBP4','UCHL5','SRSF1','FUBP3','EIF3H','SND1','SDAD1','POLR2G','SRSF9','DROSHA','DDX55','DHX30','ZNF800','PCBP1','DDX6','EXOSC5',
                  'BCCIP','LSM11','STAU2','DGCR8','FXR2','UTP18')

data <- read.csv('{}', sep='\t', header=T)
data$Group<-factor(data$Group,levels =all_k562_order )
p1=ggplot(data, aes(y=Percentage, x=Group, fill=factor(Position, levels=c0)))+geom_bar(position='stack', stat='identity') + theme(axis.text.x=element_text(angle=90)) +\
ggtitle('{}')+ theme(plot.title = element_text(face = "bold",hjust=0.5))+ guides(fill=guide_legend(title="Position"))

data <- read.csv('{}', sep='\t', header=T)
data$Group<-factor(data$Group,levels =alu_k562_order )
p2=ggplot(data, aes(y=Percentage, x=Group, fill=factor(Position, levels=c1)))+geom_bar(position='stack', stat='identity') + theme(axis.text.x=element_text(angle=90)) +\
ggtitle('{}')+ theme(plot.title = element_text(face = "bold",hjust=0.5))+ guides(fill=guide_legend(title="Position"))
pdf('{}',width=17)
ggarrange(p1,p2,heights = c(1, 1),ncol = 1, nrow = 2, align = 'v')
dev.off()

data <- read.csv('{}', sep='\t', header=T)
data$Group<-factor(data$Group,levels =all_hepg2_order )
p3=ggplot(data, aes(y=Percentage, x=Group, fill=factor(Position, levels=c0)))+geom_bar(position='stack', stat='identity') + theme(axis.text.x=element_text(angle=90)) +\
ggtitle('{}')+ theme(plot.title = element_text(face = "bold",hjust=0.5))+ guides(fill=guide_legend(title="Position"))

data <- read.csv('{}', sep='\t', header=T)
data$Group<-factor(data$Group,levels =alu_hepg2_order )
p4=ggplot(data, aes(y=Percentage, x=Group, fill=factor(Position, levels=c1)))+geom_bar(position='stack', stat='identity') + theme(axis.text.x=element_text(angle=90)) +\
ggtitle('{}')+ theme(plot.title = element_text(face = "bold",hjust=0.5))+ guides(fill=guide_legend(title="Position"))
pdf('{}',width=17)
ggarrange(p3,p4,heights = c(1, 1),ncol = 1, nrow = 2, align = 'v')
dev.off()
'''.format(TMP_DIR+'/all_peak_postion_K562.txt', 'K562_All',
           TMP_DIR + '/alu_peak_postion_K562.txt', 'K562_Alu',
           TMP_DIR + '/K562_region_dist_sort_by_order.pdf',
           TMP_DIR + '/all_peak_postion_HepG2.txt',  'HepG2_All',
           TMP_DIR + '/alu_peak_postion_HepG2.txt', 'HepG2_Alu',
           TMP_DIR + '/HepG2_region_dist_sort_by_order.pdf',))
rscript.flush()
rscript.close()

subprocess.call('Rscript '+TMP_DIR +
                '/plot_region_dist_sort_by_order.R', shell=True)

### R codes

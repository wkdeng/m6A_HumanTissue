# @author [Wankun Deng]
# @email [dengwankun@hotmail.com]
# @create date 2019-09-24 15:53:59
# @modify date 2019-11-11 15:27:32
# @desc [description]

from scipy import stats
import numpy as np
import os
import json
import math
import sys

data_path = os.path.dirname(os.path.realpath(__file__))


def count_numbers(gene_list_list, list_names, background_list=None, process_only=True, species='human',
                  id_type='UniProt', term_size=5, top_n=5, result_dir='.',width=7,height=7,text_width=50,load_data=True):
    if not load_data or not (os.path.exists(os.path.join(result_dir, 'result_array.json')) and os.path.exists(os.path.join(result_dir, 'sort_index.json'))):
        if not id_type == 'UniProt':
            id_map = {}
            for line in open(os.path.join(data_path, species + '_info.txt')):
                info = line.strip()[1:].split('|')
                if len(info) == 3:
                    id_map[info[0]] = info[2]
                    id_map[info[1]] = info[2]
            mapped_lists = []
            for i in range(len(gene_list_list)):
                mapped_lists.append([])
                for gene in gene_list_list[i]:
                    if gene in id_map.keys() and id_map[gene] not in mapped_lists[i]:
                        mapped_lists[i].append(id_map[gene])
            gene_list_list = mapped_lists
        whole_proteom_background = False
        if background_list is None:
            whole_proteom_background = True
        else:
            if not id_type == 'UniProt':
                mapped_bg = []
                for list_ in background_list:
                    mapped_bg.append([])
                    for gene in list_:
                        if gene in id_map and id_map[gene] not in mapped_bg[-1]:
                            mapped_bg[-1].append(id_map[gene])
                background_list = mapped_bg
        if not whole_proteom_background:
            all_protein = background_list
        else:
            all_protein = [[] for x in gene_list_list]
        annotated_gene = []

        go_names = {}
        for line in open(os.path.join(data_path, 'go_name.txt')):
            info = line.strip().split('\t')
            go_names[info[0]] = info[1] + '\t' + info[2]
        go_terms = []
        for i in range(len(gene_list_list)):
            go_terms.append({})
            annotated_gene.append([])

        for line in open(os.path.join(data_path, 'goa_' + species + '.gpa')):
            if not line.startswith('!'):
                info = line.strip().split('\t')
                if whole_proteom_background:
                    if info[1].strip() not in all_protein[0]:
                        [all_protein[i].append(info[1].strip())
                         for i in range(len(all_protein))]
                for i in range(len(go_terms)):
                    if info[3] not in go_terms[i]:
                        # order: m, M, n, N, eratio, pvalue
                        go_terms[i][info[3]] = [[], 0, [], 0, 0, 0]
                    if info[1] in gene_list_list[i]:
                        if info[1] not in go_terms[i][info[3]][0]:
                            go_terms[i][info[3]][0].append(info[1])
                        if info[1].strip() not in annotated_gene[i]:
                            annotated_gene[i].append(info[1].strip())
                    if info[1] not in go_terms[i][info[3]][2]:
                        go_terms[i][info[3]][2].append(info[1])
        # file_array = []
        result_array = []

        sort_index = []
        for i in range(len(gene_list_list)):
            result_array.append([])
            sort_index.append([])
            index = 0
            for id in go_terms[i]:
                go_terms[i][id][1] = len(annotated_gene[i])
                go_terms[i][id][3] = len(all_protein[i])
                go_terms[i][id][2] = len(go_terms[i][id][2])
                go_terms[i][id][0] = len(go_terms[i][id][0])
                if go_terms[i][id][0] == 0 or go_terms[i][id][1] == 0 or go_terms[i][id][2] == 0 or go_terms[i][id][
                        3] == 0:
                    continue
                go_terms[i][id][4] = (
                    go_terms[i][id][0] * go_terms[i][id][3] / float(go_terms[i][id][1] * go_terms[i][id][2]))

                if go_terms[i][id][4] > 1:

                    go_terms[i][id][5] = 1 - stats.hypergeom.cdf(go_terms[i][id][0], go_terms[i][id][3],
                                                                 go_terms[i][id][1],
                                                                 go_terms[i][id][2])
                    if go_terms[i][id][5] < 1E-22:
                        go_terms[i][id][5] = 1E-22
                    go_terms[i][id][5] = -math.log10(go_terms[i][id][5])
                    if id not in go_names:
                        go_names[id] = 'N/A\tN/A'
                    if process_only and not go_names[id].split('\t')[1] == 'biological process':
                        continue
                    if go_terms[i][id][2] < term_size:
                        continue
                    sort_index[i].append([index, go_terms[i][id][5]])
                    # print(go_terms[i][id][4])
                    index = index + 1
                    result_array[i].append(
                        [[id, go_names[id], go_terms[i][id][0]], go_terms[i][id][1], go_terms[i][id][2],
                         go_terms[i][id][3], go_terms[i][id][4], go_terms[i][id][5]])
        result_array_1 = []
        for item in result_array:
            if len(item) > 0:
                result_array_1.append(item)
        result_array = result_array_1
        sort_index_1 = []
        for item in sort_index:
            if len(item) > 0:
                sort_index_1.append(item)
        sort_index = sort_index_1
        json.dump(result_array, open(os.path.join(
            result_dir, 'result_array.json'), 'w'))
        json.dump(sort_index_1, open(os.path.join(
            result_dir, 'sort_index.json'), 'w'))
    else:
        result_array = json.load(
            open(os.path.join(result_dir, 'result_array.json')))
        sort_index = json.load(
            open(os.path.join(result_dir, 'sort_index.json')))
    to_plot = []
    y_index = []
    for i in range(len(result_array)):
        result = np.argsort(np.array(sort_index[i])[:, 1])
        for j in range(len(result)):
            term = result_array[i][result[len(result) - 1 - j]]
            term2 = term[0]
            [term2.append(x) for x in term[1:]]
            if j < top_n:
                y_index_x = '"' + str(term2[1]).split('\t')[0] + '"'
                if y_index_x not in y_index:
                    y_index.append(y_index_x)
                to_plot.append(
                    list_names[i] + '\t' + '\t'.join([str(x) for x in term2]) + '\n')
            else:
                break
    if len(to_plot) == 0:
        return 1
    r_script = open(os.path.join(result_dir, '_plot.R'), 'w')
    r_script.write("library('ggplot2')\nlibrary(stringr)\n")
    r_script.write(
        "data_table=read.delim('" + os.path.join(result_dir, "to_plot_GO.txt") + "')")
    r_script.write(
        '\ndata_table$GO_Term <- factor(data_table$GO_Term, levels = c(' + ','.join(y_index) + '))\n')
    r_script.write(
        'sp<-ggplot(data_table,aes(y=GO_Term,x=Group,color=P_value))+geom_point(aes(size=E_ratio))+'
        'scale_size_continuous(range = c(0, 10))+scale_y_discrete(labels = function(x) str_wrap(x, width = %s))\n'%text_width)
    r_script.write(
        "plotx<-sp+scale_color_gradientn(colours = rev(rainbow(5)))+ labs(color='-log10(p Value)')+theme(panel.grid.major= element_blank(),"
        "panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 60, hjust = 1))\n")
    r_script.write(
        "pdf('" + os.path.join(result_dir, "plot_result_GO.pdf") + "',width=%s,height=%s)\n"%(width,height))
    r_script.write('plot(plotx)\ndev.off()\n')
    r_script.flush()
    fp = open(os.path.join(result_dir, 'to_plot_GO.txt'), 'w')
    fp.write('Group	ID	GO_Term	Term	m	M	n	N	E_ratio	P_value\n')
    for i in range(len(to_plot)):
        fp.write(to_plot[i])
    fp.flush()
    fp.close()
    os.system('Rscript ' + os.path.join(result_dir, '_plot.R'))
    return 0


def go_analysis(args):
    gene_list_list = []
    list_names = []

    file_list = []
    species = 'human'
    top_n = 5
    process_only = True
    id_type = 'ENSG'
    term_size = 5
    result_dir = '.'

    for arg in args:
        if '=' in arg:
            info = arg.strip().split('=')
            info = [x.strip() for x in info]
            if info[0] == 'species':
                species = info[1]
            elif info[0] == 'process_only':
                process_only = info[1] == 'True'
            elif info[0] == 'term_size':
                term_size = int(info[1])
            elif info[0] == 'id_type':
                id_type = info[1]
            elif info[0] == 'top_n':
                top_n = int(info[1])
            elif info[0] == 'result_dir':
                result_dir = info[1]
        else:
            file_list.append(arg)

    for file in file_list:
        gene_list = []
        list_names.append(os.path.split(file)[1].replace(
            '.txt', '').replace('up_mGene-', ''))
        for line in open(file):
            gene_list.append(line.strip())
        gene_list_list.append(gene_list)

    count_numbers(gene_list_list, list_names, process_only=process_only, species=species, top_n=top_n,
                  term_size=term_size, id_type=id_type, result_dir=result_dir)


if __name__ == "__main__":
    # go_analysis(
    #     [
    #         r'H:/yi_lab/m6a/data/samir/up_mGene-AH-Uniq.txt', r'H:/yi_lab/m6a/data/samir/up_mGene-NH-Uniq.txt',
    #         r'H:/yi_lab/m6a/data/samir/up_mGene-NH+AH-Common.txt',
    #         r'H:/yi_lab/m6a/data/samir/up_mGene-NH+AH-Background.txt',
    #         r'H:/yi_lab/m6a/data/samir/up_mGene-AK-Uniq.txt',
    #         r'H:/yi_lab/m6a/data/samir/up_mGene-NK-Uniq.txt',
    #         r'H:/yi_lab/m6a/data/samir/up_mGene-AK+NK-Common.txt',
    #         r'H:/yi_lab/m6a/data/samir/up_mGene-NK+AK-Background.txt',
    #         'top_n=5', "term_size=10", 'species=mouse',
    #         r'result_dir=H:/yi_lab/m6a/src/scripts/GO_enrichment/batch3_Neo_Heart_Kidney_', 'id_type=UniProt'])
    go_analysis(sys.argv[1:])

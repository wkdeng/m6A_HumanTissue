# @author [Wankun Deng]
# @email [dengwankun@hotmail.com]
# @create date 2017-09-18 16:33:00
# @modify date 2019-12-11 16:33:00
# @desc [description]

import os
import re

###**--- IN-FILE CONFIG ---**###
SNAKEMAKE_FILE_DIR = config['SCRIPTS_DIR']
PROJECT_DIR = config['PROJECT_DIR']
PAIRED_END = False if 'paired_end' not in config else config['paired_end']
INCLUDE_MREAD_ANALYSIS = True if 'include_mread_an  alysis' not in config else config[
    'include_mread_analysis']

GENOME = config['genome']
GTF_PATH = config[GENOME]['gtf']
CHR_SIZE = config[GENOME]['chr_size']

SAMPLE_TYPE_DICT = config['sample_type_dict']
COMPARISON_LIST = config['sample_comparison']
MAX_TAGS = config['clam']['max_tags']
CONFIG_FILE = config['config_file']
DEG_GROUP=config['deg_group'] if 'deg_group' in config else ''
FRACTIONS=config['fractions'] if 'factions' in config else ''
IS_TRIPSEQ=config['is_tripseq']==True if 'is_tripseq' in config else False
USE_MULTIPLE_PEAK=config['use_multiple_peak']==True if 'use_multiple_peak' in config else False
PYTHON2=config['python2_path']
PYTHON3=config['default_python']
SAMPLE_MERGE=config['sample_merge']
GRP_MERGE=config['grp_merge']
REP_MERGE=config['rep_merge']

###**--- SNAKEMAKE FILE ---**##

rule all:
    input:
        "{project_dir}/summary/merged_summary.pdf".format(
            project_dir=PROJECT_DIR)

rule summary:
    input:
        # # require mapping stats
        # mapping_stats=["{project_dir}/star/{sample}/mapping_stats.txt".format(
        #     project_dir=PROJECT_DIR,
        #     sample=x)
        #     for x in SAMPLE_TYPE_DICT
        #  ],
        exp=["{project_dir}/kallisto/{sample_name}/abundance.tsv".format(
            project_dir=PROJECT_DIR,
            sample_name=x)
            for x in SAMPLE_TYPE_DICT
         ],
        # require combined homer
        homer_combined=["{project_dir}/homer/{ip_sample}__{con_sample}/clam_combined/homerResults.html".format(
            project_dir=PROJECT_DIR,
            ip_sample=x[0],
            con_sample=x[1])
            for x in COMPARISON_LIST
         ],
        # require combined topology clam
        dist_combined=[ "{project_dir}/sanity_check/{sample_merge}_combined_dist.pdf".format(
            project_dir=PROJECT_DIR,
            sample_merge=x)
            for x in SAMPLE_MERGE
         ],
        # require venn diagram
        rep_compare=["{project_dir}/sanity_check/{sample_merge}_combined_rep_compare.tiff".format(
            project_dir=PROJECT_DIR,
            sample_merge=x)
            for x in SAMPLE_MERGE
         ],
        #  peak_combined=["{project_dir}/clam/peaks/peaks-{ip_sample}__{con_sample}/narrow_peak.combined.bed".format(
        #     project_dir=PROJECT_DIR,
        #     ip_sample=x[0],
        #     con_sample=x[1])
        #     for x in COMPARISON_LIST 
        #  ],
        # # require peak instensities
        parsed_peaks="{project_dir}/parsed_peaks/peak_intensity.csv".format(
            project_dir=PROJECT_DIR),
        
        # require unmerged bw files
        # unmerged_bw=["{project_dir}/bigwig/{sample_name}/{sample_name}_multi.bw".format(
        #     project_dir=PROJECT_DIR,
        #     sample_name=x)
        #     for x in SAMPLE_TYPE_DICT
        #  ],
        # # require bw files
        # bw_files=["{project_dir}/bigwig/{sample}_{group}_merged.bw".format(
        #     project_dir=PROJECT_DIR,
        #     sample=x,
        #     group=y)
        #     for x in SAMPLE_MERGE for y in GRP_MERGE
        #  ],
         # require normalized bw files
         bw_n_files=["{project_dir}/bigwig/{sample}_normalized.bw".format(
            project_dir=PROJECT_DIR,
            sample=x)
            for x in SAMPLE_MERGE
         ],
        # # require normalized bw files for both replicates
        # bw_n_files_both=["{project_dir}/bigwig/{sample}_{rep}_normalized.bw".format(
        #     project_dir=PROJECT_DIR,
        #     sample=x,rep=y)
        #     for x in SAMPLE_MERGE for y in REP_MERGE
        #  ],
    output:
        "{project_dir}/summary/merged_summary.pdf".format(
            project_dir=PROJECT_DIR)
    params:
        out_html = "{project_dir}/reports/report_summary.html".format(
            project_dir=PROJECT_DIR),
        script_dir = SNAKEMAKE_FILE_DIR+'/scripts',
        config_file = CONFIG_FILE,
        python_path=PYTHON3
    shell:
        """
        {params.python_path} {params.script_dir}/summary/summary_m6a.py {params.config_file}
        """
rule unzip_reads:
    input:
        "{project_dir}/reads/{sample_name}.fastq.gz"
    output:
        temporary("{project_dir}/reads/{sample_name}.fastq")
    shell:
        """
        gunzip -c {input} > {output}
        """

### alignment and preprocessing
rule star_map:
    input:
        "{project_dir}/reads/{sample_name}.fastq"
    output:
        "{project_dir}/star/{sample_name}/Aligned.sortedByCoord.out.bam"
    log:
        "{project_dir}/logs/star/{sample_name}.log"
    params:
        reads = "{project_dir}/reads/{sample_name}.fastq",
        prefix = "{project_dir}/star/{sample_name}/",
        index = config[GENOME]['star_idx'],
        max_hits = 100
    threads: 20
    shell:
        """
        STAR --genomeDir {params.index} \
        --readFilesIn {params.reads}  --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {params.prefix} \
        --outFilterMultimapNmax {params.max_hits} \
        --runThreadN {threads} \
        --twopassMode Basic \
        --outStd Log >{log} 2>&1
        """

rule fastqc:
    input:
        "{project_dir}/reads/{sample_name}.fastq"
    output:
        "{project_dir}/fastqc/{sample_name}_fastqc.zip"
    log:
        "{project_dir}/fastqc/{sample_name}.log"
    params:
        out_dir = "{project_dir}/fastqc"
    shell:
        """
        fastqc -o {params.out_dir} {input} > {log} 2>&1
        """

rule mapping_stat:
    input:
        align = "{project_dir}/star/{sample_name}/Aligned.sortedByCoord.out.bam",
        qc_zip="{project_dir}/fastqc/{sample_name}_fastqc.zip"
    output:
        "{project_dir}/star/{sample_name}/mapping_stats.txt"
    params:
        scripts = SNAKEMAKE_FILE_DIR + '/scripts/summary/mapping_stat.py',
        python_path=PYTHON2
    shell:
        "{params.python_path} {params.scripts} {input.align} {input.qc_zip}> {output}"

### get gene expression
rule kallisto_quant:
    input:
        "{project_dir}/reads/{sample_name}.fastq"
    output:
        "{project_dir}/kallisto/{sample_name}/abundance.tsv"
    log:
        "{project_dir}/logs/kallisto/{sample_name}.log"
    params:
        outdir = "{project_dir}/kallisto/{sample_name}/",
        index = config[GENOME]['kallisto_idx'],
        out_log = "{project_dir}/kallisto/foo.txt"
    threads: 20
    shell:
        """
        kallisto quant -t {threads} -i {params.index} -o {params.outdir} --single -l 51 -s 5 {input}
        echo '{input} complete' >> {params.out_log}
        """


### peak prepare
rule clam_prep:
    input:
        align = "{project_dir}/star/{sample_name}/Aligned.sortedByCoord.out.bam"
    output:
        unique = "{project_dir}/clam/{sample_name}/unique.sorted.collapsed.bam" if MAX_TAGS > 0 else \
            "{project_dir}/clam/{sample_name}/unique.sorted.bam",
        multi = "{project_dir}/clam/{sample_name}/multi.sorted.collapsed.bam" if MAX_TAGS > 0 else \
            "{project_dir}/clam/{sample_name}/multi.sorted.bam"
    log:
        "{project_dir}/logs/clam/{sample_name}-prep.log"
    params:
        outdir = "{project_dir}/clam/{sample_name}",
        tagger_method = "median",
        max_tags = MAX_TAGS,
        preprocessor_script=SNAKEMAKE_FILE_DIR+'/mouse_lowinput/scripts/preprocessor1.2.py',
        python_path=PYTHON2
    shell:
        "{params.python_path} {params.preprocessor_script} {input.align} {params.outdir}  {params.tagger_method} 100 {params.max_tags} >{log} 2>&1"

rule clam_em:
    input:
        unique = "{project_dir}/clam/{sample_name}/unique.sorted.collapsed.bam" if MAX_TAGS > 0 else \
            "{project_dir}/clam/{sample_name}/unique.sorted.bam",
        multi = "{project_dir}/clam/{sample_name}/multi.sorted.collapsed.bam" if MAX_TAGS > 0 else \
            "{project_dir}/clam/{sample_name}/multi.sorted.bam"
    output:
        "{project_dir}/clam/{sample_name}/realigned.sorted.bam"
    log:
        "{project_dir}/logs/clam/{sample_name}-em.log"
    params:
        outdir = "{project_dir}/clam/{sample_name}",
        max_tags = MAX_TAGS,
        winsize = 100
    shell:
        """CLAM realigner -i {input.unique} -o {params.outdir} --winsize {params.winsize} \
        --max-tags {params.max_tags} --lib-type unstranded >{log} 2>&1"""

### peak calling

rule clam_callpeak:
    input:
        ip_prep = "{project_dir}/clam/{ip_sample}/unique.sorted.collapsed.bam" if MAX_TAGS > 0 else \
            "{project_dir}/clam/{ip_sample}/unique.sorted.bam",
        con_prep = "{project_dir}/clam/{con_sample}/unique.sorted.collapsed.bam" if MAX_TAGS > 0 else \
            "{project_dir}/clam/{con_sample}/unique.sorted.bam"
    output:
        "{project_dir}/clam/peaks/peaks-{ip_sample}__{con_sample}/narrow_peak.unique.bed"
    log:
        "{project_dir}/logs/clam/{ip_sample}__{con_sample}-callpeak.log"
    params:
        outdir = "{project_dir}/clam/peaks/peaks-{ip_sample}__{con_sample}",
        gtf = config[GENOME]['gtf'],
        binsize = 100,
        qval_cutoff = 0.5,
        fold_change = '0.1',
        threads = 20
    shell:
        """
        CLAM peakcaller -i {input.ip_prep}  -c {input.con_prep} \
            -p {params.threads} --lib-type unstranded \
            -o {params.outdir} --gtf {params.gtf} --unique-only --binsize {params.binsize} \
            --qval-cutoff {params.qval_cutoff} --fold-change {params.fold_change} >{log} 2>&1
        mv {output} {output}.all
        awk '$9<0.005 && $7>1.609' {output}.all > {output}
        """

rule clam_callpeak_mread:
    input:
        ip_prep = [
            "{project_dir}/clam/{ip_sample}/unique.sorted.collapsed.bam" if MAX_TAGS > 0 else "{project_dir}/clam/{ip_sample}/unique.sorted.bam",
            "{project_dir}/clam/{ip_sample}/realigned.sorted.bam"],
        con_prep = [
            "{project_dir}/clam/{con_sample}/unique.sorted.collapsed.bam" if MAX_TAGS > 0 else "{project_dir}/clam/{con_sample}/unique.sorted.bam",
            "{project_dir}/clam/{con_sample}/realigned.sorted.bam"]
    output:
        "{project_dir}/clam/peaks/peaks-{ip_sample}__{con_sample}/narrow_peak.combined.bed"
    log:
        "{project_dir}/logs/clam/{ip_sample}__{con_sample}-callpeak_mread.log"
    params:
        outdir = "{project_dir}/clam/peaks/peaks-{ip_sample}__{con_sample}",
        gtf = config[GENOME]['gtf'],
        binsize = 100,
        qval_cutoff = 0.5,
        fold_change = '0.1',
        threads = 4
    shell:
        """
        CLAM peakcaller -i {input.ip_prep}  -c {input.con_prep} \
            -p {params.threads} --lib-type unstranded \
            -o {params.outdir} --gtf {params.gtf} --binsize {params.binsize} \
            --qval-cutoff {params.qval_cutoff} --fold-change {params.fold_change} >{log} 2>&1
        mv {output} {output}.all
        awk '$9<0.05' {output}.all > {output}
        """

### evaluations
rule homer_motif_m:
    input:
        peak_fn = "{project_dir}/clam/peaks/peaks-{ip_sample}__{con_sample}/narrow_peak.combined.bed"
    output:
        "{project_dir}/homer/{ip_sample}__{con_sample}/clam_combined/homerResults.html"
    params:
        outdir = "{project_dir}/homer/{ip_sample}__{con_sample}/clam_combined",
        motif_len = '5,6,7',
        genome = GENOME,
        nthread = 4,
        size = 100,
        motif_num = 10
    log:
        "{project_dir}/logs/homer/log.homer.{ip_sample}__{con_sample}.txt"
    shell:
        "findMotifsGenome.pl {input.peak_fn} {params.genome} {params.outdir} "\
            " -rna -len {params.motif_len} "\
            "-p {params.nthread} -size {params.size} -S {params.motif_num} >{log} 2>&1"

rule topology_dist_mreads:
    input:
        peak_fn =["{project_dir}/clam/peaks/peaks-{sample_merge}-%s-IP__{sample_merge}-%s-Inp/narrow_peak.combined.bed"%(x,x) for x in REP_MERGE]
    output:
        dist_img = "{project_dir}/sanity_check/{sample_merge}_combined_dist.pdf"
    log:
        "{project_dir}/logs/sanity_check/{sample_merge}.log"
    params:
        outdir = "{project_dir}/sanity_check",
        count_script = SNAKEMAKE_FILE_DIR + \
            "/scripts/topological_dist/Run_dist_on_utr_multi.py",
        dist_script= SNAKEMAKE_FILE_DIR + \
            "/scripts/topological_dist/Peak_distribution_on_utr_cds.py",
        count_out="{project_dir}/sanity_check/{sample_merge}_combined.dist",
        plot_script = SNAKEMAKE_FILE_DIR + "/scripts/topological_dist/Plot_multi_dist.R",
        args=' '.join(["{sample_merge}_%s {project_dir}/clam/peaks/peaks-{sample_merge}-%s-IP__{sample_merge}-%s-Inp/narrow_peak.combined.bed"%(x,x,x) for x in REP_MERGE]),
        genome = GENOME,
        binnum = 50,
        python_path=PYTHON3,
        reps=" ".join(REP_MERGE)
    shell:
        """
        {params.python_path} {params.count_script} {params.dist_script} {params.count_out} {params.args} >{log} 2>&1
        Rscript {params.plot_script} {params.count_out} {params.binnum} {output}
        """

rule rep_compare:
    input:
        peak_fn =["{project_dir}/clam/peaks/peaks-{sample_merge}-%s-IP__{sample_merge}-%s-Inp/narrow_peak.combined.bed"%(x,x) for x in REP_MERGE]
    output:
        "{project_dir}/sanity_check/{sample_merge}_combined_rep_compare.tiff"
    log:
        "{project_dir}/logs/sanity_check/{sample_merge}.log"
    params:
        project_dir="{project_dir}",
        sample_merge="{sample_merge}",
        outdir = "{project_dir}/sanity_check",
        plot_script = SNAKEMAKE_FILE_DIR + "/scripts/summary/plot_rep_peaks_venn.R",
        reps=' '.join(REP_MERGE),
        peaks=" ".join(["{project_dir}/sanity_check/{sample_merge}_%s.bed"%x for x in REP_MERGE])
    shell:
        """
        for rep in {params.reps}
        do
            cat {params.project_dir}/clam/peaks/peaks-{params.sample_merge}-$rep-IP__{params.sample_merge}-$rep-Inp/narrow_peak.combined.bed | \
                awk '{{printf("%s:%s:%s:%s\\n",$1,$2,$3,$6)}}' > {params.outdir}/{params.sample_merge}_$rep.bed
        done
        Rscript {params.plot_script} {output} {params.peaks}  > {log} 2>&1
        """


rule parse_peaks:
    input:
        ["{project_dir}/kallisto/{sample_name}/abundance.tsv".format(
            project_dir=PROJECT_DIR,
            sample_name=x)
            for x in SAMPLE_TYPE_DICT if SAMPLE_TYPE_DICT[x] == 'con'],
        ["{project_dir}/clam/peaks/peaks-{ip_sample}__{con_sample}/narrow_peak.{unique}.bed".format(
            project_dir=PROJECT_DIR,
            ip_sample=x[0],
            con_sample=x[1],unique='combined' if USE_MULTIPLE_PEAK else 'unique')
            for x in COMPARISON_LIST]
    output:
        "{project_dir}/parsed_peaks/input_peak.RPKM.csv",
        "{project_dir}/parsed_peaks/ip_peak.RPKM.csv",
        "{project_dir}/parsed_peaks/peak_intensity.csv"
    log:  
        "{project_dir}/logs/diff_peaks/parse_peaks.log"
    params:
        config_file = CONFIG_FILE,
        script = SNAKEMAKE_FILE_DIR+"/scripts/diff_peaks/parse_peaks.py",
        python_path=PYTHON3
    shell:
        "{params.python_path} {params.script} {params.config_file} > {log} 2>&1"

rule pile_reads_for_bw:
    input:
        unique = "{project_dir}/clam/{sample_name}/unique.sorted.collapsed.bam" if MAX_TAGS > 0 else \
            "{project_dir}/clam/{sample_name}/unique.sorted.bam",
        multi = "{project_dir}/clam/{sample_name}/realigned.sorted.bam"
    output:
        "{project_dir}/bigwig/{sample_name}/{sample_name}_multi.bw"
    log:
        "{project_dir}/logs/clam/{sample_name}-bw.log"
    params:
        sample_name="{sample_name}",
        python_path=PYTHON3,
        output_dir="{project_dir}/bigwig/{sample_name}",
        script = SNAKEMAKE_FILE_DIR+"/scripts/make_bw/pileup_realigned_reads.py",
        max_tags = MAX_TAGS,
        gtf = config[GENOME]['gtf'],
        chr_size = config[GENOME]['chr_size'],
        winsize = 100
    shell:
        """{params.python_path} {params.script} -s {params.sample_name} \
            --ub {input.unique} --rb {input.multi} \
            -g {params.gtf} -c {params.chr_size} -u -o {params.output_dir} \
            -m median -t 20 > {log} 2>&1"""

rule merge_replicate_bw:
    input:
        rep1="{project_dir}/bigwig/{sample}-Rep1-{group}/{sample}-Rep1-{group}_multi.bw",
        rep2="{project_dir}/bigwig/{sample}-Rep2-{group}/{sample}-Rep2-{group}_multi.bw"
    output:
        "{project_dir}/bigwig/{sample}_{group}_merged.bw"
    log:
        "{project_dir}/logs/clam/{sample}_{group}_merged.log"
    params:
        sample="{sample}",
        group="{group}",
        python_path=PYTHON3,
        chr_size = config[GENOME]['chr_size'],
        output_dir="{project_dir}/bigwig"
    shell:
        """
        bigWigMerge {input.rep1} {input.rep2} {params.output_dir}/{params.sample}_{params.group}_merged.bedGraph

        cat {params.output_dir}/{params.sample}_{params.group}_merged.bedGraph \
            | awk '{{printf("%s\\t%s\\t%s\\t%s\\n",$1,$2,$3,$4/2)}}' \
            > {params.output_dir}/{params.sample}_{params.group}_merged_mean.bedGraph

        LC_COLLATE=C sort -k1,1 -k2,2n {params.output_dir}/{params.sample}_{params.group}_merged_mean.bedGraph \
            > {params.output_dir}/{params.sample}_{params.group}_merged_sorted.bedGraph
            
        bedGraphToBigWig {params.output_dir}/{params.sample}_{params.group}_merged_sorted.bedGraph \
            {params.chr_size} {params.output_dir}/{params.sample}_{params.group}_merged.bw
        """
rule normalize_bw:
    input:
        ip="{project_dir}/bigwig/{sample}_IP_merged.bw",
        inp="{project_dir}/bigwig/{sample}_Inp_merged.bw"
    output:
        "{project_dir}/bigwig/{sample}_normalized.bw"
    log:
        "{project_dir}/logs/clam/{sample}_normalized.log"
    params:
        sample="{sample}",
        python_path=PYTHON3,
        script = SNAKEMAKE_FILE_DIR+"/scripts/make_bw/normalize_ip_by_inp.py",
        gtf = GTF_PATH
    shell:
        """
        {params.python_path} {params.script} {input.ip} {input.inp} {params.gtf} {output}
        """

rule normalize_bw_rep:
    input:
        ip="{project_dir}/bigwig/{sample}-{rep}-IP/{sample}-{rep}-IP_multi.bw",
        inp="{project_dir}/bigwig/{sample}-{rep}-Inp/{sample}-{rep}-Inp_multi.bw"
    output:
        "{project_dir}/bigwig/{sample}_{rep}_normalized.bw"
    log:
        "{project_dir}/logs/clam/{sample}_{rep}_normalized.log"
    params:
        sample="{sample}",
        python_path=PYTHON3,
        script = SNAKEMAKE_FILE_DIR+"/scripts/make_bw/normalize_ip_by_inp.py",
        gtf = GTF_PATH
    shell:
        """
        {params.python_path} {params.script} {input.ip} {input.inp} {params.gtf} {output}
        """

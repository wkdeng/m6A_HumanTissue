###
## @author [Wankun Deng]
## @email [dengwankun@gmail.com]
## @create date 2023-03-06 01:14:05
## @modify date 2023-03-06 01:14:05
## @desc [description]
###
SNAKEMAKE_FILE_DIR = config['SCRIPTS_DIR']
PROJECT_DIR = config['PROJECT_DIR']
GENOME = config['genome']
SAMPLES = config['sample_type_dict']
SAMPLE_COMPARISON=config['sample_comparison']
T2G = config[GENOME]['t2g']
GENE_NAME = config[GENOME]['gene_name']
GTF_PATH = config[GENOME]['gtf']
READ_LENGTH = config['read_length']
rRNA_MASK = config[GENOME]['rmask']
CHR_SIZE = config[GENOME]['chr_size']
CONFIG = config['CONFIG']
MAX_TAGS=-1
PYTHON3=config['default_python']
SAMPLE_MERGE=config['sample_merge']
REP_MERGE=config['rep_merge']

rule all:
    input:
        fastqc=[f"{PROJECT_DIR}/fastqc/{sample_name}_R1.adapterTrim.round2_fastqc.zip" for sample_name in SAMPLES],
        exp=[f"{PROJECT_DIR}/kallisto/{sample_name}/abundance.tsv" for sample_name in SAMPLES if SAMPLES[sample_name] == 'con'],
        mapping_stat=[f"{PROJECT_DIR}/star/{sample_name}/mapping_stats.txt" for sample_name in SAMPLES],
        uniq_peak=[f"{PROJECT_DIR}/clam/peaks/peaks-{ip_sample}__{con_sample}/narrow_peak.unique.bed" for ip_sample, con_sample in SAMPLE_COMPARISON],
        combined_peak=[f"{PROJECT_DIR}/clam/peaks/peaks-{ip_sample}__{con_sample}/narrow_peak.combined.bed" for ip_sample, con_sample in SAMPLE_COMPARISON],
        motif_m=[f"{PROJECT_DIR}/homer/{ip_sample}__{con_sample}/clam_combined/homerResults.html" for ip_sample, con_sample in SAMPLE_COMPARISON],
        dist_combined=["{project_dir}/sanity_check/{sample_merge}_combined_dist.pdf".format(project_dir=PROJECT_DIR,sample_merge=x) for x in SAMPLE_MERGE],
        rep_compare=["{project_dir}/sanity_check/{sample_merge}_combined_rep_compare.tiff".format(project_dir=PROJECT_DIR,sample_merge=x) for x in SAMPLE_MERGE],


rule unzip_fq:
    input:
        first = "{project_dir}/reads/{sample_name}_R1.fastq.gz",
        second = "{project_dir}/reads/{sample_name}_R2.fastq.gz"
    output:
        first = temp("{project_dir}/reads/{sample_name}_R1.fastq"),
        second = temp("{project_dir}/reads/{sample_name}_R2.fastq")
    params:
        name = "{sample_name}",
        output_dir = "{project_dir}/fq_cutadapt"
    log:
        "{project_dir}/logs/archive_fq/{sample_name}.log"
    shell:
        """
gzip -c -d {input.first} > {output.first}
gzip -c -d {input.second} > {output.second}
echo `date` > {log}
        """


rule cutadapt:
    input:
        first = "{project_dir}/reads/{sample_name}_R1.fastq",
        second = "{project_dir}/reads/{sample_name}_R2.fastq"
    output:
        first = temp("{project_dir}/fq_cutadapt/{sample_name}_R1.adapterTrim.round2.fastq"),
        second = temp("{project_dir}/fq_cutadapt/{sample_name}_R2.adapterTrim.round2.fastq")
    params:
        name = "{sample_name}",
        output_dir = "{project_dir}/fq_cutadapt"
    log:
        "{project_dir}/logs/cutadapt/{sample_name}.log"
    shell:
        """
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 \
-a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT \
-A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA \
-A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA \
-A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC \
-A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT \
-o {params.output_dir}/{params.name}_R1.adapterTrim.fastq -p {params.output_dir}/{params.name}_R2.adapterTrim.fastq \
{input.first} {input.second} > {log} 2>&1
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 \
-A AACTTGTAGATCGGA -A AGGACCAAGATCGGA \
-A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG \
-A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG \
-A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG \
-A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC \
-A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT \
-o {params.output_dir}/{params.name}_R1.adapterTrim.round2.fastq -p {params.output_dir}/{params.name}_R2.adapterTrim.round2.fastq \
{params.output_dir}/{params.name}_R1.adapterTrim.fastq {params.output_dir}/{params.name}_R2.adapterTrim.fastq > {log} 2>&1
rm {params.output_dir}/{params.name}_R1.adapterTrim.fastq
rm {params.output_dir}/{params.name}_R2.adapterTrim.fastq
        """

rule fastqc:
    input:
        first = "{project_dir}/fq_cutadapt/{sample_name}_R1.adapterTrim.round2.fastq",
        second = "{project_dir}/fq_cutadapt/{sample_name}_R2.adapterTrim.round2.fastq"
    output:
        first = "{project_dir}/fastqc/{sample_name}_R1.adapterTrim.round2_fastqc.zip",
        second = "{project_dir}/fastqc/{sample_name}_R2.adapterTrim.round2_fastqc.zip"
    params:
        out_dir = "{project_dir}/fastqc"
    log:
        "{project_dir}/logs/fastqc/{sample_name}.log"
    shell:
        """
        fastqc -o {params.out_dir} {input.first} > {log} 2>&1
        fastqc -o {params.out_dir} {input.second} >> {log} 2>&1
        """

rule kallisto_quant:
    input:
        first = "{project_dir}/fq_cutadapt/{sample_name}_R1.adapterTrim.round2.fastq",
        second = "{project_dir}/fq_cutadapt/{sample_name}_R2.adapterTrim.round2.fastq"
    output:
        "{project_dir}/kallisto/{sample_name}/abundance.tsv"
    log:
        "{project_dir}/logs/kallisto/{sample_name}.log"
    params:
        reads = "{project_dir}/fq_cutadapt/{sample_name}",
        outdir = "{project_dir}/kallisto/{sample_name}",
        index = config[GENOME]['kallisto_idx'],
        out_log = "{project_dir}/kallisto/foo.txt"
    threads: 5
    shell:
        """
        kallisto quant -t {threads} -i {params.index} -o {params.outdir}/ {params.reads}_R1.adapterTrim.round2.fastq {params.reads}_R2.adapterTrim.round2.fastq > {log} 2>&1
        echo '{input} complete' >> {params.out_log}
        """

rule star_map:
    input:
        first = "{project_dir}/fq_cutadapt/{sample_name}_R1.adapterTrim.round2.fastq",
        second = "{project_dir}/fq_cutadapt/{sample_name}_R2.adapterTrim.round2.fastq"
    output:
        "{project_dir}/star/{sample_name}/Aligned.out.bam"
    log:
        "{project_dir}/logs/star/{sample_name}.log"
    params:
        reads = "{project_dir}/fq_cutadapt/{sample_name}",
        prefix = "{project_dir}/star/{sample_name}",
        index = config[GENOME]['star_idx'],
        max_hits = 100
    threads: 10
    shell:
        """
        STAR --genomeDir {params.index} \
        --readFilesIn {params.reads}_R1.adapterTrim.round2.fastq {params.reads}_R2.adapterTrim.round2.fastq  --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params.prefix}/ \
        --outFilterMultimapNmax {params.max_hits} \
        --runThreadN {threads} \
        --alignEndsProtrude 15 ConcordantPair \
        --twopassMode Basic --limitOutSJcollapsed 2000000 \
        --outStd Log > {log} 2>&1
        """

rule mapping_stat:
    input:
        align = "{project_dir}/star/{sample_name}/Aligned.out.bam",
        script = SNAKEMAKE_FILE_DIR + '/scripts/summary/mapping_stat.py',
        qc_zip='{project_dir}/fastqc/{sample_name}_R1.adapterTrim.round2_fastqc.zip'
    output:
        "{project_dir}/star/{sample_name}/mapping_stats.txt"
    log:
        "{project_dir}/logs/star/{sample_name}.mapping_stats.log"
    shell:
        "python {input.script} {input.align} {input.qc_zip}> {output}"

rule mask_rRNA:
    input:
        "{project_dir}/star/{sample_name}/Aligned.out.bam"
    output:
        "{project_dir}/star/{sample_name}/Aligned.out.rRNA_mask.bam"
    params:
        rRNA_annotation = rRNA_MASK
    log:
        "{project_dir}/logs/star/{sample_name}.mask_rna.log"
    shell:
        '''
        bedtools intersect -f 0.90 -abam {input} -b {params.rRNA_annotation} -v > {output}
        '''

rule collapse_dup:
    input:
        "{project_dir}/star/{sample_name}/Aligned.out.rRNA_mask.bam"
    output:
        "{project_dir}/star/{sample_name}/Aligned.out.rRNA_mask.dup_removed.bam"
    params:
        metric = "{project_dir}/star/{sample_name}/dup_removal.metrics.txt",
        script = SNAKEMAKE_FILE_DIR+'/scripts/prepare_reads/collapse_pcr.py'
    log:
        "{project_dir}/logs/star/{sample_name}.clps_dup.log"
    shell:
        '''
        python {params.script} -b {input} -o {output} -m {params.metric} > {log} 2>&1
        '''

rule clam_preprocess:
    input:
        "{project_dir}/star/{sample_name}/Aligned.out.rRNA_mask.dup_removed.bam"
    output:
        "{project_dir}/clam/{sample_name}/multi.sorted.bam", "{project_dir}/clam/{sample_name}/unique.sorted.bam"
    params:
        out_dir = "{project_dir}/clam/{sample_name}"
    log:
        "{project_dir}/logs/clam/{sample_name}.preprocess.log"
    shell:
        '''
        CLAM preprocessor -i {input} -o {params.out_dir} --read-tagger-method start --lib-type sense> {log} 2>&1
        '''

rule clam_realign:
    input:
        "{project_dir}/clam/{sample_name}/multi.sorted.bam", "{project_dir}/clam/{sample_name}/unique.sorted.bam"
    output:
        "{project_dir}/clam/{sample_name}/realigned.sorted.bam"
    params:
        input_dir = "{project_dir}/clam/{sample_name}"
    log:
        "{project_dir}/logs/clam/{sample_name}.realign.log"
    shell:
        '''
        CLAM realigner -i {params.input_dir}/multi.sorted.bam -o {params.input_dir} --winsize 50 --max-tags -1 --read-tagger-method start --lib-type sense > {log} 2>&1
        '''

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


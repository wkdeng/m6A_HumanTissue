###
## @author [Wankun Deng]
## @email [dengwankun@hotmail.com]
## @create date 2022-03-15 17:30:50
## @modify date 2022-03-15 17:30:50
## @desc [description]
###
from scripts.common import *
def read_gtf(fn):
    gene_annot = {}
    with open(fn, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            ele = line.strip().split('\t')
            if ele[2].lower() != 'gene':
                continue
            chr_, start, end, strand = ele[0], int(ele[3]), int(ele[4]), ele[6]
            try:
                if chr_.startswith('ERCC'):
                    gene_id=chr_
                else:
                    gene_id = re.findall(r"(\w+)", ele[-1])[1]
            except AttributeError:
                continue
            gene_annot[gene_id] = [chr_, start, end, strand, gene_id]
    return gene_annot

def scan_motif(intervals):
    sequences=intervals.sequence(fi=fasta,s=True)
    current_peak=''
    breakdown=[]
    motifs=[]
    count=0
    for line in open(sequences.seqfn):
        if line.startswith('>'):
            current_peak=line[1:].strip()
        else:
            seq=line.strip()
            peak_info=re.split(':|-|\(',current_peak)
            peak_info[1]=int(peak_info[1])
            peak_info[2]=int(peak_info[2])
            peak_info[3]=current_peak[-2]
            peak_name='%s:%s:%s:%s'%(peak_info[0],peak_info[1],peak_info[2],peak_info[3])
            peak_len=peak_info[2]-peak_info[1]
            flag=False
            for match in re.finditer('[AG][AG]AC[ACT]',seq):#[AG][AG]AC[ACT]
                if not flag:
                    count+=1
                    flag=True
                span=list(match.span())
                motif_seq=seq[span[0]:span[1]]
                residues=['R1','R2','A3','C4','H5']
                if peak_info[3]=='-':
                    residues=['H5','C4','A3','R2','R1']
                    span[0]=peak_len-span[1]
                motif_name=peak_info[0]+':%s:%s'%(peak_info[1]+span[0],peak_info[3])
                motifs.append([peak_info[0],peak_info[1]+span[0],peak_info[1]+span[0]+5,match[0],1000,peak_info[3]])
                for i in range(5):
                    breakdown.append([peak_info[0],
                                      peak_info[1]+span[0]+i,
                                      peak_info[1]+span[0]+i+1,
                                      motif_name,
                                      peak_name,
                                      residues[i],
                                      peak_info[3]])
    breakdown=pd.DataFrame(breakdown)
    breakdown[[1]]=breakdown[[1]].astype(str)
    breakdown[[2]]=breakdown[[2]].astype(str)
    breakdown.columns=['chr','start','end','motif','peak_name','residue','strand']
    breakdown['motif']=breakdown['motif']+breakdown['residue']
    breakdown=breakdown.drop_duplicates(subset=['chr','start','end'])
    breakdown=breakdown[(breakdown['chr']!='chrY')&(breakdown['chr']!='chrX')&(breakdown['chr']!='chrM')]
    
    return motifs,breakdown

def get_clean_pas(gene_list=None):
    exp=pd.read_csv('m6A_data/FA233_exp_all_tx_mean.txt',sep='\t',header=0,index_col=0)
    high_exp=[x for x in exp.index if len([y for y in exp.loc[x,] if y > 1])==66]
    pas_sites=pd.read_csv('APA/human_pas.txt',sep='\t')
    pas_sites=pas_sites[pas_sites['PAS type'].isin(["3'UTR(F)","3'UTR(L)"])]
    # pas_3utr=[x.strip() for x in open('APA/human_pas_3UTR_IDs.txt')]
    # pas_sites=pas_sites[pas_sites['PAS_ID'].isin(pas_3utr)]
    pas_sites=pas_sites[pas_sites['Ensemble ID'].isin(high_exp)]
    if gene_list is not None:
        pas_sites=pas_sites[pas_sites['Ensemble ID'].isin(gene_list)]
    pas_sites['PSE']=pas_sites['PSE'].str.rstrip('%').astype('float') / 100.0
    pas_sites=pas_sites[pas_sites['PSE']>=0.1]
    counter=Counter(list(pas_sites['Ensemble ID']))
    multi_pas=[x for x in counter.keys() if counter[x]>1]
    pas_sites=pas_sites[pas_sites['Ensemble ID'].isin(multi_pas)]
    plus_pas=pas_sites[pas_sites['Strand']=='+']
    #plus_pas.sort_values('PAS_ID',inplace=True,ignore_index=True)
    minus_pas=pas_sites[pas_sites['Strand']=='-']
    #minus_pas.sort_values('PAS_ID',inplace=True,ignore_index=True)
    distal_pas=[]
    proximal_pas=[]
    for gene in set(plus_pas['Ensemble ID']):
        sub=plus_pas[plus_pas['Ensemble ID']==gene]
        if len(sub.index)>=2:
            proximal_pas.append(list(sub.iloc[0]))
            distal_pas.append(list(sub.iloc[len(sub)-1]))

    for gene in set(minus_pas['Ensemble ID']):
        sub=minus_pas[minus_pas['Ensemble ID']==gene]
        if len(sub.index)>=2:
            proximal_pas.append(list(sub.iloc[len(sub)-1]))
            distal_pas.append(list(sub.iloc[0]))

    distal_pas=pd.DataFrame(distal_pas,columns=plus_pas.columns)
    proximal_pas=pd.DataFrame(proximal_pas,columns=plus_pas.columns)
    return proximal_pas,distal_pas

def count_rrach_around_pas(m6A_sites,distance=200,proximal_pas=None,distal_pas=None):
    if proximal_pas is None or distal_pas is None:
        proximal_pas,distal_pas=get_clean_pas(gene_list=m6A_sites['Ensembl ID'])
    proximal_gene_sum=defaultdict(list)
    proximal_position_count=defaultdict(int)
    proximal_peaks=[]
    for i in range(len(proximal_pas)):
        if not i%1000:
            print(i,end=',')
        sites=m6A_sites[m6A_sites['Ensembl ID']==proximal_pas.iloc[i]['Ensemble ID']]
        if len(sites)==0:
            continue
        sites['Distance']=np.add(sites['Position'],-proximal_pas.iloc[i]['Position'])
        sites['Distance']=sites['Distance']*sites['Strand']
        sites=sites[(sites['Distance']<=distance) &(sites['Distance']>=-distance) ]
        if len(sites)==0:
            continue
        sites['ID']=sites['Chromosome']+':'+sites['Position'].astype(str)+':'+sites['Strand'].astype(str)
        proximal_peaks.extend(list(sites['Peak']))
        proximal_gene_sum[proximal_pas.iloc[i]['Ensemble ID']].extend(list(sites['ID']))
        for val in sites['Distance']:
            proximal_position_count[val]+=1
    distal_gene_sum=defaultdict(list)
    distal_position_count=defaultdict(int)
    distal_peaks=[]
    for i in range(len(distal_pas)):
        if not i%1000:
            print(i,end=',')
        sites=m6A_sites[m6A_sites['Ensembl ID']==distal_pas.iloc[i]['Ensemble ID']]
        if len(sites)==0:
            continue
        sites['Distance']=np.add(sites['Position'],-distal_pas.iloc[i]['Position'])
        sites['Distance']=sites['Distance']*sites['Strand']
        sites=sites[(sites['Distance']<=distance) & (sites['Distance']>=-distance) ]
        sites['ID']=sites['Chromosome']+':'+sites['Position'].astype(str)+':'+sites['Strand'].astype(str)
        distal_peaks.extend(list(sites['Peak']))
        distal_gene_sum[distal_pas.iloc[i]['Ensemble ID']].append(list(sites['ID']))
        for val in sites['Distance']:
            distal_position_count[val]+=1
    return proximal_position_count,distal_position_count,proximal_peaks,distal_peaks

def get_proximal_distal_from_bed(peaks,peak2gene):
    motifs,breakdown=scan_motif(peaks)
    breakdown['gene']=[peak2gene.loc[x,'ensg'] for x in breakdown['peak_name']]
    breakdown=breakdown[breakdown['residue']=='A3']
    breakdown=breakdown[['chr','start','gene','strand','peak_name']]
    breakdown.columns=['Chromosome','Position','Ensembl ID','Strand','Peak']
    breakdown['Strand']=[1 if x=='+' else -1 for x in breakdown['Strand']]
    breakdown['Position']=breakdown['Position'].astype(int)
    return count_rrach_around_pas(breakdown)

def chunkify(a, n):
    k, m = len(a) / n, len(a) % n
    return (a[int(i * k + min(i, m)):int((i + 1) * k + min(i + 1, m))] for i in range(n))

def quantile_normalize(df):
    df_sorted = pd.DataFrame(np.sort(df.values,axis=0), 
                             index=df.index, 
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=1)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)

def write_table(file_,table):
    out_file=open(file_,'w')
    out_file.write('\n'.join(['\t'.join([str(y) for y in x]) for x in table])+'\n')
    out_file.close()

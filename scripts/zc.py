###
## @author [Wankun Deng]
## @email [dengwankun@hotmail.com]
## @create date 2022-03-15 17:30:59
## @modify date 2022-03-15 17:30:59
## @desc [description]
###
from scripts.common import *
cmd='''cat gtf/hg19.gtf | grep stop_codon| awk \
'{gsub(/\\"/,"");split($10,a,".");split($12,b,".");gsub(/;/,"");\
printf("%s\\t%s\\t%s\\t%s:%s:%s\\t.\\t%s\\n",$1,$4,$5,a[1],b[1],$18,$7)}'
'''
stop_codon=subprocess.check_output(cmd,shell=True).decode('utf-8')
stop_codon=pd.DataFrame([x.split('\t') for x in stop_codon.strip().split('\n')],
                        columns=['CHR','START','END','GENE','SCORE','STRAND'])
stop_codon['START']=stop_codon['START'].astype(int)
stop_codon['END']=stop_codon['END'].astype(int)
stop_codon['SCORE']=(stop_codon['END']+stop_codon['START'])/2
stop_codon['SCORE']=stop_codon['SCORE'].astype(int)
stop_codon['Sign']=[1 if x=='+' else -1 for x in stop_codon['STRAND']]
stop_codon.loc[stop_codon['Sign']>0,'START']=[x-2000 for x in stop_codon.loc[stop_codon['Sign']>0,'START']]
stop_codon.loc[stop_codon['Sign']<0,'END']=[x+2000 for x in stop_codon.loc[stop_codon['Sign']<0,'END']]
stop_codon.loc[stop_codon['Sign']>0,'END']=[x+8000 for x in stop_codon.loc[stop_codon['Sign']>0,'END']]
stop_codon.loc[stop_codon['Sign']<0,'START']=[max(0,x-8000) for x in stop_codon.loc[stop_codon['Sign']<0,'START']]
stop_codon.drop_duplicates(subset=['CHR','START','END','STRAND'],inplace=True)
stop_codon=BedTool(stop_codon[['CHR','START','END','GENE','SCORE','STRAND']].to_string(header=False,index=False),
                   from_string=True)

def get_peak_around_sc(peak_file,peak2gene=None):   
    peak_bed=pd.read_csv(peak_file,sep='\t',header=None,names=['CHR','START','END','GENE','SCORE','STRAND'])
    peak_bed['START']=peak_bed['START'].astype('str')
    peak_bed['END']=peak_bed['END'].astype('str')
    peak_bed['peak_id']=peak_bed['CHR']+':'+peak_bed['START']+':'+peak_bed['END']+':'+peak_bed['STRAND']
    
    p2g_cmd='''bedtools intersect -a {peak_file} -b gtf/genes.bed -f 1 -wa -wb -s | \
        awk '{{printf("%s:%s:%s:%s\\t%s\\n",$1,$2,$3,$6,$10)}}'
            '''
    if peak2gene is None:
        peak2gene=pd.read_csv(StringIO(subprocess.check_output(p2g_cmd.format(peak_file=peak_file),
                                                                   shell=True).decode('utf-8').strip()),
                                  sep='\t',header=None,names=['peak','gene'])  
        peak2gene.drop_duplicates(subset=['peak'],inplace=True)
        peak2gene.index=peak2gene['peak']
    else:
        peak2gene=pd.read_csv(peak2gene,sep='\t',header=None,index_col=0,names=['gene'])
        
    peak_bed['GENE']=list(peak2gene.reindex(peak_bed['peak_id'])['gene'])
    peak_bed=BedTool(peak_bed[['CHR','START','END','GENE','SCORE','STRAND']].to_string(header=False,index=False),
                       from_string=True)
    peak_around_sc=stop_codon.intersect(peak_bed,wa=True,wb=True)
    peak_around_sc=pd.DataFrame([x.strip().split('\t') for x in str(peak_around_sc).strip().split('\n')],
                               columns=(list('ABCDEFGHIJKL')))
    peak_around_sc=peak_around_sc.loc[[ peak_around_sc.iloc[i,9] in  peak_around_sc.iloc[i,3] \
                                     for i  in range(len(peak_around_sc))],:]
    peak_around_sc.iloc[:,1]=peak_around_sc.iloc[:,1].astype(int)
    peak_around_sc.iloc[:,2]=peak_around_sc.iloc[:,2].astype(int)
    peak_around_sc.iloc[:,4]=peak_around_sc.iloc[:,4].astype(int)
    peak_around_sc.iloc[:,7]=peak_around_sc.iloc[:,7].astype(int)
    peak_around_sc.iloc[:,8]=peak_around_sc.iloc[:,8].astype(int)
    peak_around_sc.iloc[:,5]=[1 if x=='+' else -1 for x in peak_around_sc.iloc[:,5]]
    peak_around_sc.iloc[:,11]=[1 if x=='+' else -1 for x in peak_around_sc.iloc[:,11]]
    peak_around_sc.iloc[:,10]=(peak_around_sc.iloc[:,7]+peak_around_sc.iloc[:,8])/2
    peak_around_sc.iloc[:,10]=peak_around_sc.iloc[:,10].astype(int)
    peak_around_sc['M']=(peak_around_sc.iloc[:,10]-peak_around_sc.iloc[:,4])*peak_around_sc.iloc[:,11]
    peak_around_sc['N']=[x.split(':')[1] for x in peak_around_sc['D']]
    ret=defaultdict(list)
    for i in range(len(peak_around_sc)):
        ret[peak_around_sc.iloc[i,13]].append(peak_around_sc.iloc[i,12])
    return peak_around_sc,ret

def count_ol(RBP,module,m6a_dict):
    if not 'K562' in RBP:
        return
    print(RBP,module)
    rbp_around_sc,rbp_dict=get_peak_around_sc('WGCNA/RBP_peaks/ENCODE_peaks/%s-human.sorted.bed'%RBP)
    rbp_dict={x:rbp_dict[x] for x in [y for y in rbp_dict.keys() if y in m6a_dict]}
    m6a_dict2={x:m6a_dict[x] for x in [y for y in rbp_dict.keys() if y in m6a_dict]}
    max_ele=0
    for key in rbp_dict:
        max_ele=max(len(m6a_dict2[key]),max(len(rbp_dict[key]),max_ele))
    if len(m6a_dict2)==0:
        return
    table='TX\t'
    for list_ in ['m6A','RBP']:
        for i  in range(1,max_ele+1):
            table+='%s%s\t'%(list_,i)
    table+='count1\tcount2\n'
    for tx in m6a_dict2:
        table+=tx+'\t'
        for list_ in [m6a_dict2,rbp_dict]:
            for i  in range(max_ele):
                if len(list_[tx])<=i:
                    table+='NA\t'
                else:
                    table+='%s\t'%list_[tx][i]
        table+='%s\t%s\n'%(len(m6a_dict2[tx]),len(rbp_dict[tx]))
    if not os.path.isdir('WGCNA/zc'):
        os.mkdir('WGCNA/zc')
    with open('WGCNA/zc/%s_%s.txt'%(RBP,module),'w') as fo:
        fo.write(table.strip())
        fo.close()

def count_module_ol(module):
    _,m6a_dict=get_peak_around_sc('WGCNA/'+module+'_peak.bed',peak2gene='m6A_data/peak2gene.txt')
    for file_ in os.listdir('WGCNA/RBP_peaks/ENCODE_peaks'):
        if file_.endswith('-human.sorted.bed'):
            rbp=file_.replace('-human.sorted.bed','')
            count_ol(rbp,module,m6a_dict)

pool=Pool(processes=25)
ret=pool.map(count_module_ol,['M%s'%x for x in range(1,26)])
pool.close()
pool.join()

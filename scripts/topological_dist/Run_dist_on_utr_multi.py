###
## @author [Wankun Deng]
## @email [dengwankun@hotmail.com]
## @create date 2022-04-03 23:13:43
## @modify date 2022-04-03 23:13:43
## @desc [description]
###
import sys
import pandas as pd
import subprocess
import numpy as np
bin_len=50
script=sys.argv[1]
fo=open(sys.argv[2],'w')
groups=[]
fins=[]
for i in range(int((len(sys.argv)-3)/2)):
    groups.append(sys.argv[3+i*2])
    fins.append(sys.argv[4+i*2])

for i in range(len(fins)):
    cmd='python %s %s hg19 %s'%(script,fins[i],50)
    info=subprocess.check_output(cmd,shell=True).decode('utf-8').strip().split('\n')
    for j in range(len(info)):
        fo.write('%s\t%s\t%s\n'%(groups[i],j+1,info[j].split('\t')[0]))
fo.close()

df=pd.read_csv(sys.argv[2],sep='\t',header=None,names=['Group','X','Y'])
df['Accumulate']=np.zeros(bin_len*3*len(groups))
for i in range(len(groups)):
    df.loc[df['Group']==groups[i],'Y']=list(df[df['Group']==groups[i]]['Y']/sum(df[df['Group']==groups[i]]['Y']))
    df.iloc[i*bin_len*3:(i+1)*bin_len*3,3]=[sum(df.iloc[i*bin_len*3:i*bin_len*3+x+1,2]) for x in range(bin_len*3)]
df.to_csv(sys.argv[2],sep='\t',index=False)
###
## @author [Wankun Deng]
## @email [dengwankun@hotmail.com]
## @create date 2022-03-15 17:30:42
## @modify date 2022-03-15 17:30:42
## @desc [description]
###
import os
import pandas as pd
import subprocess
import re
from multiprocessing.pool import Pool
import multiprocessing
from collections import defaultdict
import statistics
import random
from functools import partial
os.chdir('../')
peak_intensity=pd.read_csv('m6A_data/peak_intensity_m_1_renamed.csv',sep='\t',index_col=0)
#colnames=list(set(['_'.join(x.split('_')[:-1]) for x in list(peak_intensity.columns)]))

valid_idx=[]
[valid_idx.append(x) for x in peak_intensity.columns if 'Fetal' not in x]

colnames=list(set(['_'.join(x.split('_')[:-1]) for x in valid_idx]))
print(colnames)
random.seed(100)
def tau_index(site_in):
    return sum([1-(x/max(site_in)) for x in site_in])/len(site_in)

def select(colnames):
    return [x+'_'+str(random.randrange(1,3,1)) for x in colnames]


def chunkify(a, n):
    k, m = len(a) / n, len(a) % n
    return [a[int(i * k + min(i, m)):int((i + 1) * k + min(i + 1, m))] for i in range(n)]

count=0
def print_progress():
    count+=100
    if not count%1000:
        print(count)
        
def calc_a_list(args):
    (t_,list_)=args
    tau_={}
    count_=0
    for site in list_:
        count_+=1
        tau_list=[]
        for i in range(101):
            tau_list.append(tau_index(list(peak_intensity.loc[site,select(colnames)])))
        tau_[site]=statistics.median(tau_list)
        if not count_%1000:
            print('Thread %s:%s'%(t_,count_))
#             lock.acquire()
#             print_progress()
#             lock.release()
    return tau_

tau={}
count=0
site_list_list=chunkify(list(peak_intensity.index),20)

to_map=[(x,site_list_list[x]) for x in range(len(site_list_list))]

pool=Pool(processes=20)
# m = multiprocessing.Manager()
# l = m.Lock()
# func = partial(calc_a_list, l)

dict_list=pool.map(calc_a_list,to_map)
pool.terminate()
pool.join()

for dict_ in dict_list:
    tau.update(dict_)
pd.Series(tau).to_frame('tau_index').to_csv('ts_data/tau_index_rand.txt')

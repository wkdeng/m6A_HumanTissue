###
## @author [Wankun Deng]
## @email [dengwankun@gmail.com]
## @create date 2023-03-06 01:13:00
## @modify date 2023-03-06 01:13:00
## @desc [description]
###
import os
import sys
from collections import defaultdict

unique_rate=defaultdict(float)
for dir in os.listdir(sys.argv[1]):
    file_path=os.path.join(sys.argv[1], dir, 'mapping_stats.txt')
    if os.path.exists(file_path):
        for line in open(file_path):
            if line.startswith('% of reads unmapped: too short'):
                info = line.strip().split('\t')[1][:-1]
                unique_rate['_'.join(dir.split('_')[:2])] = float(info) / 100
                break
for key in unique_rate:
    print( key+'\t' + str(unique_rate[key]))
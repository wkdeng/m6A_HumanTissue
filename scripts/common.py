###
## @author [Wankun Deng]
## @email [dengwankun@gmail.com]
## @create date 2023-03-06 01:13:11
## @modify date 2023-03-06 01:13:11
## @desc [description]
###
import re
import os
import csv
import sys
import math
import json
import scipy
import pysam
import shutil
import random
import pyBigWig
import warnings
import subprocess
import statistics
import multiprocessing
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as mt
from io import StringIO
from scipy import stats
from scripts import hyper_geometry
from functools import partial
from itertools import product
from pybedtools import BedTool
from collections import defaultdict
from collections import Counter
from scipy.stats import chisquare
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
from scipy.stats.mstats import gmean
from multiprocessing.pool import Pool
from sklearn.linear_model import LogisticRegression

random.seed(8923)
tissues=['AAT', 'AMO', 'FWH', 'ACS', 'AMG', 'FLU', 'ATR', 'FBR', 'ASC', 'ATY', 'ALE', 'ATA', 'ABL', 'ADU', 'AKI', 'AAD', 
         'ASO', 'AAC', 'FLI', 'AAO', 'ASG', 'FKI', 'APR', 'ABR', 'AHI', 'ALV', 'FAG', 'AEP', 'ACO', 'APA', 'AST', 'APE', 
         'APL', 'APG', 'AUT', 'ALN', 'ATE', 'ARV', 'ALI', 'AIM', 'ATL', 'FSP', 'ACL', 'ACE', 'FSC', 'FHE', 'FTH', 'ACN', 
         'AAG', 'ACA', 'AES', 'ACC', 'APO', 'ASM', 'ASI', 'ATN', 'AOV', 'ASP', 'AHE', 'ATG', 'AFL', 'AAP', 'AAS', 'AIC', 
         'ALU', 'AJE']
matched_gtex=['Brain_Cerebellum','Brain_Cortex','Brain_Hippocampus','Pituitary','Spleen','Thyroid','Adrenal_Gland',
                'Pancreas','Lung','Uterus','Ovary','Artery_Aorta','Heart_Left_Ventricle','Stomach','Prostate','Testis',
                'Liver','Muscle_Skeletal']
matched_ht=['NervousSystem_Adult_Cerebellum','NervousSystem_Adult_CerebralCortex','NervousSystem_Adult_Hippocampus',
              'GlandularTissues_Adult_PituitaryGland','GlandularTissues_Adult_Spleen','GlandularTissues_Adult_Thyroid',
              'GlandularTissues_Adult_AdrenalGland','GlandularTissues_Adult_Pancreas','RespiratorySystem_Adult_Lung',
              'ReproductiveSystem_Adult_Uterus','ReproductiveSystem_Adult_Ovary','CirculatorySystem_Adult_Aorta',
              'CirculatorySystem_Adult_LeftVentricle','DigestiveSystem_Adult_Stomach', 'ReproductiveSystem_Adult_Prostate',
              'ReproductiveSystem_Adult_Testis','Other_Adult_Liver','Other_Adult_SkeletalMuscle']
sample_tissue=['ABR','ACE','ACC','ACN','AHI','APO','ASC','AHE','ARV','AAD','AAO','AST','ADU','ASI','ACL']
tissue_full={'AAO':'CirculatorySystem_Adult_Aorta', 'AAD':'CirculatorySystem_Adult_AuricleDextra', 
             'AAS':'CirculatorySystem_Adult_AuricleSinistra', 'AHE':'CirculatorySystem_Adult_Heart', 
             'ALV':'CirculatorySystem_Adult_LeftVentricle', 'APE':'CirculatorySystem_Adult_Pericardium', 
             'ARV':'CirculatorySystem_Adult_RightVentricle', 'FHE':'CirculatorySystem_Fetal_Heart', 
             'AAP':'DigestiveSystem_Adult_Appendix', 'ACA':'DigestiveSystem_Adult_Ascendingcolon', 
             'ACL':'DigestiveSystem_Adult_Colon', 'ADU':'DigestiveSystem_Adult_Duodenum', 
             'AES':'DigestiveSystem_Adult_Esophagus', 'AIC':'DigestiveSystem_Adult_Ileocecum', 
             'AIM':'DigestiveSystem_Adult_Ileum', 'AJE':'DigestiveSystem_Adult_Jejunum', 
             'ASI':'DigestiveSystem_Adult_SmallIntestine', 'AST':'DigestiveSystem_Adult_Stomach', 
             'ACS':'DigestiveSystem_Adult_StomachCardia', 'ASO':'DigestiveSystem_Adult_StomachCorpus', 
             'ATG':'DigestiveSystem_Adult_Tongue', 'AAC':'GlandularTissues_Adult_AdrenalCortex', 
             'AAG':'GlandularTissues_Adult_AdrenalGland', 'AMG':'GlandularTissues_Adult_MammaryGland', 
             'APA':'GlandularTissues_Adult_Pancreas', 'APG':'GlandularTissues_Adult_PituitaryGland', 
             'ASG':'GlandularTissues_Adult_SalivaryGland', 'ASP':'GlandularTissues_Adult_Spleen', 
             'ATY':'GlandularTissues_Adult_Thyroid', 'ATN':'GlandularTissues_Adult_Tonsil', 
             'FAG':'GlandularTissues_Fetal_AdrenalGland', 'FSP':'GlandularTissues_Fetal_Spleen', 
             'FTH':'GlandularTissues_Fetal_Thymus', 'ALN':'ImmuneSystem_Adult_LymphNode', 
             'ALE':'ImmuneSystem_Adult_PeripheralLeukocytes', 'ABR':'NervousSystem_Adult_Brain', 
             'ACN':'NervousSystem_Adult_CaudateNucleus', 'ACE':'NervousSystem_Adult_Cerebellum', 
             'ACC':'NervousSystem_Adult_CerebralCortex', 'ACO':'NervousSystem_Adult_CorpusCallosum', 
             'AFL':'NervousSystem_Adult_FrontalLobe', 'AHI':'NervousSystem_Adult_Hippocampus', 
             'AMO':'NervousSystem_Adult_MedullaOblongata', 'APO':'NervousSystem_Adult_Pons', 
             'ASC':'NervousSystem_Adult_SpinalCord', 'ATL':'NervousSystem_Adult_TemporalLobe', 
             'ATA':'NervousSystem_Adult_Thalamus', 'FBR':'NervousSystem_Fetal_Brain', 
             'FSC':'NervousSystem_Fetal_SpinalCord', 'AAT':'Other_Adult_AdiposeTissue', 
             'ALI':'Other_Adult_Liver', 'ASM':'Other_Adult_SkeletalMuscle', 'FWH':'Other_Fetal_Fetus', 
             'FLI':'Other_Fetal_Liver', 'AEP':'ReproductiveSystem_Adult_Epididymis', 
             'AOV':'ReproductiveSystem_Adult_Ovary', 'APL':'ReproductiveSystem_Adult_Placenta', 
             'APR':'ReproductiveSystem_Adult_Prostate', 'ATE':'ReproductiveSystem_Adult_Testis', 
             'AUT':'ReproductiveSystem_Adult_Uterus', 'ALU':'RespiratorySystem_Adult_Lung', 
             'ATR':'RespiratorySystem_Adult_Trachea', 'FLU':'RespiratorySystem_Fetal_Lung', 
             'ABL':'UrinarySystem_Adult_Bladder', 'AKI':'UrinarySystem_Adult_Kidney', 'FKI':'UrinarySystem_Fetal_Kidney'}
tissue_abbr={tissue_full[x]:x for x in tissues}
fasta='gtf/hg19.fasta'

from scripts.utils import *

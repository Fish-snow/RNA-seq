import pandas as pd

v1 = pd.read_csv("/Users/yuxueqian/Downloads/GSE175396_RAW/GSM5332640_Vehicle_S1.counts.txt", sep='\t',keep_default_na=False,na_values=['YXQ'])
v2 = pd.read_csv("/Users/yuxueqian/Downloads/GSE175396_RAW/GSM5332641_Vehicle_S2.counts.txt", sep='\t')
a1 = pd.read_csv("/Users/yuxueqian/Downloads/GSE175396_RAW/GSM5332642_Ethanol_08_S1.counts.txt", sep='\t')
a2 = pd.read_csv("/Users/yuxueqian/Downloads/GSE175396_RAW/GSM5332643_Ethanol_08_S2.counts.txt", sep='\t')
a3 = pd.read_csv("/Users/yuxueqian/Downloads/GSE175396_RAW/GSM5332644_Ethanol_08_S3.counts.txt", sep='\t')
a4 = pd.read_csv("/Users/yuxueqian/Downloads/GSE175396_RAW/GSM5332645_Ethanol_16_S1.counts.txt", sep='\t')
a5 = pd.read_csv("/Users/yuxueqian/Downloads/GSE175396_RAW/GSM5332646_Ethanol_16_S2.counts.txt", sep='\t')
m1 = pd.merge(v1, v2, how='outer', on='Geneid')
m1 = pd.merge(m1, a1, how='outer', on='Geneid')
m1 = pd.merge(m1, a2, how='outer', on='Geneid')
m1 = pd.merge(m1, a3, how='outer', on='Geneid')
m1 = pd.merge(m1, a4, how='outer', on='Geneid')
m1 = pd.merge(m1, a5, how='outer', on='Geneid')

del a1,a2,a3,a4,a5,v1,v2
m1['sum']=m1.iloc[:,1:9].sum(axis=1)
filter_m1=m1.query('sum!=0')
filter_m1=filter_m1.drop('sum',axis=1)

filter_m1=filter_m1.T
m1=filter_m1
del filter_m1
m1.columns=m1.loc['Geneid'].tolist()
m1=m1.drop('Geneid',axis=0)

import os
import pickle as pkl

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

SAVE = True
if SAVE:
    OUTPUT_PATH ="/Users/yuxueqian/Desktop/"
    os.makedirs(OUTPUT_PATH, exist_ok=True)

metadata=pd.read_excel("/Users/yuxueqian/Desktop/metadata.xlsx")
metadata=metadata.set_index('Unnamed: 0')


#单因素分析
inference = DefaultInference(n_cpus=8)
dds = DeseqDataSet(counts=m1,
                   metadata=metadata,
                   design_factors="group",
                   refit_cooks=True,
                   inference=inference)


dds.deseq2()
if SAVE:
    with open(os.path.join(OUTPUT_PATH, "dds.pkl"), "wb") as f:
        pkl.dump(dds, f)

stat_res = DeseqStats(dds, inference=inference)
stat_res.summary()
if SAVE:
    with open(os.path.join(OUTPUT_PATH, "stat_results.pkl"), "wb")as f:
        pkl.dump(stat_res, f)

stat_res.lfc_shrink(coeff="group_v_vs_t")
if SAVE:
    with open(os.path.join(OUTPUT_PATH, "stat_results.pkl"), "wb")as f:
        pkl.dump(stat_res, f)
stat_res.results_df.to_csv(os.path.join(OUTPUT_PATH, "stat_results.csv"))

#匹配
match=pd.read_csv("/Users/yuxueqian/Desktop/mart_export.txt",sep='\t')
result=pd.merge(stat_res.results_df,match,left_index=True,right_on='Gene stable ID')
result.to_csv(os.path.join(OUTPUT_PATH, "results.csv"))


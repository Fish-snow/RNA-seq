import numpy as np
import pandas as pd
import sns

#制作差异分析结果数据框
genearray=np.asarray(pvalue)
result=pd.DataFrame({'pvalue':genearray,'FoldChange':fold})
result['log(pvalue)']=-np.log10(result['pvalue'])

#准备
#选定差异基因标准I.差异倍数绝对值大于1，II.差异分析P值小于0.05
result['sig']='normal'
result['size']=np.abs(result['FoldChange'])/10

result.loc[(result.FoldChange>1)&(result.pvalue<0.05),'sig']='up'
result.loc[(result.FoldChange<-1)&(result.pvalue<0.05),'sig']='down'

ax=sns.scatterplot(x="FoldChange",y="log(pvalue)"),
                   hue='sig',
                   hue_order=('down','normal','up'),
                   palette=("#377EB8",'grey',"#E41A1C"),
                   data=result)
ax.set_ylabel('-log(pvalue)',fontweight='bold')
ax.set_xlabel('FoldChange',fontweight='bold')

#筛选差异基因
fold_cutoff=1
pvalue_cutoff=0.05

filterd_ids=[]
for i in range(0,number_of_genes):
    if(abs(fold[i]>=fold_cutoff)and(pvalue[i]<=pvalue_cutoff):
       filtered_ids.append(i)

filtered=data2.iloc[filtered_ids,:]
print("Number of DE genes:")
print(len(filterd_idex))

#heatmap
sns.clustermap(filtered,cmap='RdYlGn_r',standard_scale=0)
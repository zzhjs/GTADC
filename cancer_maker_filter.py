import pandas as pd
from pathlib import Path
from _api import ScMags
adata_cortex = pd.read_csv('./data/scRNA_data.csv', index_col=0).T
adata_cortex_meta = pd.read_csv('./data/scRNA_meta.csv', index_col=0)
data=adata_cortex.to_numpy()
gene_ann=adata_cortex.columns.values
labels=adata_cortex_meta['type'].to_numpy()
labels = labels.reshape(data.shape[0])
gene_ann = gene_ann.reshape(data.shape[1])
obj = ScMags(data, labels, gene_ann, verbose = False)
obj.filter_genes(nof_sel=200)
obj.sel_clust_marker(
        nof_markers = 30, 
        n_cores=20
    ) 
mark_ann = obj.get_markers(ind_return=False)
mark_ann.to_csv('genelist.csv')
markers_res=pd.read_csv('genelist.csv', index_col=0).T
markers_res_=[]
for column in markers_res: 
    markers_res_.extend(markers_res[column].tolist())
markers_res_ = list(set(markers_res_))
with open('list.txt', 'w') as f:
    # 将列表中的每个元素转换为字符串，并逐行写入文件
    for item in set(markers_res_):
        f.write(str(item) + '\n')
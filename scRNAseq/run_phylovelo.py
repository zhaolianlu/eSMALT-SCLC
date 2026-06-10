import anndata as ad
import phylovelo as pv
import numpy as np
import pandas as pd
import scanpy as sc
import re
import pickle

adata = ad.read_h5ad('/data3/wangkun/zhoalian_jlab/adata_ent.h5ad')
drop_genes = []
for i in adata.var_names:
    for num in re.finditer('[0-9]+', i):
        if len(num.group(0)) >= 4:
            drop_genes.append(i)
            continue
sel_genes = adata.var_names[~np.isin(adata.var_names, drop_genes)]
sel_genes = sel_genes[:-3]

adata = adata[:, sel_genes]
adata = adata[adata.obs['numMut']!='Others']

from tqdm.notebook import tqdm
for c in tqdm(set(adata.obs['clone'])):
    tmp = adata[adata.obs['clone']==c]
    tmp.write_h5ad(f'/data3/wangkun/zhoalian_jlab/clone_spec/adata_ent_{c}.h5ad')

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
count = pd.DataFrame(adata.X.A, index=adata.obs_names, columns=adata.var_names)
sd = pv.scData()
sd.count = count
sd.count = sd.count[sd.count.columns[sd.count.sum()>0]]
sd.drop_duplicate_genes()
sd.x_normed = sd.count
mutnum = np.array([int(i) for i in adata.obs['numMut']])
pv.velocity_inference(sd, mutnum, cutoff=0.97, target='x_normed')
v_megs = pd.DataFrame(data=sd.velocity, index=sd.x_normed.columns, columns=['velocity']).loc[sd.megs]
v_megs['p_value']=sd.pvals.T
v_megs['q_value']=sd.qvals.T
v_megs.to_csv(f'./megs_ent.csv')

sd.Xdr = adata.obs[['umapharmony_1','umapharmony_2']]
pv.velocity_embedding(sd, target='x_normed', n_neigh=min(1000, sd.Xdr.shape[0]//3))
pv.calc_phylo_pseudotime(sd, r_sample=0.03)
sd =  pickle.load(open('/data3/wangkun/zhoalian_jlab/sd_ent.pkl', 'rb'))

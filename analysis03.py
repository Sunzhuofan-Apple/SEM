import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cellxgene_census
import scanpy as sc

data = sc.read_h5ad("5d871206-9489-4d9f-8106-94305ccb1c3a.h5ad")  
print(data)
gene_list = ["ENSG00000091138", "ENSG00000138670", "ENSG00000257647", "ENSG00000169783", 
             "ENSG00000178104", "ENSG00000237837", "ENSG00000224363", 
             "ENSG00000266302", "ENSG00000253103", "ENSG00000168490", "ENSG00000251555"]
data_subset = data[:, gene_list]
print(data_subset)


male_group = data[data.obs['sex'] == 'male']
female_group = data[data.obs['sex'] == 'female']


sc.tl.rank_genes_groups(data, 'sex', groups=['male'], reference='female', method='t-test')


result = data.uns['rank_genes_groups']
genes = result['names']  
scores = result['scores']  


top_genes = genes[:10]
print("Top differentially expressed genes between male and female:", top_genes)


sc.pl.rank_genes_groups_heatmap(data, n_genes=10, groupby='sex', show_gene_labels=True)

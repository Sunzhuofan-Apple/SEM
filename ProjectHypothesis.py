import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import cellxgene_census
import scanpy as sc
import matplotlib.pyplot as plt

data = sc.read_h5ad("5d871206-9489-4d9f-8106-94305ccb1c3a.h5ad")  
print(data)
gene_list = ["ENSG00000091138", "ENSG00000138670", "ENSG00000257647", "ENSG00000169783", 
             "ENSG00000178104", "ENSG00000237837", "ENSG00000224363", 
             "ENSG00000266302", "ENSG00000253103", "ENSG00000168490", "ENSG00000251555"]
data_subset = data[:, gene_list]
print(data_subset)
# 将 Age 列转换为数值格式
data.obs['Age'] = pd.to_numeric(data.obs['Age'], errors='coerce')

# 创建年龄组：将 <= 65 的样本标记为 'young'，其他为 'old'
data.obs['age_group'] = data.obs['Age'].apply(lambda x: 'young' if x <= 65 else 'old')
data.obs['sex'] = data.obs['sex'].astype(str)  # 转换为字符串，确保数据类型一致
data.obs['sex'] = data.obs['sex'].astype('category')  # 转换为 categorical 类型

# 检查是否包含 'male' 和 'female' 两个类别
print(data.obs['sex'].cat.categories)

# 差异表达分析：性别（male vs female）和年龄（young vs old）比较
sc.tl.rank_genes_groups(data, 'sex', groups=['male'], reference='female', method='t-test')
sc.tl.rank_genes_groups(data, 'age_group', groups=['young'], reference='old', method='t-test')

# 打印前10个性别差异表达基因
top_genes_sex = data.uns['rank_genes_groups']['names'][:10]
print("Top differentially expressed genes by sex:", top_genes_sex)

# 打印前10个年龄差异表达基因
top_genes_age = data.uns['rank_genes_groups']['names'][:10]
print("Top differentially expressed genes by age:", top_genes_age)

# 可视化性别和年龄差异基因的热图
sc.tl.dendrogram(data, groupby='sex')
sc.tl.dendrogram(data, groupby='age_group')

# 创建图像并绘制两个热图和一个小提琴图
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# 绘制性别差异的热图
sc.pl.rank_genes_groups_heatmap(data, n_genes=10, groupby='sex', show_gene_labels=True, ax=axes[0])
axes[0].set_title('Top Genes by Sex')

# 绘制年龄组差异的热图
sc.pl.rank_genes_groups_heatmap(data, n_genes=10, groupby='age_group', show_gene_labels=True, ax=axes[1])
axes[1].set_title('Top Genes by Age Group')

# 绘制额外的小提琴图来展示特定基因的表达情况
sc.pl.violin(data, keys=['ENSG00000134569'], groupby='sex', ax=axes[2])
axes[2].set_title('Expression of ENSG00000134569 by Sex')

# 调整布局
plt.tight_layout()
plt.show()
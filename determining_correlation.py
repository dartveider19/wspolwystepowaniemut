import GEOparse
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np


files_folder_dir = "C:/Users/micle/Desktop/Python/dyplom/Practice/"



gse72217 = GEOparse.get_GEO(filepath="./GSE72217_family.soft.gz")
samples = list(gse72217.gsms.keys())
samples_expression = gse72217.pivot_samples('VALUE')[samples]
annotation = gse72217.gpls['GPL5175'].table
annotation = annotation.rename(columns={'ID':'ID_REF'})
gene_ids = list(samples_expression.index)
annotation = annotation[annotation['ID_REF'].isin(gene_ids)]
data = pd.concat([annotation, samples_expression], axis=1)
data.to_csv('C:/Users/micle/Desktop/Python/dyplom/GSE72217_expression_data.csv')

data = pd.read_csv("C:/Users/micle/Desktop/Python/dyplom/Practice/GSE72217_expression_data.csv", low_memory=False)

main_data = data.loc[data['category'] == 'main']

gene_names = []

for gene_name in main_data['gene_assignment']:
    if len(gene_name) > 6:
        gene_names.append(gene_name.rsplit('//')[1])
    else:
        gene_names.append(gene_name)

main_data['gene_name'] = gene_names
#main_data.insert(location, col_name, col_values)

main_data.to_csv(files_folder_dir + 'data_with_gene_names.csv')
main_data['gene_name'] = main_data['gene_name'].str.strip()

#clean_data = main_data.drop(main_data.columns[range(1, 13)], axis = 1)

clean_data = main_data.drop(main_data.columns[range(2,13)], axis = 1)
clean_data = clean_data.drop(clean_data.columns[[0]], axis = 1)
clean_data = clean_data.set_index("ID_REF")
clean_data.index = clean_data.index.map(str)
clean_data = clean_data.T



correlation_data = pd.DataFrame()
correlation_data['gene_name'] = main_data['gene_name']


pearson = []
spearman = []
kendall = []


egfr = clean_data['3002640']
TTLL10 = clean_data['2315554']

for gene_name in clean_data:
    a = clean_data[gene_name]
    rho = egfr.corr(a, method = 'pearson')
    cor = egfr.corr(a, method = 'spearman')
    tau = egfr.corr(a, method = 'kendall')
    pearson.append(rho)
    spearman.append(cor)
    kendall.append(tau)
    print(gene_name)
    print(rho, cor, tau)

correlation_data['pearson'] = pearson
correlation_data['spearman'] = spearman
correlation_data['kendall'] = kendall
correlation_data['ID_REF'] = main_data['ID_REF']

correlation_data.to_csv(files_folder_dir + 'correlation_data.csv')

#Visualization example
slope, intercept, r, p, stderr = scipy.stats.linregress(egfr, TTLL10)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fig, ax = plt.subplots()
ax.plot(TTLL10, egfr, linewidth=0, marker='s', label='Data points')
ax.plot(TTLL10, intercept + slope * TTLL10, label=line)
ax.set_xlabel('TTLL10 gene expression')
ax.set_ylabel('EGFR gene expression')
ax.legend(facecolor='white')
plt.show()

#Visualization of the 5 most correlated genes to EGFR
SEC61G = clean_data['3051395']
NCOA1 = clean_data['2473149']
INNPP4A = clean_data['2495446']
CACHD1 = clean_data['2340078']
INO80D = clean_data['2596162']



slope, intercept, r, p, stderr = scipy.stats.linregress(egfr, SEC61G)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fig, ax = plt.subplots()
ax.plot(SEC61G, egfr, linewidth=0, marker='s', label='Data points')
ax.plot(SEC61G, intercept + slope * SEC61G, label=line)
ax.set_xlabel('SEC61G gene expression')
ax.set_ylabel('EGFR gene expression')
ax.legend(facecolor='white')
plt.show()



slope, intercept, r, p, stderr = scipy.stats.linregress(egfr, NCOA1)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fig, ax = plt.subplots()
ax.plot(NCOA1, egfr, linewidth=0, marker='s', label='Data points')
ax.plot(NCOA1, intercept + slope * NCOA1, label=line)
ax.set_xlabel('NCOA1 gene expression')
ax.set_ylabel('EGFR gene expression')
ax.legend(facecolor='white')
plt.show()


slope, intercept, r, p, stderr = scipy.stats.linregress(egfr, INNPP4A)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fig, ax = plt.subplots()
ax.plot(INNPP4A, egfr, linewidth=0, marker='s', label='Data points')
ax.plot(INNPP4A, intercept + slope * INNPP4A, label=line)
ax.set_xlabel('INNPP4A gene expression')
ax.set_ylabel('EGFR gene expression')
ax.legend(facecolor='white')
plt.show()



slope, intercept, r, p, stderr = scipy.stats.linregress(egfr, CACHD1)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fig, ax = plt.subplots()
ax.plot(CACHD1, egfr, linewidth=0, marker='s', label='Data points')
ax.plot(CACHD1, intercept + slope * CACHD1, label=line)
ax.set_xlabel('CACHD1 gene expression')
ax.set_ylabel('EGFR gene expression')
ax.legend(facecolor='white')
plt.show()



slope, intercept, r, p, stderr = scipy.stats.linregress(egfr, INO80D)
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

fig, ax = plt.subplots()
ax.plot(INO80D, egfr, linewidth=0, marker='s', label='Data points')
ax.plot(INO80D, intercept + slope * INO80D, label=line)
ax.set_xlabel('INO80D gene expression')
ax.set_ylabel('EGFR gene expression')
ax.legend(facecolor='white')
plt.show()

correlation_data_improved = pd.DataFrame()
correlation_data_improved['gene_name'] = main_data['gene_name']
correlation_data_improved['ID_REF'] = main_data['ID_REF']

pearson_improved = []
pearson_p_value = []
spearman_improved = []
spearman_p_value = []
kendall_improved = []
kendall_p_value = []

for id_ref in clean_data:
    y = egfr
    x = clean_data[id_ref]
    pearson_improved.append(scipy.stats.pearsonr(x, y)[0])
    pearson_p_value.append(scipy.stats.pearsonr(x, y)[1])
    spearman_improved.append(scipy.stats.spearmanr(x, y)[0])
    spearman_p_value.append(scipy.stats.spearmanr(x, y)[1])
    kendall_improved.append(scipy.stats.kendalltau(x, y)[0])
    kendall_p_value.append(scipy.stats.kendalltau(x, y)[1])
    print(id_ref)
    print(scipy.stats.pearsonr(x, y)[0], scipy.stats.pearsonr(x, y)[1])
    print(scipy.stats.spearmanr(x, y)[0], scipy.stats.spearmanr(x, y)[1])
    print(scipy.stats.kendalltau(x, y)[0], scipy.stats.kendalltau(x, y)[1])

correlation_data_improved['Pearson r'] = pearson_improved
correlation_data_improved['Pearson p-value'] = pearson_p_value
correlation_data_improved['Spearman rho'] = spearman_improved
correlation_data_improved['Spearman p-value'] = spearman_p_value
correlation_data_improved['Kendall tau'] = kendall_improved
correlation_data_improved['Kendall p-value'] = kendall_p_value

correlation_data_improved.to_csv(files_folder_dir + 'correlation_data_improved_p_value.csv')
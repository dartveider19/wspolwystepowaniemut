import GEOparse
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np


##data_unpack("GSE72217_family.soft.gz", "C:/Users/micle/Desktop/Python/dyplom")

#path = the folder where is your file
##name = the name of the file in that folder
###change backslashes to slashes


#extracts data from the GEO file (only downloaded file)
## name = name of the file
## path = direct path to the folder where file is located
def data_unpack(name, path):
    global data

    gse = GEOparse.get_GEO(filepath= path + '/' + name)
    samples = list(gse.gsms.keys())
    samples_expression = gse.pivot_samples('VALUE')[samples]
    for k in gse.gpls:
        gpls_number = k
    annotation = gse.gpls[gpls_number].table
    annotation = annotation.rename(columns={'ID': 'ID_REF'})
    gene_ids = list(samples_expression.index)
    annotation = annotation[annotation['ID_REF'].isin(gene_ids)]
    data = annotation.join(samples_expression, on = 'ID_REF')
    data.to_csv(path + '/' + name + '_unpacked_total_data.csv')

    return data

#filters control genes out of the main df
def data_controls ():
    global controls

    controls = data.loc[data['SPOT_ID'] == 'control']

    return controls

#filters positive controls out of the controls df
def data_controls_pos ():
    global controls_pos

    controls_pos = controls[controls['mrna_assignment'].str.contains('pos')]

    return controls_pos

#filters negative controls out of the controls df
def data_controls_neg():
    global controls_neg

    controls_neg = controls[controls['mrna_assignment'].str.contains('neg')]

    return controls_neg

#returns statistical description of all negative controls
def controls_neg_stat_desc():
    global neg_contr_st_desc
    global neg_contr_clean

    controls_neg.set_index('ID_REF', inplace=True)
    neg_contr_clean = controls_neg[controls_neg.columns[11:58]]
    neg_contr_st_desc = neg_contr_clean.describe(percentiles=[.25, .5, .75, .95])

    return neg_contr_clean, neg_contr_st_desc


#returns statistical description of the positive controls
def controls_pos_stat_desc():
    global pos_contr_st_desc
    global pos_contr_clean

    controls_pos.set_index('ID_REF', inplace=True)
    pos_contr_clean = controls_pos[controls_pos.columns[11:58]]
    pos_contr_clean = pos_contr_clean.transpose()
    pos_contr_st_desc = pos_contr_clean.describe(percentiles=[.25, .5, .75, .95])

    return pos_contr_clean, pos_contr_st_desc

#extracts data for the main genes (filter out all non main genes)
def data_main ():
    global main_expr

    main_expr = data.loc[data['category'] == 'main']

    return main_expr

#returns only genes with positive expression lvl compared to control genes
def neg_sample_comparison ():

    global st_sig_genes
    global st_sig_genes_clean

    st_sig_genes = data.loc[:, ['ID_REF']]
    for (gsm, gsm_data) in neg_contr_st_desc.iteritems():
        local_main_expr = main_expr.loc[:, ['ID_REF', gsm]]
        st_sig_genes_local = local_main_expr.loc[local_main_expr[gsm] > neg_contr_st_desc.at['95%', gsm]]
        #st_sig_genes = pd.concat([st_sig_genes, st_sig_genes_local], axis=1)
        st_sig_genes= st_sig_genes.merge(st_sig_genes_local, on=['ID_REF'], how='outer')

    st_sig_genes_clean = st_sig_genes.dropna()

    return st_sig_genes, st_sig_genes_clean


#takes argument 'target' = every data that need to be assigned by ID_REF
##mozna dodac jeszcze zmienna by = zeby mozna bylo wybrac wg czego dopasowywac
def gene_assignment (target):
    global gene_data_assigned

    gene_data = data.loc[:, ['ID_REF', 'gene_assignment']]
    gene_data['gene_assignment'] = gene_data['gene_assignment'].values.astype('str')
    gene_names = []
    for gene_name in gene_data['gene_assignment']:
        if len(gene_name) > 6:
            gene_names.append(gene_name.rsplit('//')[1])
        else:
            gene_names.append(gene_name)
    gene_data['gene_name'] = gene_names
    gene_data['gene_name'] = gene_data['gene_name'].str.strip()
    gene_data = gene_data.drop('gene_assignment', axis='columns')
    gene_data_assigned = target.merge(gene_data, on=['ID_REF'], how='left')
    first_index = gene_data_assigned.pop('gene_name')
    gene_data_assigned.insert(1, 'gene_name', first_index)

    return gene_data_assigned

#returns on;y genes with negative expression
def neg_sample_comparison_neg_expr ():

    global st_sig_genes_neg_expr
    global st_sig_genes_clean_neg_expr

    st_sig_genes_neg_expr = data.loc[:, ['ID_REF']]
    for (gsm, gsm_data) in neg_contr_st_desc.iteritems():
        local_main_expr = main_expr.loc[:, ['ID_REF', gsm]]
        st_sig_genes_local = local_main_expr.loc[local_main_expr[gsm] < neg_contr_st_desc.at['95%', gsm]]
        #st_sig_genes = pd.concat([st_sig_genes, st_sig_genes_local], axis=1)
        st_sig_genes_neg_expr= st_sig_genes_neg_expr.merge(st_sig_genes_local, on=['ID_REF'], how='outer')

    st_sig_genes_clean_neg_expr = st_sig_genes_neg_expr.dropna()

    return st_sig_genes_neg_expr, st_sig_genes_clean_neg_expr


#takes all genes and assignes their epression levels to the 95% of control expression level
def determining_expr_lvl ():

    global data_expr_lvl

    data_expr_lvl = main_expr.loc[:, ['ID_REF']]
    for (gsm, gsm_data) in neg_contr_st_desc.iteritems():
        local_main_expr = main_expr.loc[:, ['ID_REF', gsm]]
        values = []
        for a in local_main_expr[gsm]:
            a -= neg_contr_st_desc.at['95%', gsm]
            values.append(a)
        local_main_expr[gsm] = values
        #st_sig_genes = pd.concat([st_sig_genes, st_sig_genes_local], axis=1)

        data_expr_lvl= data_expr_lvl.merge(local_main_expr, on=['ID_REF'], how='outer')

    return data_expr_lvl

def correlation (gene_name):
    global correlation_data

    coor_test_data = data_expr_lvl
    coor_test_data = coor_test_data.set_index('ID_REF').T
    gene_id = gene_data_assigned.loc[gene_data_assigned['gene_name'] == gene_name, 'ID_REF'].iloc[0]
    gene = coor_test_data[gene_id]
    pearson_improved = []
    pearson_p_value = []
    spearman_improved = []
    spearman_p_value = []
    kendall_improved = []
    kendall_p_value = []
    for a in coor_test_data:
        i = coor_test_data[a]
        pearson_improved.append(scipy.stats.pearsonr(i, gene)[0])
        pearson_p_value.append(scipy.stats.pearsonr(i, gene)[1])
        spearman_improved.append(scipy.stats.spearmanr(i, gene)[0])
        spearman_p_value.append(scipy.stats.spearmanr(i, gene)[1])
        kendall_improved.append(scipy.stats.kendalltau(i, gene)[0])
        kendall_p_value.append(scipy.stats.kendalltau(i, gene)[1])

    correlation_data = pd.DataFrame()
    correlation_data['ID_REF'] = data_expr_lvl['ID_REF']
    correlation_data['Pearson r'] = pearson_improved
    correlation_data['Pearson p-value'] = pearson_p_value
    correlation_data['Spearman rho'] = spearman_improved
    correlation_data['Spearman p-value'] = spearman_p_value
    correlation_data['Kendall tau'] = kendall_improved
    correlation_data['Kendall p-value'] = kendall_p_value
    gene_assignment_local = gene_data_assigned.loc[:, ['ID_REF', 'gene_name']]
    correlation_data = correlation_data.merge(gene_assignment_local, on='ID_REF', how='left')
    first_index = correlation_data.pop('gene_name')
    correlation_data.insert(1, 'gene_name', first_index)


    return correlation_data

def only_significant_corr(corr_data):
    global only_significant_correlation

    only_significant_correlation = corr_data.loc[corr_data['Pearson p-value'] <= 0.05]
    only_significant_correlation = only_significant_correlation.loc[only_significant_correlation['Spearman p-value'] <= 0.05]
    only_significant_correlation = only_significant_correlation.loc[only_significant_correlation['Kendall p-value'] <= 0.05]

def determine_correlation(corr):
    if corr > -0.3 and corr < 0.3:
        return 'negligible correlation'
    elif corr >= 0.3 and corr < 0.5:
        return 'low positive correlation'
    elif corr >= 0.5 and corr < 0.7:
        return 'moderate positive correlation'
    elif corr >= 0.7 and corr < 0.9:
        return 'high positive correlation'
    elif corr >= 0.9:
        return 'very high positive correlation'
    elif corr <= -0.3 and corr > -0.5:
        return 'low negative correlation'
    elif corr <= -0.5 and corr > -0.7:
        return 'moderate positive correlation'
    elif corr <= -0.7 and corr > -0.9:
        return 'high negative correlation'
    elif corr <= 0.9:
        return 'very high negative correlation'

only_significant_correlation['correlation size'] = only_significant_correlation['Pearson r'].apply(determine_correlation)





data_unpack()
data_controls()
data_controls_neg()
controls_neg_stat_desc()
data_main()
gene_assignment(main_expr)
determining_expr_lvl()
correlation('EGFR')
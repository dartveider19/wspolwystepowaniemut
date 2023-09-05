import GEOparse
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np


##data_unpack("GSE72217_family.soft.gz", "C:/Users/micle/Desktop/Python/dyplom")

#path = the folder where is your file
##name = the name of the file in that folder
###change backslashes to slashes

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

def data_controls ():
    global controls

    controls = data.loc[data['SPOT_ID'] == 'control']

    return controls

def data_controls_pos ():
    global controls_pos

    controls_pos = controls[controls['mrna_assignment'].str.contains('pos')]

    return controls_pos

def data_controls_neg():
    global controls_neg

    controls_neg = controls[controls['mrna_assignment'].str.contains('neg')]

    return controls_neg

def controls_neg_stat_desc():
    global neg_contr_st_desc
    global neg_contr_clean

    controls_neg.set_index('ID_REF', inplace=True)
    neg_contr_clean = controls_neg[controls_neg.columns[11:58]]
    neg_contr_st_desc = neg_contr_clean.describe(percentiles=[.25, .5, .75, .95])

    return neg_contr_clean, neg_contr_st_desc

def controls_pos_stat_desc():
    global pos_contr_st_desc
    global pos_contr_clean

    controls_pos.set_index('ID_REF', inplace=True)
    pos_contr_clean = controls_pos[controls_pos.columns[11:58]]
    pos_contr_clean = pos_contr_clean.transpose()
    pos_contr_st_desc = pos_contr_clean.describe(percentiles=[.25, .5, .75, .95])

    return pos_contr_clean, pos_contr_st_desc

def data_main ():
    global main_expr

    main_expr = data.loc[data['category'] == 'main']

    return main_expr

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



def gene_assignment ():
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
    gene_data_assigned = st_sig_genes_clean.merge(gene_data, on=['ID_REF'], how='left')
    first_index = gene_data_assigned.pop('gene_name')
    gene_data_assigned.insert(1, 'gene_name', first_index)

    return gene_data_assigned



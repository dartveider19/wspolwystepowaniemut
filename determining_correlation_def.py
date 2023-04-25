import GEOparse
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import numpy as np



#path = the folder where is your file
##name = the name of the file in that folder
###change backslashes to slashes

def data_unpack(name, path):
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
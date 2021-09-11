import GEOparse
import pandas as pd
import matplotlib.pyplot as plt


"""
To kod, ktory wyciaga 5 genow ref z okreslonych chromosomow, usunalbym to bo i tak mam kod, ktory wyciaga te wszystkie z hrt atlas
"""
gse = GEOparse.get_GEO(filepath="./GSE72217_family.soft.gz")

print ()
print ("GSM example: ")
for gsm_name in gse.gsms.items ():
    print ("Name: ", gsm_name[0])
    print ('Patient number: ', gsm_name[1].metadata.get('title'))
    for gpl_name in gse.gpls.items ():
        print ("Dataset: ", gpl_name)
        name = str(gpl_name[0])
        gpl = gse.gpls['GPL5175'].table
    gpl_egfr = gpl[gpl['seqname']=='chr7']
    gpl_egfr = gpl_egfr[gpl_egfr['RANGE_START'].astype(int) >= 55086576]
    gpl_egfr = gpl_egfr[gpl_egfr['RANGE_STOP'].astype(int) <= 55324279]
    gpl_egfr.sort_values(by=['ID'], inplace=True)
    gpl_egfr.reset_index(drop=True, inplace=True)
    
    sample_name = str(gsm_name[0])
    
    t = gse.gsms[sample_name].table
    out_egfr = gpl_egfr.set_index('ID').join(t.set_index('ID_REF'))
    out_egfr.reset_index(drop=False, inplace=True)
    out_egfr = out_egfr[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_egfr.dropna(subset=['VALUE'], inplace=True)
    
    '''find a way to forge dataframes in one''' #see below

    gpl_hprt1 = gpl[gpl['seqname']=='chrX']
    gpl_hprt1 = gpl_hprt1[gpl_hprt1['RANGE_START'].astype(int) >= 133594194]
    gpl_hprt1 = gpl_hprt1[gpl_hprt1['RANGE_STOP'].astype(int) <= 133669219]
    gpl_hprt1.sort_values(by=['ID'], inplace=True)
    gpl_hprt1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_hprt1 = gpl_hprt1.set_index('ID').join(t.set_index('ID_REF'))
    out_hprt1.reset_index(drop=False, inplace=True)
    out_hprt1 = out_hprt1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_hprt1.dropna(subset=['VALUE'], inplace=True)
    out_hprt1 = out_hprt1.append (out_egfr, ignore_index=True)



    gpl_nestin = gpl[gpl['seqname']=='chr1']
    gpl_nestin = gpl_nestin[gpl_nestin['RANGE_START'].astype(int) >= 156638593]
    gpl_nestin = gpl_nestin[gpl_nestin['RANGE_STOP'].astype(int) <= 156647203]
    gpl_nestin.sort_values(by=['ID'], inplace=True)
    gpl_nestin.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin = gpl_nestin.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin.reset_index(drop=False, inplace=True)
    out_nestin = out_nestin[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin.dropna(subset=['VALUE'], inplace=True)
    out_nestin = out_nestin.append (out_hprt1, ignore_index=True)
    



    gpl_b2m = gpl[gpl['seqname']=='chr15']
    gpl_b2m = gpl_b2m[gpl_b2m['RANGE_START'].astype(int) >= 45003685]
    gpl_b2m = gpl_b2m[gpl_b2m['RANGE_STOP'].astype(int) <= 45011055]
    gpl_b2m.sort_values(by=['ID'], inplace=True)
    gpl_b2m.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_b2m = gpl_b2m.set_index('ID').join(t.set_index('ID_REF'))
    out_b2m.reset_index(drop=False, inplace=True)
    out_b2m = out_b2m[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_b2m.dropna(subset=['VALUE'], inplace=True)
    out_b2m = out_b2m.append (out_nestin, ignore_index=True)
    


    gpl_actin = gpl[gpl['seqname'] == 'chr7']
    gpl_actin = gpl_actin[gpl_actin['RANGE_START'].astype(int) >= 5566907]
    gpl_actin = gpl_actin[gpl_actin['RANGE_STOP'].astype(int) <= 5570485]
    gpl_actin.sort_values(by=['ID'], inplace=True)
    gpl_actin.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_actin = gpl_actin.set_index('ID').join(t.set_index('ID_REF'))
    out_actin.reset_index(drop=False, inplace=True)
    out_actin = out_actin[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_actin.dropna(subset=['VALUE'], inplace=True)
    out_actin = out_actin.append (out_b2m, ignore_index=True)
    


    gpl_tata = gpl[gpl['seqname'] == 'chr6']
    gpl_tata = gpl_tata[gpl_tata['RANGE_START'].astype(int) >= 170863275]
    gpl_tata = gpl_tata[gpl_tata['RANGE_STOP'].astype(int) <= 170882020]
    gpl_tata.sort_values(by=['ID'], inplace=True)
    gpl_tata.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_tata = gpl_tata.set_index('ID').join(t.set_index('ID_REF'))
    out_tata.reset_index(drop=False, inplace=True)
    out_tata = out_tata[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_tata.dropna(subset=['VALUE'], inplace=True)
    out_tata = out_tata.append (out_actin, ignore_index=True)
    print(out_tata)
    out_tata.to_csv ('exp_results_' + str(gsm_name[1].metadata.get('title')) + '.csv')


    #convert columns to rows and rows to columns
    #first row  == empty name
    #add row row_name == patient ('title')
    #save as one file
    #done


    

#async def ref_gene_searcher ():
    #break

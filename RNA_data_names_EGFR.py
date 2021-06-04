import GEOparse
import pandas as pd
import matplotlib.pyplot as plt



gse = GEOparse.get_GEO(filepath="./GSE72218_family.soft.gz")

print ()
print ("GSM example: ")
for gsm_name in gse.gsms.items ():
    print ("Name: ", gsm_name)
    print (type(gsm_name))
    print (gsm_name[1].table)
    break
for gpl_name in gse.gpls.items ():
    print ("Dataset: ", gpl_name)
    print('gpl_name(type): ',type(gpl_name))
    print(gpl_name [1], type(gpl_name[1]))
    print (gpl_name[1].table)
    #gpl_name[1].table.to_csv ('Dataset_RNA_test1.csv')
    break
gpl = gse.gpls['GPL17586'].table
chromosome_7 = gpl[gpl['seqname']=='chr7']
#chromosome_9.to_csv('Chromosome_9_first_take_data.csv')
gpl_egfr = gpl[gpl['seqname']=='chr7']
gpl_egfr = gpl_egfr[gpl_egfr['start'].str.contains ("1|2|3|4|5|6|7|8|9|0", na=False)]
#gpl_egfr.to_csv ('EGFR_first_take.csv')
gpl_egfr = gpl_egfr[gpl_egfr['start'].astype(int) >= 53586678]
gpl_egfr = gpl_egfr[gpl_egfr['stop'].astype(int) <= 56779262]
gpl_egfr =  gpl_egfr[['ID', 'seqname','start', 'SPOT_ID']]
gpl_egfr.sort_values(by=['start'], inplace=True)
gpl_egfr.reset_index(drop=True, inplace=True)


t = gse.gsms['GSM1857830'].table
#out = chromosome_7.set_index('ID').join(t.set_index('ID_REF'))
#out.reset_index(drop=False, inplace=True)
#out = out[['ID', 'VALUE', 'seqname', 'SPOT_ID']]
#print(out)
#print (gpl_egfr)

select_id = t[t['ID_REF'].isin(gpl_egfr['ID'])]
select_id.reset_index(drop=True, inplace=True)
physical_position = gpl_egfr['start'].astype(int)

combine_by_key = select_id.assign(**{
    "Chromosome": gpl_egfr['seqname'], "Physical Position": physical_position
})
#combine_by_key.sort_values(by=['Physical Position'], inplace=True)
#combine_by_key = combine_by_key.assign(Classification=intron_classification)
print (combine_by_key)


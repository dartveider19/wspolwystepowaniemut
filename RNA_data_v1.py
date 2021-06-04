import GEOparse
import pandas as pd
import matplotlib.pyplot as plt



gse = GEOparse.get_GEO(filepath="./GSE72218_family.soft.gz")

print ()
print ("GSM example: ")
for gsm_name in gse.gsms.items ():
    print ("Name: ", gsm_name)
    print ("Metadata: ",)

    break
for key in gse.gsms.items ():
    print (key)
    print (type (key))
    print (len(key))
    print (type(key[0]))
    print (type (key[1]))
    print (key[1].show_metadata())
    #print( key[1]) #patient data inside
    #print (key[1].metadata)
    print (key[1].table)
    #key[1].table.plot()
    #plt.show()
    print (key[1].columns)
    key[1].table.to_csv ('output_data_RNA_test1.csv')
    break
 

#for key in gse.metadata.items ():
    #print (key)
    #print (type(key))
#print (gse.gpls['GPL16131'].columns)
#print (gse.gsms["GSM1857567"].columns)
    

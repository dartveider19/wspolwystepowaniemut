import GEOparse
import pandas



gse = GEOparse.get_GEO(geo="GSE72209", destdir="./")
print()
print("GSM example:")
gsm_name = 'GSM1857567'
for gsm_name, gsm in gse.gsms.items():
    print("Name: ", gsm_name)
    print("Metadata:",)
    for key, value in gsm.metadata.items():
        print(" - %s : %s" % (key, ", ".join(value)))
    print ("Table data:",)
    print (gsm.table.head())
    break

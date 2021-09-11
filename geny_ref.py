import GEOparse
import pandas as pd
import csv

"""
To jest kod, ktory wybiera geny referencyjne dla kazdego pacjenta i pozniej skleja wszystkie geny w jeden plik .csv
W kometarzach sa jakies mysli na temat jak to powinno dzialac tylkko ze nie umiem tego doprowadzic do takiej postaci, mozna to zrobic pozniej
"""

# gse = GEOparse.get_GEO (geoparse_link)
#geoparse_link = input('please enter link to GEOdataset of study, that you need')

gse = GEOparse.get_GEO(filepath="./GSE72217_family.soft.gz")

print ("GSM example: ")
for gsm_name in gse.gsms.items ():
    print ('Name:', gsm_name[0])
    print ('Patient:', gsm_name[1].metadata.get('title'))
    for gpl_name in gse.gpls.items():
        print ('Dataset:', gpl_name[0])
        gpl = gse.gpls[gpl_name[0]].table
        
    #input to find a gene in GPL
    #take ID from GPL
    #compare ID in GSM to find expression level
    '''gene_name = str(print (input ('Write the name of gene you are looking for ')))'''

    sample_name = str(gsm_name[0])
#PEA15
    gpl_nestin1_1 = gpl[gpl['seqname']=='chr1']
    gpl_nestin1_1 = gpl_nestin1_1[gpl_nestin1_1['RANGE_START'].astype(int) >= 160175125]
    gpl_nestin1_1 = gpl_nestin1_1[gpl_nestin1_1['RANGE_STOP'].astype(int) <= 160185154]
    gpl_nestin1_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin1_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin1_1 = gpl_nestin1_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin1_1.reset_index(drop=False, inplace=True)
    out_nestin1_1 = out_nestin1_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin1_1.dropna(subset=['VALUE'], inplace=True)
    #out_nestin1_1 = out_nestin1_1.append (out_nestin2, ignore_index=True)


#GNB
    gpl_nestin1_2 = gpl[gpl['seqname']=='chr1']
    gpl_nestin1_2 = gpl_nestin1_2[gpl_nestin1_2['RANGE_START'].astype(int) >= 1714732]
    gpl_nestin1_2 = gpl_nestin1_2[gpl_nestin1_2['RANGE_STOP'].astype(int) <= 1822637]
    gpl_nestin1_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin1_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin1_2 = gpl_nestin1_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin1_2.reset_index(drop=False, inplace=True)
    out_nestin1_2 = out_nestin1_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin1_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin1_2 = out_nestin1_1.append (out_nestin1_2, ignore_index=True)

#PARK7(?)
    gpl_nestin1_3 = gpl[gpl['seqname']=='chr1']
    gpl_nestin1_3 = gpl_nestin1_3[gpl_nestin1_3['RANGE_START'].astype(int) >= 8021733]
    gpl_nestin1_3 = gpl_nestin1_3[gpl_nestin1_3['RANGE_STOP'].astype(int) <= 8045545]
    gpl_nestin1_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin1_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin1_3 = gpl_nestin1_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin1_3.reset_index(drop=False, inplace=True)
    out_nestin1_3 = out_nestin1_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin1_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin1_3 = out_nestin1_2.append (out_nestin1_3, ignore_index=True)


#S100A6
    gpl_nestin1_4 = gpl[gpl['seqname']=='chr1']
    gpl_nestin1_4 = gpl_nestin1_4[gpl_nestin1_4['RANGE_START'].astype(int) >= 153507115]
    gpl_nestin1_4 = gpl_nestin1_4[gpl_nestin1_4['RANGE_STOP'].astype(int) <= 153509224]
    gpl_nestin1_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin1_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin1_4 = gpl_nestin1_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin1_4.reset_index(drop=False, inplace=True)
    out_nestin1_4 = out_nestin1_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin1_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin1_4 = out_nestin1_3.append (out_nestin1_4, ignore_index=True)


#SSR2(?)
    gpl_nestin1_5 = gpl[gpl['seqname']=='chr1']
    gpl_nestin1_5 = gpl_nestin1_5[gpl_nestin1_5['RANGE_START'].astype(int) >= 153507115]
    gpl_nestin1_5 = gpl_nestin1_5[gpl_nestin1_5['RANGE_STOP'].astype(int) <= 153509224]
    gpl_nestin1_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin1_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin1_5 = gpl_nestin1_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin1_5.reset_index(drop=False, inplace=True)
    out_nestin1_5 = out_nestin1_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin1_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin1_5 = out_nestin1_4.append (out_nestin1_5, ignore_index=True)

    
#NCL (nucleoline)
    gpl_nestin2_1 = gpl[gpl['seqname']=='chr2']
    gpl_nestin2_1 = gpl_nestin2_1[gpl_nestin2_1['RANGE_START'].astype(int) >= 232319469]
    gpl_nestin2_1 = gpl_nestin2_1[gpl_nestin2_1['RANGE_STOP'].astype(int) <= 232372334]
    gpl_nestin2_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin2_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin2_1 = gpl_nestin2_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin2_1.reset_index(drop=False, inplace=True)
    out_nestin2_1 = out_nestin2_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin2_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin2_1 = out_nestin1_5.append (out_nestin2_1, ignore_index=True)

#CAB39
    gpl_nestin2_2 = gpl[gpl['seqname']=='chr2']
    gpl_nestin2_2 = gpl_nestin2_2[gpl_nestin2_2['RANGE_START'].astype(int) >= 231577537]
    gpl_nestin2_2 = gpl_nestin2_2[gpl_nestin2_2['RANGE_STOP'].astype(int) <= 231685783]
    gpl_nestin2_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin2_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin2_2 = gpl_nestin2_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin2_2.reset_index(drop=False, inplace=True)
    out_nestin2_2 = out_nestin2_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin2_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin2_2 = out_nestin2_1.append (out_nestin2_2, ignore_index=True)


#AGFG1
    gpl_nestin2_3 = gpl[gpl['seqname']=='chr2']
    gpl_nestin2_3 = gpl_nestin2_3[gpl_nestin2_3['RANGE_START'].astype(int) >= 228336878]
    gpl_nestin2_3 = gpl_nestin2_3[gpl_nestin2_3['RANGE_STOP'].astype(int) <= 228425923]
    gpl_nestin2_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin2_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin2_3 = gpl_nestin2_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin2_3.reset_index(drop=False, inplace=True)
    out_nestin2_3 = out_nestin2_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin2_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin2_3 = out_nestin2_2.append (out_nestin2_3, ignore_index=True)


#MRPL44
    gpl_nestin2_4 = gpl[gpl['seqname']=='chr2']
    gpl_nestin2_4 = gpl_nestin2_4[gpl_nestin2_4['RANGE_START'].astype(int) >= 224822134]
    gpl_nestin2_4 = gpl_nestin2_4[gpl_nestin2_4['RANGE_STOP'].astype(int) <= 224832422]
    gpl_nestin2_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin2_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin2_4 = gpl_nestin2_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin2_4.reset_index(drop=False, inplace=True)
    out_nestin2_4 = out_nestin2_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin2_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin2_4 = out_nestin2_3.append (out_nestin2_4, ignore_index=True)


#DNAJB2 (?)
    gpl_nestin2_5 = gpl[gpl['seqname']=='chr2']
    gpl_nestin2_5 = gpl_nestin2_5[gpl_nestin2_5['RANGE_START'].astype(int) >= 220144061]
    gpl_nestin2_5 = gpl_nestin2_5[gpl_nestin2_5['RANGE_STOP'].astype(int) <= 220152056]
    gpl_nestin2_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin2_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin2_5 = gpl_nestin2_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin2_5.reset_index(drop=False, inplace=True)
    out_nestin2_5 = out_nestin2_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin2_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin2_5 = out_nestin2_4.append (out_nestin2_5, ignore_index=True)


#RPL35A
    gpl_nestin3_1 = gpl[gpl['seqname']=='chr3']
    gpl_nestin3_1 = gpl_nestin3_1[gpl_nestin3_1['RANGE_START'].astype(int) >= 197676879]
    gpl_nestin3_1 = gpl_nestin3_1[gpl_nestin3_1['RANGE_STOP'].astype(int) <= 197683466]
    gpl_nestin3_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin3_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin3_1 = gpl_nestin3_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin3_1.reset_index(drop=False, inplace=True)
    out_nestin3_1 = out_nestin3_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin3_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin3_1 = out_nestin2_5.append (out_nestin3_1, ignore_index=True)


#NCBP2
    gpl_nestin3_2 = gpl[gpl['seqname']=='chr3']
    gpl_nestin3_2 = gpl_nestin3_2[gpl_nestin3_2['RANGE_START'].astype(int) >= 196662280]
    gpl_nestin3_2 = gpl_nestin3_2[gpl_nestin3_2['RANGE_STOP'].astype(int) <= 196670366]
    gpl_nestin3_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin3_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin3_2 = gpl_nestin3_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin3_2.reset_index(drop=False, inplace=True)
    out_nestin3_2 = out_nestin3_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin3_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin3_2 = out_nestin3_1.append (out_nestin3_2, ignore_index=True)


#TMEM41A
    gpl_nestin3_3 = gpl[gpl['seqname']=='chr3']
    gpl_nestin3_3 = gpl_nestin3_3[gpl_nestin3_3['RANGE_START'].astype(int) >= 185197971]
    gpl_nestin3_3 = gpl_nestin3_3[gpl_nestin3_3['RANGE_STOP'].astype(int) <= 185216821]
    gpl_nestin3_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin3_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin3_3 = gpl_nestin3_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin3_3.reset_index(drop=False, inplace=True)
    out_nestin3_3 = out_nestin3_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin3_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin3_3 = out_nestin3_2.append (out_nestin3_3, ignore_index=True)


#MRPL47
    gpl_nestin3_4 = gpl[gpl['seqname']=='chr3']
    gpl_nestin3_4 = gpl_nestin3_4[gpl_nestin3_4['RANGE_START'].astype(int) >= 179301073]
    gpl_nestin3_4 = gpl_nestin3_4[gpl_nestin3_4['RANGE_STOP'].astype(int) <= 179322739]
    gpl_nestin3_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin3_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin3_4 = gpl_nestin3_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin3_4.reset_index(drop=False, inplace=True)
    out_nestin3_4 = out_nestin3_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin3_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin3_4 = out_nestin3_3.append (out_nestin3_4, ignore_index=True)


#KPNA1
    gpl_nestin3_5 = gpl[gpl['seqname']=='chr3']
    gpl_nestin3_5 = gpl_nestin3_5[gpl_nestin3_5['RANGE_START'].astype(int) >= 122135254]
    gpl_nestin3_5 = gpl_nestin3_5[gpl_nestin3_5['RANGE_STOP'].astype(int) <= 122233889]
    gpl_nestin3_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin3_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin3_5 = gpl_nestin3_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin3_5.reset_index(drop=False, inplace=True)
    out_nestin3_5 = out_nestin3_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin3_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin3_5 = out_nestin3_4.append (out_nestin3_5, ignore_index=True)


#FRG1
    gpl_nestin4_1 = gpl[gpl['seqname']=='chr4']
    gpl_nestin4_1 = gpl_nestin4_1[gpl_nestin4_1['RANGE_START'].astype(int) >= 190814874]
    gpl_nestin4_1 = gpl_nestin4_1[gpl_nestin4_1['RANGE_STOP'].astype(int) <= 190901312]
    gpl_nestin4_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin4_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin4_1 = gpl_nestin4_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin4_1.reset_index(drop=False, inplace=True)
    out_nestin4_1 = out_nestin4_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin4_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin4_1 = out_nestin3_5.append (out_nestin4_1, ignore_index=True)


    #TMEM184C
    gpl_nestin4_2 = gpl[gpl['seqname']=='chr4']
    gpl_nestin4_2 = gpl_nestin4_2[gpl_nestin4_2['RANGE_START'].astype(int) >= 148507714]
    gpl_nestin4_2 = gpl_nestin4_2[gpl_nestin4_2['RANGE_STOP'].astype(int) <= 148593195]
    gpl_nestin4_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin4_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin4_2 = gpl_nestin4_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin4_2.reset_index(drop=False, inplace=True)
    out_nestin4_2 = out_nestin4_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin4_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin4_2 = out_nestin4_1.append (out_nestin4_2, ignore_index=True)


    #INTS12
    gpl_nestin4_3 = gpl[gpl['seqname']=='chr4']
    gpl_nestin4_3 = gpl_nestin4_3[gpl_nestin4_3['RANGE_START'].astype(int) >= 106603183]
    gpl_nestin4_3 = gpl_nestin4_3[gpl_nestin4_3['RANGE_STOP'].astype(int) <= 106630432]
    gpl_nestin4_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin4_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin4_3 = gpl_nestin4_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin4_3.reset_index(drop=False, inplace=True)
    out_nestin4_3 = out_nestin4_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin4_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin4_3 = out_nestin4_2.append (out_nestin4_3, ignore_index=True)


    #TMEM165
    gpl_nestin4_4 = gpl[gpl['seqname']=='chr4']
    gpl_nestin4_4 = gpl_nestin4_4[gpl_nestin4_4['RANGE_START'].astype(int) >= 56257019]
    gpl_nestin4_4 = gpl_nestin4_4[gpl_nestin4_4['RANGE_STOP'].astype(int) <= 56319563]
    gpl_nestin4_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin4_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin4_4 = gpl_nestin4_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin4_4.reset_index(drop=False, inplace=True)
    out_nestin4_4 = out_nestin4_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin4_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin4_4 = out_nestin4_3.append (out_nestin4_4, ignore_index=True)


    #RPL9
    gpl_nestin4_5 = gpl[gpl['seqname']=='chr4']
    gpl_nestin4_5 = gpl_nestin4_5[gpl_nestin4_5['RANGE_START'].astype(int) >= 39454785]
    gpl_nestin4_5 = gpl_nestin4_5[gpl_nestin4_5['RANGE_STOP'].astype(int) <= 39460550]
    gpl_nestin4_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin4_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin4_5 = gpl_nestin4_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin4_5.reset_index(drop=False, inplace=True)
    out_nestin4_5 = out_nestin4_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin4_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin4_5 = out_nestin4_4.append (out_nestin4_5, ignore_index=True)


     #HNRNPAB
    gpl_nestin5_1 = gpl[gpl['seqname']=='chr5']
    gpl_nestin5_1 = gpl_nestin5_1[gpl_nestin5_1['RANGE_START'].astype(int) >= 177631518]
    gpl_nestin5_1 = gpl_nestin5_1[gpl_nestin5_1['RANGE_STOP'].astype(int) <= 177638158]
    gpl_nestin5_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin5_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin5_1 = gpl_nestin5_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin5_1.reset_index(drop=False, inplace=True)
    out_nestin5_1 = out_nestin5_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin5_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin5_1 = out_nestin4_5.append (out_nestin5_1, ignore_index=True)


    #LMAN2
    gpl_nestin5_2 = gpl[gpl['seqname']=='chr5']
    gpl_nestin5_2 = gpl_nestin5_2[gpl_nestin5_2['RANGE_START'].astype(int) >= 176743801]
    gpl_nestin5_2 = gpl_nestin5_2[gpl_nestin5_2['RANGE_STOP'].astype(int) <= 176778833]
    gpl_nestin5_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin5_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin5_2 = gpl_nestin5_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin5_2.reset_index(drop=False, inplace=True)
    out_nestin5_2 = out_nestin5_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin5_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin5_2 = out_nestin5_1.append (out_nestin5_2, ignore_index=True)


    #HIGD2A
    gpl_nestin5_3 = gpl[gpl['seqname']=='chr5']
    gpl_nestin5_3 = gpl_nestin5_3[gpl_nestin5_3['RANGE_START'].astype(int) >= 175815821]
    gpl_nestin5_3 = gpl_nestin5_3[gpl_nestin5_3['RANGE_STOP'].astype(int) <= 175816738]
    gpl_nestin5_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin5_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin5_3 = gpl_nestin5_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin5_3.reset_index(drop=False, inplace=True)
    out_nestin5_3 = out_nestin5_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin5_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin5_3 = out_nestin5_2.append (out_nestin5_3, ignore_index=True)


    #SLU7
    gpl_nestin5_4 = gpl[gpl['seqname']=='chr5']
    gpl_nestin5_4 = gpl_nestin5_4[gpl_nestin5_4['RANGE_START'].astype(int) >= 159828648]
    gpl_nestin5_4 = gpl_nestin5_4[gpl_nestin5_4['RANGE_STOP'].astype(int) <= 159846381]
    gpl_nestin5_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin5_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin5_4 = gpl_nestin5_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin5_4.reset_index(drop=False, inplace=True)
    out_nestin5_4 = out_nestin5_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin5_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin5_4 = out_nestin5_3.append (out_nestin5_4, ignore_index=True)


    #MRPL22
    gpl_nestin5_5 = gpl[gpl['seqname']=='chr5']
    gpl_nestin5_5 = gpl_nestin5_5[gpl_nestin5_5['RANGE_START'].astype(int) >= 154320628]
    gpl_nestin5_5 = gpl_nestin5_5[gpl_nestin5_5['RANGE_STOP'].astype(int) <= 154346528]
    gpl_nestin5_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin5_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin5_5 = gpl_nestin5_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin5_5.reset_index(drop=False, inplace=True)
    out_nestin5_5 = out_nestin5_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin5_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin5_5 = out_nestin5_4.append (out_nestin5_5, ignore_index=True)


    #PSMB1
    gpl_nestin6_1 = gpl[gpl['seqname']=='chr6']
    gpl_nestin6_1 = gpl_nestin6_1[gpl_nestin6_1['RANGE_START'].astype(int) >= 170764272]
    gpl_nestin6_1 = gpl_nestin6_1[gpl_nestin6_1['RANGE_STOP'].astype(int) <= 170864476]
    gpl_nestin6_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin6_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin6_1 = gpl_nestin6_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin6_1.reset_index(drop=False, inplace=True)
    out_nestin6_1 = out_nestin6_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin6_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin6_1 = out_nestin5_5.append (out_nestin6_1, ignore_index=True)


    #MRPL18
    gpl_nestin6_2 = gpl[gpl['seqname']=='chr6']
    gpl_nestin6_2 = gpl_nestin6_2[gpl_nestin6_2['RANGE_START'].astype(int) >= 160209889]
    gpl_nestin6_2 = gpl_nestin6_2[gpl_nestin6_2['RANGE_STOP'].astype(int) <= 160219740]
    gpl_nestin6_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin6_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin6_2 = gpl_nestin6_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin6_2.reset_index(drop=False, inplace=True)
    out_nestin6_2 = out_nestin6_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin6_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin6_2 = out_nestin6_1.append (out_nestin6_2, ignore_index=True)


    #DYNLT1
    gpl_nestin6_3 = gpl[gpl['seqname']=='chr6']
    gpl_nestin6_3 = gpl_nestin6_3[gpl_nestin6_3['RANGE_START'].astype(int) >= 159057514]
    gpl_nestin6_3 = gpl_nestin6_3[gpl_nestin6_3['RANGE_STOP'].astype(int) <= 159086613]
    gpl_nestin6_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin6_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin6_3 = gpl_nestin6_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin6_3.reset_index(drop=False, inplace=True)
    out_nestin6_3 = out_nestin6_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin6_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin6_3 = out_nestin6_2.append (out_nestin6_3, ignore_index=True)


    #IFNGR1
    gpl_nestin6_4 = gpl[gpl['seqname']=='chr6']
    gpl_nestin6_4 = gpl_nestin6_4[gpl_nestin6_4['RANGE_START'].astype(int) >= 137518624]
    gpl_nestin6_4 = gpl_nestin6_4[gpl_nestin6_4['RANGE_STOP'].astype(int) <= 137540586]
    gpl_nestin6_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin6_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin6_4 = gpl_nestin6_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin6_4.reset_index(drop=False, inplace=True)
    out_nestin6_4 = out_nestin6_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin6_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin6_4 = out_nestin6_3.append (out_nestin6_4, ignore_index=True)


    #RPS12
    gpl_nestin6_5 = gpl[gpl['seqname']=='chr6']
    gpl_nestin6_5 = gpl_nestin6_5[gpl_nestin6_5['RANGE_START'].astype(int) >= 39421198]
    gpl_nestin6_5 = gpl_nestin6_5[gpl_nestin6_5['RANGE_STOP'].astype(int) <= 39423782]
    gpl_nestin6_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin6_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin6_5 = gpl_nestin6_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin6_5.reset_index(drop=False, inplace=True)
    out_nestin6_5 = out_nestin6_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin6_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin6_5 = out_nestin6_4.append (out_nestin6_5, ignore_index=True)


    #ABCF2
    gpl_nestin7_1 = gpl[gpl['seqname']=='chr7']
    gpl_nestin7_1 = gpl_nestin7_1[gpl_nestin7_1['RANGE_START'].astype(int) >= 150888698]
    gpl_nestin7_1 = gpl_nestin7_1[gpl_nestin7_1['RANGE_STOP'].astype(int) <= 150924646]
    gpl_nestin7_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin7_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin7_1 = gpl_nestin7_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin7_1.reset_index(drop=False, inplace=True)
    out_nestin7_1 = out_nestin7_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin7_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin7_1 = out_nestin6_5.append (out_nestin7_1, ignore_index=True)


    #NDUFB2
    gpl_nestin7_2 = gpl[gpl['seqname']=='chr7']
    gpl_nestin7_2 = gpl_nestin7_2[gpl_nestin7_2['RANGE_START'].astype(int) >= 140396481]
    gpl_nestin7_2 = gpl_nestin7_2[gpl_nestin7_2['RANGE_STOP'].astype(int) <= 140426118]
    gpl_nestin7_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin7_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin7_2 = gpl_nestin7_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin7_2.reset_index(drop=False, inplace=True)
    out_nestin7_2 = out_nestin7_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin7_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin7_2 = out_nestin7_1.append (out_nestin7_2, ignore_index=True)


    #POLR2J
    gpl_nestin7_3 = gpl[gpl['seqname']=='chr7']
    gpl_nestin7_3 = gpl_nestin7_3[gpl_nestin7_3['RANGE_START'].astype(int) >= 102113575]
    gpl_nestin7_3 = gpl_nestin7_3[gpl_nestin7_3['RANGE_STOP'].astype(int) <= 102119454]
    gpl_nestin7_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin7_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin7_3 = gpl_nestin7_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin7_3.reset_index(drop=False, inplace=True)
    out_nestin7_3 = out_nestin7_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin7_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin7_3 = out_nestin7_2.append (out_nestin7_3, ignore_index=True)


    #PLOD3
    gpl_nestin7_4 = gpl[gpl['seqname']=='chr7']
    gpl_nestin7_4 = gpl_nestin7_4[gpl_nestin7_4['RANGE_START'].astype(int) >= 100849285]
    gpl_nestin7_4 = gpl_nestin7_4[gpl_nestin7_4['RANGE_STOP'].astype(int) <= 100861142]
    gpl_nestin7_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin7_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin7_4 = gpl_nestin7_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin7_4.reset_index(drop=False, inplace=True)
    out_nestin7_4 = out_nestin7_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin7_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin7_4 = out_nestin7_3.append (out_nestin7_4, ignore_index=True)


    #GNB2
    gpl_nestin7_5 = gpl[gpl['seqname']=='chr7']
    gpl_nestin7_5 = gpl_nestin7_5[gpl_nestin7_5['RANGE_START'].astype(int) >= 100271375]
    gpl_nestin7_5 = gpl_nestin7_5[gpl_nestin7_5['RANGE_STOP'].astype(int) <= 100276792]
    gpl_nestin7_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin7_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin7_5 = gpl_nestin7_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin7_5.reset_index(drop=False, inplace=True)
    out_nestin7_5 = out_nestin7_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin7_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin7_5 = out_nestin7_4.append (out_nestin7_5, ignore_index=True)



    #RPL8
    gpl_nestin8_1 = gpl[gpl['seqname']=='chr8']
    gpl_nestin8_1 = gpl_nestin8_1[gpl_nestin8_1['RANGE_START'].astype(int) >= 146015190]
    gpl_nestin8_1 = gpl_nestin8_1[gpl_nestin8_1['RANGE_STOP'].astype(int) <= 146017811]
    gpl_nestin8_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin8_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin8_1 = gpl_nestin8_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin8_1.reset_index(drop=False, inplace=True)
    out_nestin8_1 = out_nestin8_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin8_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin8_1 = out_nestin7_5.append (out_nestin8_1, ignore_index=True)


    #CYC1
    gpl_nestin8_2 = gpl[gpl['seqname']=='chr8']
    gpl_nestin8_2 = gpl_nestin8_2[gpl_nestin8_2['RANGE_START'].astype(int) >= 145149540]
    gpl_nestin8_2 = gpl_nestin8_2[gpl_nestin8_2['RANGE_STOP'].astype(int) <= 145152398]
    gpl_nestin8_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin8_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin8_2 = gpl_nestin8_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin8_2.reset_index(drop=False, inplace=True)
    out_nestin8_2 = out_nestin8_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin8_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin8_2 = out_nestin8_1.append (out_nestin8_2, ignore_index=True)


    #AZIN1
    gpl_nestin8_3 = gpl[gpl['seqname']=='chr8']
    gpl_nestin8_3 = gpl_nestin8_3[gpl_nestin8_3['RANGE_START'].astype(int) >= 103838549]
    gpl_nestin8_3 = gpl_nestin8_3[gpl_nestin8_3['RANGE_STOP'].astype(int) <= 103884791]
    gpl_nestin8_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin8_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin8_3 = gpl_nestin8_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin8_3.reset_index(drop=False, inplace=True)
    out_nestin8_3 = out_nestin8_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin8_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin8_3 = out_nestin8_2.append (out_nestin8_3, ignore_index=True)


    #POLR2K
    gpl_nestin8_4 = gpl[gpl['seqname']=='chr8']
    gpl_nestin8_4 = gpl_nestin8_4[gpl_nestin8_4['RANGE_START'].astype(int) >= 101162832]
    gpl_nestin8_4 = gpl_nestin8_4[gpl_nestin8_4['RANGE_STOP'].astype(int) <= 101166217]
    gpl_nestin8_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin8_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin8_4 = gpl_nestin8_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin8_4.reset_index(drop=False, inplace=True)
    out_nestin8_4 = out_nestin8_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin8_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin8_4 = out_nestin8_3.append (out_nestin8_4, ignore_index=True)


    #MRPL41
    gpl_nestin9_1 = gpl[gpl['seqname']=='chr9']
    gpl_nestin9_1 = gpl_nestin9_1[gpl_nestin9_1['RANGE_START'].astype(int) >= 140445671]
    gpl_nestin9_1 = gpl_nestin9_1[gpl_nestin9_1['RANGE_STOP'].astype(int) <= 140446992]
    gpl_nestin9_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin9_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin9_1 = gpl_nestin9_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin9_1.reset_index(drop=False, inplace=True)
    out_nestin9_1 = out_nestin9_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin9_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin9_1 = out_nestin8_4.append (out_nestin9_1, ignore_index=True)


    #GOLGA2
    gpl_nestin9_2 = gpl[gpl['seqname']=='chr9']
    gpl_nestin9_2 = gpl_nestin9_2[gpl_nestin9_2['RANGE_START'].astype(int) >= 131018115]
    gpl_nestin9_2 = gpl_nestin9_2[gpl_nestin9_2['RANGE_STOP'].astype(int) <= 131039180]
    gpl_nestin9_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin9_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin9_2 = gpl_nestin9_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin9_2.reset_index(drop=False, inplace=True)
    out_nestin9_2 = out_nestin9_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin9_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin9_2 = out_nestin9_1.append (out_nestin9_2, ignore_index=True)


    #NDUFA8
    gpl_nestin9_3 = gpl[gpl['seqname']=='chr9']
    gpl_nestin9_3 = gpl_nestin9_3[gpl_nestin9_3['RANGE_START'].astype(int) >= 124906365]
    gpl_nestin9_3 = gpl_nestin9_3[gpl_nestin9_3['RANGE_STOP'].astype(int) <= 124922098]
    gpl_nestin9_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin9_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin9_3 = gpl_nestin9_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin9_3.reset_index(drop=False, inplace=True)
    out_nestin9_3 = out_nestin9_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin9_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin9_3 = out_nestin9_2.append (out_nestin9_3, ignore_index=True)


    #POLR1E
    gpl_nestin9_4 = gpl[gpl['seqname']=='chr9']
    gpl_nestin9_4 = gpl_nestin9_4[gpl_nestin9_4['RANGE_START'].astype(int) >= 37485945]
    gpl_nestin9_4 = gpl_nestin9_4[gpl_nestin9_4['RANGE_STOP'].astype(int) <= 37510379]
    gpl_nestin9_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin9_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin9_4 = gpl_nestin9_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin9_4.reset_index(drop=False, inplace=True)
    out_nestin9_4 = out_nestin9_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin9_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin9_4 = out_nestin9_3.append (out_nestin9_4, ignore_index=True)


    #CDKN2A
    gpl_nestin9_5 = gpl[gpl['seqname']=='chr9']
    gpl_nestin9_5 = gpl_nestin9_5[gpl_nestin9_5['RANGE_START'].astype(int) >= 21939019]
    gpl_nestin9_5 = gpl_nestin9_5[gpl_nestin9_5['RANGE_STOP'].astype(int) <= 21995248]
    gpl_nestin9_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin9_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin9_5 = gpl_nestin9_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin9_5.reset_index(drop=False, inplace=True)
    out_nestin9_5 = out_nestin9_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin9_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin9_5 = out_nestin9_4.append (out_nestin9_5, ignore_index=True)


    #DCTN3
    gpl_nestin9_6 = gpl[gpl['seqname']=='chr9']
    gpl_nestin9_6 = gpl_nestin9_6[gpl_nestin9_6['RANGE_START'].astype(int) >= 34613417]
    gpl_nestin9_6 = gpl_nestin9_6[gpl_nestin9_6['RANGE_STOP'].astype(int) <= 34624091]
    gpl_nestin9_6.sort_values(by=['ID'], inplace=True)
    gpl_nestin9_6.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin9_6 = gpl_nestin9_6.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin9_6.reset_index(drop=False, inplace=True)
    out_nestin9_6 = out_nestin9_6[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin9_6.dropna(subset=['VALUE'], inplace=True)
    out_nestin9_6 = out_nestin9_5.append (out_nestin9_6, ignore_index=True)


    #PSAP(?)
    gpl_nestin10_1 = gpl[gpl['seqname']=='chr10']
    gpl_nestin10_1 = gpl_nestin10_1[gpl_nestin10_1['RANGE_START'].astype(int) >= 73576066]
    gpl_nestin10_1 = gpl_nestin10_1[gpl_nestin10_1['RANGE_STOP'].astype(int) <= 73611106]
    gpl_nestin10_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin10_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin10_1 = gpl_nestin10_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin10_1.reset_index(drop=False, inplace=True)
    out_nestin10_1 = out_nestin10_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin10_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin10_1 = out_nestin9_6.append (out_nestin10_1, ignore_index=True)


    #PPA1(?)
    gpl_nestin10_2 = gpl[gpl['seqname']=='chr10']
    gpl_nestin10_2 = gpl_nestin10_2[gpl_nestin10_2['RANGE_START'].astype(int) >= 71962598]
    gpl_nestin10_2 = gpl_nestin10_2[gpl_nestin10_2['RANGE_STOP'].astype(int) <= 71993647]
    gpl_nestin10_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin10_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin10_2 = gpl_nestin10_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin10_2.reset_index(drop=False, inplace=True)
    out_nestin10_2 = out_nestin10_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin10_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin10_2 = out_nestin10_1.append (out_nestin10_2, ignore_index=True)


    #ACTR1A(?)
    gpl_nestin10_3 = gpl[gpl['seqname']=='chr10']
    gpl_nestin10_3 = gpl_nestin10_3[gpl_nestin10_3['RANGE_START'].astype(int) >= 104239030]
    gpl_nestin10_3 = gpl_nestin10_3[gpl_nestin10_3['RANGE_STOP'].astype(int) <= 104274458]
    gpl_nestin10_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin10_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin10_3 = gpl_nestin10_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin10_3.reset_index(drop=False, inplace=True)
    out_nestin10_3 = out_nestin10_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin10_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin10_3 = out_nestin10_2.append (out_nestin10_3, ignore_index=True)


    #GSTO1(?)
    gpl_nestin10_4 = gpl[gpl['seqname']=='chr10']
    gpl_nestin10_4 = gpl_nestin10_4[gpl_nestin10_4['RANGE_START'].astype(int) >= 105995133]
    gpl_nestin10_4 = gpl_nestin10_4[gpl_nestin10_4['RANGE_STOP'].astype(int) <= 106027899]
    gpl_nestin10_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin10_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin10_4 = gpl_nestin10_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin10_4.reset_index(drop=False, inplace=True)
    out_nestin10_4 = out_nestin10_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin10_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin10_4 = out_nestin10_3.append (out_nestin10_4, ignore_index=True)


    #MRPS16
    gpl_nestin10_5 = gpl[gpl['seqname']=='chr10']
    gpl_nestin10_5 = gpl_nestin10_5[gpl_nestin10_5['RANGE_START'].astype(int) >= 75010510]
    gpl_nestin10_5 = gpl_nestin10_5[gpl_nestin10_5['RANGE_STOP'].astype(int) <= 75012440]
    gpl_nestin10_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin10_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin10_5 = gpl_nestin10_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin10_5.reset_index(drop=False, inplace=True)
    out_nestin10_5 = out_nestin10_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin10_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin10_5 = out_nestin10_4.append (out_nestin10_5, ignore_index=True)


    #COX8A
    gpl_nestin11_1 = gpl[gpl['seqname']=='chr11']
    gpl_nestin11_1 = gpl_nestin11_1[gpl_nestin11_1['RANGE_START'].astype(int) >= 63742055]
    gpl_nestin11_1 = gpl_nestin11_1[gpl_nestin11_1['RANGE_STOP'].astype(int) <= 63743986]
    gpl_nestin11_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin11_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin11_1 = gpl_nestin11_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin11_1.reset_index(drop=False, inplace=True)
    out_nestin11_1 = out_nestin11_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin11_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin11_1 = out_nestin10_5.append (out_nestin11_1, ignore_index=True)


    #AP2A2
    gpl_nestin11_2 = gpl[gpl['seqname']=='chr11']
    gpl_nestin11_2 = gpl_nestin11_2[gpl_nestin11_2['RANGE_START'].astype(int) >= 922387]
    gpl_nestin11_2 = gpl_nestin11_2[gpl_nestin11_2['RANGE_STOP'].astype(int) <= 1012226]
    gpl_nestin11_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin11_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin11_2 = gpl_nestin11_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin11_2.reset_index(drop=False, inplace=True)
    out_nestin11_2 = out_nestin11_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin11_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin11_2 = out_nestin11_1.append (out_nestin11_2, ignore_index=True)


    #SF3B2
    gpl_nestin11_3 = gpl[gpl['seqname']=='chr11']
    gpl_nestin11_3 = gpl_nestin11_3[gpl_nestin11_3['RANGE_START'].astype(int) >= 65816923]
    gpl_nestin11_3 = gpl_nestin11_3[gpl_nestin11_3['RANGE_STOP'].astype(int) <= 65836767]
    gpl_nestin11_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin11_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin11_3 = gpl_nestin11_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin11_3.reset_index(drop=False, inplace=True)
    out_nestin11_3 = out_nestin11_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin11_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin11_3 = out_nestin11_2.append (out_nestin11_3, ignore_index=True)


    #EIF4G2(?)
    gpl_nestin11_4 = gpl[gpl['seqname']=='chr11']
    gpl_nestin11_4 = gpl_nestin11_4[gpl_nestin11_4['RANGE_START'].astype(int) >= 10530998]
    gpl_nestin11_4 = gpl_nestin11_4[gpl_nestin11_4['RANGE_STOP'].astype(int) <= 10832170]
    gpl_nestin11_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin11_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin11_4 = gpl_nestin11_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin11_4.reset_index(drop=False, inplace=True)
    out_nestin11_4 = out_nestin11_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin11_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin11_4 = out_nestin11_3.append (out_nestin11_4, ignore_index=True)


    #CCT2(?)
    gpl_nestin12_1 = gpl[gpl['seqname']=='chr12']
    gpl_nestin12_1 = gpl_nestin12_1[gpl_nestin12_1['RANGE_START'].astype(int) >= 69979241]
    gpl_nestin12_1 = gpl_nestin12_1[gpl_nestin12_1['RANGE_STOP'].astype(int) <= 69995345]
    gpl_nestin12_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin12_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin12_1 = gpl_nestin12_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin12_1.reset_index(drop=False, inplace=True)
    out_nestin12_1 = out_nestin12_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin12_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin12_1 = out_nestin11_4.append (out_nestin12_1, ignore_index=True)


    #SNRNP35
    gpl_nestin12_2 = gpl[gpl['seqname']=='chr12']
    gpl_nestin12_2 = gpl_nestin12_2[gpl_nestin12_2['RANGE_START'].astype(int) >= 123942201]
    gpl_nestin12_2 = gpl_nestin12_2[gpl_nestin12_2['RANGE_STOP'].astype(int) <= 123956909]
    gpl_nestin12_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin12_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin12_2 = gpl_nestin12_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin12_2.reset_index(drop=False, inplace=True)
    out_nestin12_2 = out_nestin12_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin12_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin12_2 = out_nestin12_1.append (out_nestin12_2, ignore_index=True)


    #NDUFA12
    gpl_nestin12_3 = gpl[gpl['seqname']=='chr12']
    gpl_nestin12_3 = gpl_nestin12_3[gpl_nestin12_3['RANGE_START'].astype(int) >= 95290851]
    gpl_nestin12_3 = gpl_nestin12_3[gpl_nestin12_3['RANGE_STOP'].astype(int) <= 95397511]
    gpl_nestin12_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin12_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin12_3 = gpl_nestin12_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin12_3.reset_index(drop=False, inplace=True)
    out_nestin12_3 = out_nestin12_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin12_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin12_3 = out_nestin12_2.append (out_nestin12_3, ignore_index=True)


    #XPOT
    gpl_nestin12_4 = gpl[gpl['seqname']=='chr12']
    gpl_nestin12_4 = gpl_nestin12_4[gpl_nestin12_4['RANGE_START'].astype(int) >= 64798148]
    gpl_nestin12_4 = gpl_nestin12_4[gpl_nestin12_4['RANGE_STOP'].astype(int) <= 64844895]
    gpl_nestin12_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin12_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin12_4 = gpl_nestin12_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin12_4.reset_index(drop=False, inplace=True)
    out_nestin12_4 = out_nestin12_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin12_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin12_4 = out_nestin12_3.append (out_nestin12_4, ignore_index=True)


    #DCTN2
    gpl_nestin12_5 = gpl[gpl['seqname']=='chr12']
    gpl_nestin12_5 = gpl_nestin12_5[gpl_nestin12_5['RANGE_START'].astype(int) >= 57923338]
    gpl_nestin12_5 = gpl_nestin12_5[gpl_nestin12_5['RANGE_STOP'].astype(int) <= 57941011]
    gpl_nestin12_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin12_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin12_5 = gpl_nestin12_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin12_5.reset_index(drop=False, inplace=True)
    out_nestin12_5 = out_nestin12_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin12_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin12_5 = out_nestin12_4.append (out_nestin12_5, ignore_index=True)


    #IPO5
    gpl_nestin13_1 = gpl[gpl['seqname']=='chr13']
    gpl_nestin13_1 = gpl_nestin13_1[gpl_nestin13_1['RANGE_START'].astype(int) >= 98572093]
    gpl_nestin13_1 = gpl_nestin13_1[gpl_nestin13_1['RANGE_STOP'].astype(int) <= 98688931]
    gpl_nestin13_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin13_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin13_1 = gpl_nestin13_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin13_1.reset_index(drop=False, inplace=True)
    out_nestin13_1 = out_nestin13_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin13_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin13_1 = out_nestin12_5.append (out_nestin13_1, ignore_index=True)


    #GTF2F2
    gpl_nestin13_2 = gpl[gpl['seqname']=='chr13']
    gpl_nestin13_2 = gpl_nestin13_2[gpl_nestin13_2['RANGE_START'].astype(int) >= 45694617]
    gpl_nestin13_2 = gpl_nestin13_2[gpl_nestin13_2['RANGE_STOP'].astype(int) <= 45858234]
    gpl_nestin13_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin13_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin13_2 = gpl_nestin13_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin13_2.reset_index(drop=False, inplace=True)
    out_nestin13_2 = out_nestin13_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin13_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin13_2 = out_nestin13_1.append (out_nestin13_2, ignore_index=True)


    #EXOSC8
    gpl_nestin13_3 = gpl[gpl['seqname']=='chr13']
    gpl_nestin13_3 = gpl_nestin13_3[gpl_nestin13_3['RANGE_START'].astype(int) >= 37572953]
    gpl_nestin13_3 = gpl_nestin13_3[gpl_nestin13_3['RANGE_STOP'].astype(int) <= 37583750]
    gpl_nestin13_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin13_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin13_3 = gpl_nestin13_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin13_3.reset_index(drop=False, inplace=True)
    out_nestin13_3 = out_nestin13_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin13_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin13_3 = out_nestin13_2.append (out_nestin13_3, ignore_index=True)


    #POLR1D
    gpl_nestin13_4 = gpl[gpl['seqname']=='chr13']
    gpl_nestin13_4 = gpl_nestin13_4[gpl_nestin13_4['RANGE_START'].astype(int) >= 28194903]
    gpl_nestin13_4 = gpl_nestin13_4[gpl_nestin13_4['RANGE_STOP'].astype(int) <= 28271155]
    gpl_nestin13_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin13_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin13_4 = gpl_nestin13_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin13_4.reset_index(drop=False, inplace=True)
    out_nestin13_4 = out_nestin13_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin13_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin13_4 = out_nestin13_3.append (out_nestin13_4, ignore_index=True)


    #PAPOLA
    gpl_nestin14_1 = gpl[gpl['seqname']=='chr14']
    gpl_nestin14_1 = gpl_nestin14_1[gpl_nestin14_1['RANGE_START'].astype(int) >= 96968285]
    gpl_nestin14_1 = gpl_nestin14_1[gpl_nestin14_1['RANGE_STOP'].astype(int) <= 97034874]
    gpl_nestin14_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin14_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin14_1 = gpl_nestin14_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin14_1.reset_index(drop=False, inplace=True)
    out_nestin14_1 = out_nestin14_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin14_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin14_1 = out_nestin13_4.append (out_nestin14_1, ignore_index=True)


    #NDUFB1
    gpl_nestin14_2 = gpl[gpl['seqname']=='chr14']
    gpl_nestin14_2 = gpl_nestin14_2[gpl_nestin14_2['RANGE_START'].astype(int) >= 92582472]
    gpl_nestin14_2 = gpl_nestin14_2[gpl_nestin14_2['RANGE_STOP'].astype(int) <= 92588153]
    gpl_nestin14_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin14_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin14_2 = gpl_nestin14_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin14_2.reset_index(drop=False, inplace=True)
    out_nestin14_2 = out_nestin14_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin14_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin14_2 = out_nestin14_1.append (out_nestin14_2, ignore_index=True)


    #PSMA3
    gpl_nestin14_3 = gpl[gpl['seqname']=='chr14']
    gpl_nestin14_3 = gpl_nestin14_3[gpl_nestin14_3['RANGE_START'].astype(int) >= 58517788]
    gpl_nestin14_3 = gpl_nestin14_3[gpl_nestin14_3['RANGE_STOP'].astype(int) <= 58740629]
    gpl_nestin14_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin14_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin14_3 = gpl_nestin14_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin14_3.reset_index(drop=False, inplace=True)
    out_nestin14_3 = out_nestin14_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin14_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin14_3 = out_nestin14_2.append (out_nestin14_3, ignore_index=True)


    #TINF2
    gpl_nestin14_4 = gpl[gpl['seqname']=='chr14']
    gpl_nestin14_4 = gpl_nestin14_4[gpl_nestin14_4['RANGE_START'].astype(int) >= 24706625]
    gpl_nestin14_4 = gpl_nestin14_4[gpl_nestin14_4['RANGE_STOP'].astype(int) <= 24712110]
    gpl_nestin14_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin14_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin14_4 = gpl_nestin14_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin14_4.reset_index(drop=False, inplace=True)
    out_nestin14_4 = out_nestin14_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin14_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin14_4 = out_nestin14_3.append (out_nestin14_4, ignore_index=True)


    #PSME2
    gpl_nestin14_5 = gpl[gpl['seqname']=='chr14']
    gpl_nestin14_5 = gpl_nestin14_5[gpl_nestin14_5['RANGE_START'].astype(int) >= 24612580]
    gpl_nestin14_5 = gpl_nestin14_5[gpl_nestin14_5['RANGE_STOP'].astype(int) <= 24616532]
    gpl_nestin14_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin14_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin14_5 = gpl_nestin14_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin14_5.reset_index(drop=False, inplace=True)
    out_nestin14_5 = out_nestin14_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin14_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin14_5 = out_nestin14_4.append (out_nestin14_5, ignore_index=True)


    #MRPL46
    gpl_nestin15_1 = gpl[gpl['seqname']=='chr15']
    gpl_nestin15_1 = gpl_nestin15_1[gpl_nestin15_1['RANGE_START'].astype(int) >= 88994011]
    gpl_nestin15_1 = gpl_nestin15_1[gpl_nestin15_1['RANGE_STOP'].astype(int) <= 89010596]
    gpl_nestin15_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin15_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin15_1 = gpl_nestin15_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin15_1.reset_index(drop=False, inplace=True)
    out_nestin15_1 = out_nestin15_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin15_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin15_1 = out_nestin14_5.append (out_nestin15_1, ignore_index=True)


    #RPL4
    gpl_nestin15_2 = gpl[gpl['seqname']=='chr15']
    gpl_nestin15_2 = gpl_nestin15_2[gpl_nestin15_2['RANGE_START'].astype(int) >= 185135240]
    gpl_nestin15_2 = gpl_nestin15_2[gpl_nestin15_2['RANGE_STOP'].astype(int) <= 185136616]
    gpl_nestin15_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin15_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin15_2 = gpl_nestin15_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin15_2.reset_index(drop=False, inplace=True)
    out_nestin15_2 = out_nestin15_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin15_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin15_2 = out_nestin15_1.append (out_nestin15_2, ignore_index=True)


    #RPL13
    gpl_nestin16_1 = gpl[gpl['seqname']=='chr16']
    gpl_nestin16_1 = gpl_nestin16_1[gpl_nestin16_1['RANGE_START'].astype(int) >= 89627076]
    gpl_nestin16_1 = gpl_nestin16_1[gpl_nestin16_1['RANGE_STOP'].astype(int) <= 89631049]
    gpl_nestin16_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin16_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin16_1 = gpl_nestin16_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin16_1.reset_index(drop=False, inplace=True)
    out_nestin16_1 = out_nestin16_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin16_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin16_1 = out_nestin15_2.append (out_nestin16_1, ignore_index=True)


    #MAP1LC3B
    gpl_nestin16_2 = gpl[gpl['seqname']=='chr16']
    gpl_nestin16_2 = gpl_nestin16_2[gpl_nestin16_2['RANGE_START'].astype(int) >= 87425416]
    gpl_nestin16_2 = gpl_nestin16_2[gpl_nestin16_2['RANGE_STOP'].astype(int) <= 87439363]
    gpl_nestin16_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin16_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin16_2 = gpl_nestin16_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin16_2.reset_index(drop=False, inplace=True)
    out_nestin16_2 = out_nestin16_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin16_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin16_2 = out_nestin16_1.append (out_nestin16_2, ignore_index=True)


    #COX4I1
    gpl_nestin16_3 = gpl[gpl['seqname']=='chr16']
    gpl_nestin16_3 = gpl_nestin16_3[gpl_nestin16_3['RANGE_START'].astype(int) >= 85822600]
    gpl_nestin16_3 = gpl_nestin16_3[gpl_nestin16_3['RANGE_STOP'].astype(int) <= 85840592]
    gpl_nestin16_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin16_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin16_3 = gpl_nestin16_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin16_3.reset_index(drop=False, inplace=True)
    out_nestin16_3 = out_nestin16_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin16_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin16_3 = out_nestin16_2.append (out_nestin16_3, ignore_index=True)


    #AARS
    gpl_nestin16_4 = gpl[gpl['seqname']=='chr16']
    gpl_nestin16_4 = gpl_nestin16_4[gpl_nestin16_4['RANGE_START'].astype(int) >= 70286212]
    gpl_nestin16_4 = gpl_nestin16_4[gpl_nestin16_4['RANGE_STOP'].astype(int) <= 70323421]
    gpl_nestin16_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin16_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin16_4 = gpl_nestin16_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin16_4.reset_index(drop=False, inplace=True)
    out_nestin16_4 = out_nestin16_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin16_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin16_4 = out_nestin16_3.append (out_nestin16_4, ignore_index=True)


    #DYNC1LI2
    gpl_nestin16_5 = gpl[gpl['seqname']=='chr16']
    gpl_nestin16_5 = gpl_nestin16_5[gpl_nestin16_5['RANGE_START'].astype(int) >= 66754803]
    gpl_nestin16_5 = gpl_nestin16_5[gpl_nestin16_5['RANGE_STOP'].astype(int) <= 66785711]
    gpl_nestin16_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin16_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin16_5 = gpl_nestin16_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin16_5.reset_index(drop=False, inplace=True)
    out_nestin16_5 = out_nestin16_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin16_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin16_5 = out_nestin16_4.append (out_nestin16_5, ignore_index=True)


    #ACTG1
    gpl_nestin17_1 = gpl[gpl['seqname']=='chr17']
    gpl_nestin17_1 = gpl_nestin17_1[gpl_nestin17_1['RANGE_START'].astype(int) >= 79476861]
    gpl_nestin17_1 = gpl_nestin17_1[gpl_nestin17_1['RANGE_STOP'].astype(int) <= 79491202]
    gpl_nestin17_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin17_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin17_1 = gpl_nestin17_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin17_1.reset_index(drop=False, inplace=True)
    out_nestin17_1 = out_nestin17_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin17_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin17_1 = out_nestin16_5.append (out_nestin17_1, ignore_index=True)


    #MRPL38
    gpl_nestin17_2 = gpl[gpl['seqname']=='chr17']
    gpl_nestin17_2 = gpl_nestin17_2[gpl_nestin17_2['RANGE_START'].astype(int) >= 73894736]
    gpl_nestin17_2 = gpl_nestin17_2[gpl_nestin17_2['RANGE_STOP'].astype(int) <= 73905654]
    gpl_nestin17_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin17_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin17_2 = gpl_nestin17_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin17_2.reset_index(drop=False, inplace=True)
    out_nestin17_2 = out_nestin17_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin17_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin17_2 = out_nestin17_1.append (out_nestin17_2, ignore_index=True)


    #MRPS7
    gpl_nestin17_3 = gpl[gpl['seqname']=='chr17']
    gpl_nestin17_3 = gpl_nestin17_3[gpl_nestin17_3['RANGE_START'].astype(int) >= 73257779]
    gpl_nestin17_3 = gpl_nestin17_3[gpl_nestin17_3['RANGE_STOP'].astype(int) <= 73263031]
    gpl_nestin17_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin17_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin17_3 = gpl_nestin17_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin17_3.reset_index(drop=False, inplace=True)
    out_nestin17_3 = out_nestin17_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin17_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin17_3 = out_nestin17_2.append (out_nestin17_3, ignore_index=True)


    #DYNLL2
    gpl_nestin17_4 = gpl[gpl['seqname']=='chr17']
    gpl_nestin17_4 = gpl_nestin17_4[gpl_nestin17_4['RANGE_START'].astype(int) >= 56124176]
    gpl_nestin17_4 = gpl_nestin17_4[gpl_nestin17_4['RANGE_STOP'].astype(int) <= 56192845]
    gpl_nestin17_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin17_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin17_4 = gpl_nestin17_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin17_4.reset_index(drop=False, inplace=True)
    out_nestin17_4 = out_nestin17_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin17_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin17_4 = out_nestin17_3.append (out_nestin17_4, ignore_index=True)


    #PHB
    gpl_nestin17_5 = gpl[gpl['seqname']=='chr17']
    gpl_nestin17_5 = gpl_nestin17_5[gpl_nestin17_5['RANGE_START'].astype(int) >= 47481220]
    gpl_nestin17_5 = gpl_nestin17_5[gpl_nestin17_5['RANGE_STOP'].astype(int) <= 47492243]
    gpl_nestin17_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin17_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin17_5 = gpl_nestin17_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin17_5.reset_index(drop=False, inplace=True)
    out_nestin17_5 = out_nestin17_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin17_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin17_5 = out_nestin17_4.append (out_nestin17_5, ignore_index=True)


    #RPL27A
    gpl_nestin17_6 = gpl[gpl['seqname']=='chr17']
    gpl_nestin17_6 = gpl_nestin17_6[gpl_nestin17_6['RANGE_START'].astype(int) >= 48445252]
    gpl_nestin17_6 = gpl_nestin17_6[gpl_nestin17_6['RANGE_STOP'].astype(int) <= 48450562]
    gpl_nestin17_6.sort_values(by=['ID'], inplace=True)
    gpl_nestin17_6.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin17_6 = gpl_nestin17_6.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin17_6.reset_index(drop=False, inplace=True)
    out_nestin17_6 = out_nestin17_6[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin17_6.dropna(subset=['VALUE'], inplace=True)
    out_nestin17_6 = out_nestin17_5.append (out_nestin17_6, ignore_index=True)


    #LMAN1
    gpl_nestin18_1 = gpl[gpl['seqname']=='chr18']
    gpl_nestin18_1 = gpl_nestin18_1[gpl_nestin18_1['RANGE_START'].astype(int) >= 56994862]
    gpl_nestin18_1 = gpl_nestin18_1[gpl_nestin18_1['RANGE_STOP'].astype(int) <= 57044113]
    gpl_nestin18_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin18_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin18_1 = gpl_nestin18_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin18_1.reset_index(drop=False, inplace=True)
    out_nestin18_1 = out_nestin18_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin18_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin18_1 = out_nestin17_6.append (out_nestin18_1, ignore_index=True)


    #USP14
    gpl_nestin18_2 = gpl[gpl['seqname']=='chr18']
    gpl_nestin18_2 = gpl_nestin18_2[gpl_nestin18_2['RANGE_START'].astype(int) >= 158498]
    gpl_nestin18_2 = gpl_nestin18_2[gpl_nestin18_2['RANGE_STOP'].astype(int) <= 214716]
    gpl_nestin18_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin18_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin18_2 = gpl_nestin18_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin18_2.reset_index(drop=False, inplace=True)
    out_nestin18_2 = out_nestin18_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin18_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin18_2 = out_nestin18_1.append (out_nestin18_2, ignore_index=True)


    #PSMG2
    gpl_nestin18_3 = gpl[gpl['seqname']=='chr18']
    gpl_nestin18_3 = gpl_nestin18_3[gpl_nestin18_3['RANGE_START'].astype(int) >= 12658062]
    gpl_nestin18_3 = gpl_nestin18_3[gpl_nestin18_3['RANGE_STOP'].astype(int) <= 12750532]
    gpl_nestin18_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin18_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin18_3 = gpl_nestin18_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin18_3.reset_index(drop=False, inplace=True)
    out_nestin18_3 = out_nestin18_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin18_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin18_3 = out_nestin18_2.append (out_nestin18_3, ignore_index=True)


    #EPN1
    gpl_nestin19_1 = gpl[gpl['seqname']=='chr19']
    gpl_nestin19_1 = gpl_nestin19_1[gpl_nestin19_1['RANGE_START'].astype(int) >= 56186592]
    gpl_nestin19_1 = gpl_nestin19_1[gpl_nestin19_1['RANGE_STOP'].astype(int) <= 56221204]
    gpl_nestin19_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin19_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin19_1 = gpl_nestin19_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin19_1.reset_index(drop=False, inplace=True)
    out_nestin19_1 = out_nestin19_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin19_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin19_1 = out_nestin18_3.append (out_nestin19_1, ignore_index=True)


    #RPL28
    gpl_nestin19_2 = gpl[gpl['seqname']=='chr19']
    gpl_nestin19_2 = gpl_nestin19_2[gpl_nestin19_2['RANGE_START'].astype(int) >= 55895814]
    gpl_nestin19_2 = gpl_nestin19_2[gpl_nestin19_2['RANGE_STOP'].astype(int) <= 55914690]
    gpl_nestin19_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin19_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin19_2 = gpl_nestin19_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin19_2.reset_index(drop=False, inplace=True)
    out_nestin19_2 = out_nestin19_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin19_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin19_2 = out_nestin19_1.append (out_nestin19_2, ignore_index=True)


    #ACTN4
    gpl_nestin19_3 = gpl[gpl['seqname']=='chr19']
    gpl_nestin19_3 = gpl_nestin19_3[gpl_nestin19_3['RANGE_START'].astype(int) >= 39138290]
    gpl_nestin19_3 = gpl_nestin19_3[gpl_nestin19_3['RANGE_STOP'].astype(int) <= 39221163]
    gpl_nestin19_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin19_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin19_3 = gpl_nestin19_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin19_3.reset_index(drop=False, inplace=True)
    out_nestin19_3 = out_nestin19_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin19_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin19_3 = out_nestin19_2.append (out_nestin19_3, ignore_index=True)


    #TBCB
    gpl_nestin19_4 = gpl[gpl['seqname']=='chr19']
    gpl_nestin19_4 = gpl_nestin19_4[gpl_nestin19_4['RANGE_START'].astype(int) >= 36605438]
    gpl_nestin19_4 = gpl_nestin19_4[gpl_nestin19_4['RANGE_STOP'].astype(int) <= 36616825]
    gpl_nestin19_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin19_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin19_4 = gpl_nestin19_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin19_4.reset_index(drop=False, inplace=True)
    out_nestin19_4 = out_nestin19_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin19_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin19_4 = out_nestin19_3.append (out_nestin19_4, ignore_index=True)


    #SIRT2
    gpl_nestin19_5 = gpl[gpl['seqname']=='chr19']
    gpl_nestin19_5 = gpl_nestin19_5[gpl_nestin19_5['RANGE_START'].astype(int) >= 39369211]
    gpl_nestin19_5 = gpl_nestin19_5[gpl_nestin19_5['RANGE_STOP'].astype(int) <= 39390361]
    gpl_nestin19_5.sort_values(by=['ID'], inplace=True)
    gpl_nestin19_5.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin19_5 = gpl_nestin19_5.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin19_5.reset_index(drop=False, inplace=True)
    out_nestin19_5 = out_nestin19_5[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin19_5.dropna(subset=['VALUE'], inplace=True)
    out_nestin19_5 = out_nestin19_4.append (out_nestin19_5, ignore_index=True)


    #POLR2I
    gpl_nestin19_6 = gpl[gpl['seqname']=='chr19']
    gpl_nestin19_6 = gpl_nestin19_6[gpl_nestin19_6['RANGE_START'].astype(int) >= 36604621]
    gpl_nestin19_6 = gpl_nestin19_6[gpl_nestin19_6['RANGE_STOP'].astype(int) <= 36606238]
    gpl_nestin19_6.sort_values(by=['ID'], inplace=True)
    gpl_nestin19_6.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin19_6 = gpl_nestin19_6.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin19_6.reset_index(drop=False, inplace=True)
    out_nestin19_6 = out_nestin19_6[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin19_6.dropna(subset=['VALUE'], inplace=True)
    out_nestin19_6 = out_nestin19_5.append (out_nestin19_6, ignore_index=True)


    #RPS21
    gpl_nestin20_1 = gpl[gpl['seqname']=='chr20']
    gpl_nestin20_1 = gpl_nestin20_1[gpl_nestin20_1['RANGE_START'].astype(int) >= 60951896]
    gpl_nestin20_1 = gpl_nestin20_1[gpl_nestin20_1['RANGE_STOP'].astype(int) <= 60963561]
    gpl_nestin20_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin20_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin20_1 = gpl_nestin20_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin20_1.reset_index(drop=False, inplace=True)
    out_nestin20_1 = out_nestin20_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin20_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin20_1 = out_nestin19_6.append (out_nestin20_1, ignore_index=True)


    #TOP1
    gpl_nestin20_2 = gpl[gpl['seqname']=='chr20']
    gpl_nestin20_2 = gpl_nestin20_2[gpl_nestin20_2['RANGE_START'].astype(int) >= 39657477]
    gpl_nestin20_2 = gpl_nestin20_2[gpl_nestin20_2['RANGE_STOP'].astype(int) <= 39753121]
    gpl_nestin20_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin20_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin20_2 = gpl_nestin20_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin20_2.reset_index(drop=False, inplace=True)
    out_nestin20_2 = out_nestin20_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin20_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin20_2 = out_nestin20_1.append (out_nestin20_2, ignore_index=True)


    #MRPS26
    gpl_nestin20_3 = gpl[gpl['seqname']=='chr20']
    gpl_nestin20_3 = gpl_nestin20_3[gpl_nestin20_3['RANGE_START'].astype(int) >= 3026591]
    gpl_nestin20_3 = gpl_nestin20_3[gpl_nestin20_3['RANGE_STOP'].astype(int) <= 3028889]
    gpl_nestin20_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin20_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin20_3 = gpl_nestin20_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin20_3.reset_index(drop=False, inplace=True)
    out_nestin20_3 = out_nestin20_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin20_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin20_3 = out_nestin20_2.append (out_nestin20_3, ignore_index=True)


    #MRPS6
    gpl_nestin21_1 = gpl[gpl['seqname']=='chr21']
    gpl_nestin21_1 = gpl_nestin21_1[gpl_nestin21_1['RANGE_START'].astype(int) >= 35322421]
    gpl_nestin21_1 = gpl_nestin21_1[gpl_nestin21_1['RANGE_STOP'].astype(int) <= 35516060]
    gpl_nestin21_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin21_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin21_1 = gpl_nestin21_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin21_1.reset_index(drop=False, inplace=True)
    out_nestin21_1 = out_nestin21_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin21_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin21_1 = out_nestin20_3.append (out_nestin21_1, ignore_index=True)


    #MRPL39
    gpl_nestin21_2 = gpl[gpl['seqname']=='chr21']
    gpl_nestin21_2 = gpl_nestin21_2[gpl_nestin21_2['RANGE_START'].astype(int) >= 26870083]
    gpl_nestin21_2 = gpl_nestin21_2[gpl_nestin21_2['RANGE_STOP'].astype(int) <= 26980158]
    gpl_nestin21_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin21_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin21_2 = gpl_nestin21_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin21_2.reset_index(drop=False, inplace=True)
    out_nestin21_2 = out_nestin21_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin21_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin21_2 = out_nestin21_1.append (out_nestin21_2, ignore_index=True)


    #CERK
    gpl_nestin22_1 = gpl[gpl['seqname']=='chr22']
    gpl_nestin22_1 = gpl_nestin22_1[gpl_nestin22_1['RANGE_START'].astype(int) >= 47075679]
    gpl_nestin22_1 = gpl_nestin22_1[gpl_nestin22_1['RANGE_STOP'].astype(int) <= 47158685]
    gpl_nestin22_1.sort_values(by=['ID'], inplace=True)
    gpl_nestin22_1.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin22_1 = gpl_nestin22_1.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin22_1.reset_index(drop=False, inplace=True)
    out_nestin22_1 = out_nestin22_1[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin22_1.dropna(subset=['VALUE'], inplace=True)
    out_nestin22_1 = out_nestin21_2.append (out_nestin22_1, ignore_index=True)


    #PMM1
    gpl_nestin22_2 = gpl[gpl['seqname']=='chr22']
    gpl_nestin22_2 = gpl_nestin22_2[gpl_nestin22_2['RANGE_START'].astype(int) >= 41972934]
    gpl_nestin22_2 = gpl_nestin22_2[gpl_nestin22_2['RANGE_STOP'].astype(int) <= 41986451]
    gpl_nestin22_2.sort_values(by=['ID'], inplace=True)
    gpl_nestin22_2.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin22_2 = gpl_nestin22_2.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin22_2.reset_index(drop=False, inplace=True)
    out_nestin22_2 = out_nestin22_2[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin22_2.dropna(subset=['VALUE'], inplace=True)
    out_nestin22_2 = out_nestin22_1.append (out_nestin22_2, ignore_index=True)


    #TOMM22
    gpl_nestin22_3 = gpl[gpl['seqname']=='chr22']
    gpl_nestin22_3 = gpl_nestin22_3[gpl_nestin22_3['RANGE_START'].astype(int) >= 39076559]
    gpl_nestin22_3 = gpl_nestin22_3[gpl_nestin22_3['RANGE_STOP'].astype(int) <= 39080763]
    gpl_nestin22_3.sort_values(by=['ID'], inplace=True)
    gpl_nestin22_3.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin22_3 = gpl_nestin22_3.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin22_3.reset_index(drop=False, inplace=True)
    out_nestin22_3 = out_nestin22_3[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin22_3.dropna(subset=['VALUE'], inplace=True)
    out_nestin22_3 = out_nestin22_2.append (out_nestin22_3, ignore_index=True)


    #POLR2F
    gpl_nestin22_4 = gpl[gpl['seqname']=='chr22']
    gpl_nestin22_4 = gpl_nestin22_4[gpl_nestin22_4['RANGE_START'].astype(int) >= 38349684]
    gpl_nestin22_4 = gpl_nestin22_4[gpl_nestin22_4['RANGE_STOP'].astype(int) <= 38422669]
    gpl_nestin22_4.sort_values(by=['ID'], inplace=True)
    gpl_nestin22_4.reset_index(drop=True, inplace=True)
    t = gse.gsms[sample_name].table
    out_nestin22_4 = gpl_nestin22_4.set_index('ID').join(t.set_index('ID_REF'))
    out_nestin22_4.reset_index(drop=False, inplace=True)
    out_nestin22_4 = out_nestin22_4[['ID', 'VALUE', 'seqname', 'RANGE_START', 'RANGE_STOP','gene_assignment']]
    out_nestin22_4.dropna(subset=['VALUE'], inplace=True)
    out_nestin22_4 = out_nestin22_3.append (out_nestin22_4, ignore_index=True)

    out_nestin22_4.to_csv (str(gsm_name[1].metadata.get('title')) + '_referencja_v0.csv')
    print ('EXPORTED')


    
    

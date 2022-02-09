import pandas as pd
import sys

tab_merged=pd.read_table(sys.argv[1], header=None)
tab_merged=tab_merged.rename(columns={0:"chromosome" , 1:"start" , 2:"end" , 3:"samples"})
tab_merged['SI'] = tab_merged.apply(lambda row : len(set(row["samples"].split(","))), axis = 1)

tab_merged= tab_merged.drop("samples",1)
#tab_merged.to_csv("/home/stefano/Documents/IFO/folgiero/data/multiinter/SI_tab.bed", sep="\t", index=None)
tab_merged.to_csv(sys.argv[2], sep="\t", index=None, header=None)
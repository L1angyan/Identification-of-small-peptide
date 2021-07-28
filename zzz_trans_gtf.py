def trans(i):
    i = str(i)
    j='transcript_id "peptide_maize'+i+'"; '+'gene_id "peptide_maize'+i+'"'
    return(j)

def add_stop_codon(df):
    strand = df.iloc[0,6]
    length=len(df)
    if strand=="+":
        df.iloc[length-1,4] = df.iloc[length-1,4]+3
    else:
        df.iloc[length-1,3] = df.iloc[length-1,3]-3
    return(df)

import pandas as pd
import csv
gtf = pd.read_table("zzz_merge.gtf",sep="\t",header=None)
del(gtf[8])
gtf.iloc[:,8] = gtf[9].map(trans)
#number -> annotation
x = gtf.groupby(9)
#divide lines into gene models
gtf = x.apply(add_stop_codon)
#positive strand:end_position + 3
#negative strand:start_position - 3

gtf.to_csv("zzz_peptide_ok.gtf",sep="\t",header=None,index=None,quoting=csv.QUOTE_NONE)

import pandas as pd
import glob

res_file = glob.glob("res_merge/*txt")
res = pd.DataFrame()
for i in res_file:
    df = pd.read_table(i,sep="\t")
    df["length"] = df.ORF_ID.str.split("_").str.get(3).astype("int")
    df["tissue"] = i.split(".")[2]
    res = pd.concat([res,df],axis="index")

# filtering
res = res[(res.pval_combined<0.05) & (res.length<=100)]

res = res.loc[:,
    ['ORF_ID','ORF_type','transcript_type','gene_id','gene_name','gene_type',
     'chrom','strand','pval_combined','AAseq','length','tissue']].drop_duplicates().sort_values("ORF_ID")
res.to_csv("1.merge_res.csv",header=True,index=False)

res = res.groupby(
    ['ORF_ID','ORF_type','transcript_type','gene_id','gene_name','gene_type',
     'chrom','strand','AAseq','length'])['tissue'].agg('/'.join).reset_index()
res.to_csv("1.merge_res1.csv",header=True,index=False)

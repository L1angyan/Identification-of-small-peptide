# 千万不要再傻乎乎地在pandas里面写for循环，记住用map、apply、applymap三大函数！！！！！！


import pandas as pd
import csv

txt = pd.read_table("zzz_all.txt",sep="\t")
txt["number"] = txt.index
txt = txt.drop_duplicates(subset=['ORF_ID', 'strand', 'ORF_length', 'ORF_gstart', 'ORF_gstop', 'AAseq'],keep="first")
txt = txt[txt["ORF_length"]<=300]
#给txt文件加上编号，并去重复行
txt.to_csv("zzz_merge.txt",sep="\t",header=True,index=False)

gtf = pd.read_table("zzz_all.gtf",sep="\t",header=None)
gtf["number"] = 0
gtf.iloc[gtf[gtf[2]=="ORF"].index,9] = 1
gtf["Number"] = gtf["number"].cumsum()
gtf["Number"] = gtf["Number"] - 1
#number = -1
#for i in gtf.index:
#    if gtf.iloc[i,2] == "ORF":number +=1
#    gtf.iloc[i,9] = number
#给gtf文件加上编号
f = lambda e,l : (e in l)
overlap = gtf["Number"].apply(f,args=(txt["number"],))
#gtf的number列的行号是否在txt的number列中，输出True or False到overlap列表里面
gtf = gtf[overlap]
#gtf = gtf.iloc[:,0:9]
gtf.drop("number",axis=1,inplace=True)
gtf.to_csv("zzz_merge.gtf",sep="\t",header=False,index=False,quoting=csv.QUOTE_NONE)

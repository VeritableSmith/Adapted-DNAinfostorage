import glob, os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
sns.set_color_codes("bright")

from pylab import rcParams
from re import sub

# Local Imports:
from Classes.helloworld_ONT_decoding  import Decode_ONT



## ONT_Decoding:
data = pd.read_csv("H01-H12/Data/H01-H12_dataTable_ONT.csv",sep=",", header=0,names=["template_ID","match","template","strandC","strandR","time","read"])
data["time"] = pd.to_datetime(data["time"])

# start is defined as the time stamp of the first read

data["hr_from_start"] = (data["time"].subtract(min(data["time"])))/np.timedelta64(1, 'h')
data["min_from_start"] = (data["time"].subtract(min(data["time"])))/np.timedelta64(1, 'm')
data["sec_from_start"] = (data["time"].subtract(min(data["time"])))/np.timedelta64(1, 's')
data.fillna('',inplace=True)

data["possiblehit"] = data.apply(Decode_ONT.findpossiblehits,axis=1)

strands = [data.loc[data["template_ID"]=="H01"].shape[0],
           data.loc[data["template_ID"]=="H02"].shape[0],
           data.loc[data["template_ID"]=="H03"].shape[0],
           data.loc[data["template_ID"]=="H04"].shape[0],
           data.loc[data["template_ID"]=="H05"].shape[0],
           data.loc[data["template_ID"]=="H06"].shape[0],
           data.loc[data["template_ID"]=="H07"].shape[0],
           data.loc[data["template_ID"]=="H08"].shape[0],
           data.loc[data["template_ID"]=="H09"].shape[0],
           data.loc[data["template_ID"]=="H10"].shape[0],
           data.loc[data["template_ID"]=="H11"].shape[0],
           data.loc[data["template_ID"]=="H12"].shape[0]]
print(strands)

# set the number of strands per template to be 500
minsample=500

H01 = data.loc[data["template_ID"]=="H01"]
H02 = data.loc[data["template_ID"]=="H02"]
H03 = data.loc[data["template_ID"]=="H03"]
H04 = data.loc[data["template_ID"]=="H04"]
H05 = data.loc[data["template_ID"]=="H05"]
H06 = data.loc[data["template_ID"]=="H06"]
H07 = data.loc[data["template_ID"]=="H07"]
H08 = data.loc[data["template_ID"]=="H08"]
H09 = data.loc[data["template_ID"]=="H09"]
H10 = data.loc[data["template_ID"]=="H10"]
H11 = data.loc[data["template_ID"]=="H11"]
H12 = data.loc[data["template_ID"]=="H12"]

data_downsampled=pd.concat([H01.sample(n=minsample),H02.sample(n=minsample),H03.sample(n=minsample),H04.sample(n=minsample),H05.sample(n=minsample),
          H06.sample(n=minsample),H07.sample(n=minsample),H08.sample(n=minsample),H09.sample(n=minsample),H10.sample(n=minsample),
          H11.sample(n=minsample),H12.sample(n=minsample)],ignore_index=True)

# Some cool graph or something
plt.rcParams['figure.figsize'] = 10,10
sns.set(style="white")

colorscheme = ["#ececec"]
g=sns.stripplot(data = data_downsampled, x="hr_from_start",y="template_ID", alpha = 1,palette = colorscheme,jitter = 0.35,size=4)

data_possibles = data_downsampled.loc[data_downsampled["possiblehit"]==1]
colorscheme = ["#acacac"]
g=sns.stripplot(data = data_possibles, x="hr_from_start",y="template_ID", alpha = 1, palette = colorscheme,jitter = 0.35,size=4)

data_hitsonly = data_downsampled.loc[data_downsampled["match"]==1]
colorscheme = ["black"]
g= sns.stripplot(data = data_hitsonly, x="hr_from_start",y="template_ID", alpha = 0.8,palette = colorscheme, jitter = 0.35,size=4)
g.set(ylabel='template ID')
g.set(xlabel='sequencing time (hrs)')
sns.despine(offset=10,trim=True)


hrs = list(range(0,50,2))
hrs.pop(0)
print(hrs[1:])


iters = 50
# iters = 10000 #trials used in publication - long run time! may want to parallelize

H01_decoding = Decode_ONT.decodingprobs(H01,hrs[1:],iters)
H02_decoding = Decode_ONT.decodingprobs(H02,hrs[1:],iters)
H03_decoding = Decode_ONT.decodingprobs(H03,hrs[1:],iters)
H04_decoding = Decode_ONT.decodingprobs(H04,hrs[1:],iters)
H05_decoding = Decode_ONT.decodingprobs(H05,hrs[1:],iters)
H06_decoding = Decode_ONT.decodingprobs(H06,hrs[1:],iters)
H07_decoding = Decode_ONT.decodingprobs(H07,hrs[1:],iters)
H08_decoding = Decode_ONT.decodingprobs(H08,hrs[1:],iters)
H09_decoding = Decode_ONT.decodingprobs(H09,hrs[1:],iters)
H10_decoding = Decode_ONT.decodingprobs(H10,hrs[1:],iters)
H11_decoding = Decode_ONT.decodingprobs(H11,hrs[1:],iters)
H12_decoding = Decode_ONT.decodingprobs(H12,hrs[1:],iters)


pd.DataFrame(list(zip(hrs, H01_decoding, H02_decoding, H03_decoding, H04_decoding,
                     H05_decoding, H06_decoding, H07_decoding, H08_decoding,
                     H09_decoding, H10_decoding, H11_decoding, H12_decoding,)),
              columns=['time (hrs)','H01','H02', 'H03','H04','H05','H06', 'H07','H08',
                      'H09','H10', 'H11','H12'])


# plot the data for previously calculated results
# Set 48 hours as 1.0 of sequencing time

sns.set(font_scale=1.5)
sns.set_style("white")


strand_decoding = {'H01':[24,48],'H02':[8,48],'H03':[12,48],'H04':[6,48],'H05':[6,48],'H06':[12,48],
                   'H07':[4,48],'H08':[6,48],'H09':[4,48],'H10':[4,48],'H11':[4,48],'H12':[8,48]}

strand_decoding_df = pd.DataFrame.from_dict(strand_decoding)
strand_decoding_df = strand_decoding_df/48
strand_decoding_df["method"] = pd.Series(['ONT','MiSeq'])
strand_decoding_df_melted = pd.DataFrame.melt(strand_decoding_df,id_vars='method')


fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize=(7,3.5))
strand_decoding_df_melted_ONT = strand_decoding_df_melted.loc[strand_decoding_df_melted["method"]=='ONT']
strand_decoding_df_melted_MiSeq = strand_decoding_df_melted.loc[strand_decoding_df_melted["method"]=='MiSeq']

g = sns.pointplot(x="value",y='variable',data =strand_decoding_df_melted_ONT,markers=['D'],color="black",ax=ax,join=False,scale = 1)
g = sns.pointplot(x="value",y='variable',data =strand_decoding_df_melted_MiSeq,markers=['o'],color="black",ax=ax,join=False,scale = 1)

g.set(xlabel='Fraction of sequencing time')
g.set(ylabel='Templates')


plt.savefig("decodingtimes.pdf")

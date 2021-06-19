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
from helloworld_Encode_Decode import EncodeDecode
from helloworld_ONT_decoding  import Decode_ONT
from helloworld_MiSeq_analyses import MiSeq_Analysis

## Encode_Decode:

# Demonstration of Encoding from Ternary Values

hello_world = ["00010212","00110202","00211000","01011000","01111010","01201012",
               "02011102","02111010","02211020","10011000","10110201","10201020"]

for string in hello_world:
    print(string,''.join(EncodeDecode.encodeStr(string)))

####
#### Decoding by 2-step filter rule:
#### 1. Template architecture (length, terminal base)
#### 2. Most frequent
####
data = pd.read_csv("H01-H12/H01-H12_dataTable.csv",sep=",", header=0,
           dtype={"template_ID":str,"match":int,"template":str,"strandC":str,"strandR":str,"strandR_len":int,"strandC_len":int,"template_align":str,"strand_align":str})

data["possiblehit"] = data.apply(Decode_ONT.findpossiblehits,axis=1)

templates_to_decode = ['H01','H02','H03','H04','H05','H06','H07','H08','H09','H10','H11','H12']

for template in templates_to_decode:
    decoding(data,template)






## ONT_Decoding:
data = pd.read_csv("H01-H12_dataTable_ONT.csv",sep=",", header=0,names=["template_ID","match","template","strandC","strandR","time","read"])
data["time"] = pd.to_datetime(data["time"])

# start is defined as the time stamp of the first read

data["hr_from_start"] = (data["time"].subtract(min(data["time"])))/np.timedelta64(1, 'h')
data["min_from_start"] = (data["time"].subtract(min(data["time"])))/np.timedelta64(1, 'm')
data["sec_from_start"] = (data["time"].subtract(min(data["time"])))/np.timedelta64(1, 's')
data.fillna('',inplace=True)

data["possiblehit"] = data.apply(findpossiblehits,axis=1)

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

H01_decoding = decodingprobs(H01,hrs[1:],iters)
H02_decoding = decodingprobs(H02,hrs[1:],iters)
H03_decoding = decodingprobs(H03,hrs[1:],iters)
H04_decoding = decodingprobs(H04,hrs[1:],iters)
H05_decoding = decodingprobs(H05,hrs[1:],iters)
H06_decoding = decodingprobs(H06,hrs[1:],iters)
H07_decoding = decodingprobs(H07,hrs[1:],iters)
H08_decoding = decodingprobs(H08,hrs[1:],iters)
H09_decoding = decodingprobs(H09,hrs[1:],iters)
H10_decoding = decodingprobs(H10,hrs[1:],iters)
H11_decoding = decodingprobs(H11,hrs[1:],iters)
H12_decoding = decodingprobs(H12,hrs[1:],iters)


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





## MiSeq Analysis:
data = pd.read_csv("H01-H12_dataTable.csv",sep=",", header=0,
           dtype={"template_ID":str,"match":int,"template":str,"strandC":str,"strandR":str,"strandR_len":int,"strandC_len":int,"template_align":str,"strand_align":str})

#extract the perfect strands for each template

H01_perfect = data.loc[(data["template_ID"]=='H01') & (data["match"]==1)]
H02_perfect = data.loc[(data["template_ID"]=='H02') & (data["match"]==1)]
H03_perfect = data.loc[(data["template_ID"]=='H03') & (data["match"]==1)]
H04_perfect = data.loc[(data["template_ID"]=='H04') & (data["match"]==1)]
H05_perfect = data.loc[(data["template_ID"]=='H05') & (data["match"]==1)]
H06_perfect = data.loc[(data["template_ID"]=='H06') & (data["match"]==1)]
H07_perfect = data.loc[(data["template_ID"]=='H07') & (data["match"]==1)]
H08_perfect = data.loc[(data["template_ID"]=='H08') & (data["match"]==1)]
H09_perfect = data.loc[(data["template_ID"]=='H09') & (data["match"]==1)]
H10_perfect = data.loc[(data["template_ID"]=='H10') & (data["match"]==1)]
H11_perfect = data.loc[(data["template_ID"]=='H11') & (data["match"]==1)]
H12_perfect = data.loc[(data["template_ID"]=='H12') & (data["match"]==1)]


# run length encode each perfect strand individually

H01_perfectruns = rle(H01_perfect)
H02_perfectruns = rle(H02_perfect)
H03_perfectruns = rle(H03_perfect)
H04_perfectruns = rle(H04_perfect)
H05_perfectruns = rle(H05_perfect)
H06_perfectruns = rle(H06_perfect)
H07_perfectruns = rle(H07_perfect)
H08_perfectruns = rle(H08_perfect)
H09_perfectruns = rle(H09_perfect)
H10_perfectruns = rle(H10_perfect)
H11_perfectruns = rle(H11_perfect)
H12_perfectruns = rle(H12_perfect)


#plot the extension length for each nucleotide, for all perfect strands of H01-H12

plt.rcParams['figure.figsize'] = 4,6
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['MyriadPro']
sns.set(font_scale=3)

sns.set_style("white")
strands = {"H01":"AGACTGTG","H02":"AGTAGCTG", "H03":"AGCGTCTC", "H04":"ACTACTCT", "H05":"ACGTAGTC", "H06":"ACAGTCGC",
        "H07":"ATCGTAGC","H08":"ATACGACT","H09":"ATGTAGCT","H10":"TCTACTCT","H11":"TCGTCAGT","H12":"TCAGTCAG"}

basecolor = {"A":"#d7191c","C":"#2c7bb6","G":"#fdae61","T":"#abd9e9"}


sorder = ['H01','H02','H03','H04','H05','H06','H07','H08','H09','H10','H11','H12']

fig, axes = plt.subplots(nrows=2, ncols=6, sharex=False, sharey=True, figsize=(20,10))
axes_list = axes.flatten()

for (s,a) in zip(sorder,axes_list):

    #create the color palette based on the bases
    num = 0
    my_pal = {}
    colors = []
    currStrand = strands[s]
    for char in currStrand:
        my_pal[str(num)] = basecolor[char]
        colors.append(basecolor[char])
        num+=1

    currData = eval(s+"_perfectruns")

    if(len(currData.columns)>8):
        currData.drop(currData.columns[[-1,]], axis=1, inplace=True)

    g=sns.lvplot(data=currData,ax=a,palette=my_pal)
    g.set(xticklabels=list(currStrand))
    g.set(title=s)

    [t.set_color(i) for (i,t) in zip(colors,a.xaxis.get_ticklabels())]

    a.tick_params(which='both',bottom='off',left='off',right='off',top='off')
    a.spines['left'].set_visible(False)
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)

fig.text(0.5, -0.01, 'Strands', ha='center')
fig.text(-0.01, 0.5, 'Extension length (bases)', va='center', rotation='vertical')
fig.suptitle('Extension lengths for perfect strands of H01-H12', fontsize=24, y=1)

plt.tight_layout()
plt.savefig("helloworld_perfect_extensionLengths.pdf")


# trim for the sequences that are 9 long

trimmed_data = {'trimmedlen':[],'trimmedcompressed':[]}
data_noNA = data.dropna()
data_noNA.reset_index(inplace=True)

for index, row in data_noNA.iterrows():
    (curr_rawlen,curr_compressedlen) = snipTrailing(row)
    trimmed_data['trimmedlen'].append(curr_rawlen)
    trimmed_data['trimmedcompressed'].append(curr_compressedlen)

trimmed_data_df = pd.DataFrame.from_dict(trimmed_data)
data_withTrim = pd.concat([data_noNA,trimmed_data_df],axis=1)


# extract only the perfect strands and tabulate extension lengths for each transition

data_withTrim_perfects = data_withTrim.loc[data_withTrim["match"]==1]
trans = countTransitions(data_withTrim_perfects)

df = pd.DataFrame.from_dict(trans, orient='index')
df_melt=pd.melt(df.transpose())
df_melt_dropna = df_melt.dropna()
df_melt_dropna.reset_index(inplace=True)
df_melt_dropna.drop(columns=['index'],inplace=True)


# Sweet fig
plt.rcParams['figure.figsize'] = 7.5,5
sns.set(font_scale=1.5)
sns.set_style("white")

g = sns.lvplot(x="variable", y="value", data=df_melt_dropna, color='#D8D8D8')
g.set(xlabel='transitions (previous>next)')
g.set(ylabel='Extension length (bases)')
g.set(title='Extension lengths for each transition in H01-H12')
sns.despine()

plt.savefig('perfectTransitions.pdf')


# visualize the increase in raw lengths as a function of template length (or cycles)

data_withTrim.drop(columns=['index'],inplace=True)
data_withTrim_compressed8 = data_withTrim.loc[(data_withTrim["trimmedcompressed"]<=8)&(data_withTrim["trimmedcompressed"]>0)]

sns.set(font_scale=1)
sns.set_style("white")

fig, axes = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize=(5,5))
g=sns.lvplot(x="trimmedcompressed",y="trimmedlen",data=data_withTrim_compressed8,color='#D8D8D8',ax=axes)
g.set(xlabel='Compressed length (bases)')
g.set(ylabel='Raw length (bases)')
g.set(title='Length of all synthesized strands over cycles ')

plt.savefig('Length_raw_to_cycles.pdf')


# raw lengths of perfect strands vs all synthesized strands

plt.rcParams['figure.figsize'] = 6,3
sns.set(font_scale=1)
sns.set_style("white")

rawlengths = data_withTrim["trimmedlen"].astype("float")
matches = data_withTrim.loc[data_withTrim["match"]==1, "trimmedlen"].astype("float")

g=sns.kdeplot(matches,color="black",bw=0.75)
g=sns.kdeplot(rawlengths,color='gray',bw=0.75,shade=True)
g.set(xlim=(0, 60))
g.set(xlabel='Raw length (bases)')
g.set(ylabel='Density')
g.set(title='Raw lengths for all (gray) and perfect (black) strands for H01-H12')
g.legend_.remove()

plt.tight_layout()

plt.savefig("helloworld_bulklengths.pdf")


#plot the trimmed length distributions

sns.set(font_scale=1)
sns.set_style("white")

sorder = ['H01','H02','H03','H04','H05','H06','H07','H08','H09','H10','H11','H12']

fig, axes = plt.subplots(nrows=4, ncols=3, sharex=False, sharey=True, figsize=(10,8))
axes_list = axes.flatten()

for (s,a) in zip(sorder,axes_list):
    currSample = data_withTrim.loc[data_withTrim["template_ID"]==s]
    rawlengths = currSample["trimmedlen"].astype("float")
    matches = currSample.loc[currSample["match"]==1, "trimmedlen"].astype("float")

    g=sns.kdeplot(rawlengths,shade=True,bw=0.75,color="gray",ax=a)
    g=sns.kdeplot(matches, bw=0.75,color="black",linestyle='dashed',ax=a)
    g.set(title=s)
    g.legend_.remove()
    g.set(xlim=(0, 75),ylim=(0, 0.1))
    a.tick_params(which='both',bottom='off',left='off',right='off',top='off')
    a.spines['top'].set_visible(False)
    a.spines['right'].set_visible(False)

    fig.text(0.5, -0.01, 'Raw length (bases)', ha='center')
    fig.text(-0.01, 0.5, 'Density', va='center', rotation='vertical')
    fig.suptitle('Raw lengths for all (gray) and perfect (dotted) strands for each of H01-H12', fontsize=16, y=1.01)

    plt.tight_layout()
plt.savefig("helloworld_lengthdist.pdf")



# tabulate errors for all synthesized strands

numBases = list(range(0,9))
all_errors = bulkErrors_count(data)

# plot errors for all synthesized strands in bulk

plt.rcParams['figure.figsize'] = 8,8
sns.set(font_scale=1.3)
sns.set_style("white")

numBases = list(range(0,9))

mismatches = all_errors['mismatches']
mismatches_normalized = [entry/all_errors['total'] for entry in mismatches]

insertions = all_errors['insertions']
insertions_normalized = [entry/all_errors['total'] for entry in insertions]

deletions = all_errors['deletions']
deletions_normalized = [entry/all_errors['total'] for entry in deletions]

fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=False, figsize=(4,8))
axes_list = axes.flatten()

axes_list[0].hist(numBases,9, weights = mismatches_normalized, edgecolor='black', linewidth=1,color='#D8D8D8')
axes_list[0].set_ylim(0,1)
axes_list[0].set_xlim(0,9)
sns.despine(ax=axes_list[0])
axes_list[0].set_title('Mismatch')

axes_list[1].hist(numBases,9, weights = insertions_normalized,edgecolor='black', linewidth=1,color='#D8D8D8')
axes_list[1].set_ylim(0,1)
axes_list[1].set_xlim(0,9)
sns.despine(ax=axes_list[1])
axes_list[1].set_title('Insertion')

axes_list[2].hist(numBases,9, weights = deletions_normalized,edgecolor='black', linewidth=1,color='#D8D8D8')
axes_list[2].set_ylim(0,1)
axes_list[2].set_xlim(0,9)
sns.despine(ax=axes_list[2])
axes_list[2].set_title('Missing')
axes_list[2].set_xlabel('Bases in strandC')


fig.text(-0.01, 0.5, 'Fraction of all strands', va='center', rotation='vertical')
fig.suptitle('Errors for all H01-H12 strands', fontsize=14, y=1)

plt.tight_layout()
plt.savefig("helloworld_bulkerrors.pdf")


# select all strands for each template sequence individually

H01 = data.loc[data["template_ID"]=='H01']
H02 = data.loc[data["template_ID"]=='H02']
H03 = data.loc[data["template_ID"]=='H03']
H04 = data.loc[data["template_ID"]=='H04']
H05 = data.loc[data["template_ID"]=='H05']
H06 = data.loc[data["template_ID"]=='H06']
H07 = data.loc[data["template_ID"]=='H07']
H08 = data.loc[data["template_ID"]=='H08']
H09 = data.loc[data["template_ID"]=='H09']
H10 = data.loc[data["template_ID"]=='H10']
H11 = data.loc[data["template_ID"]=='H11']
H12 = data.loc[data["template_ID"]=='H12']

# tabulate errors for strands synthesized per template

H01_errors = bulkErrors_count(H01)
H02_errors = bulkErrors_count(H02)
H03_errors = bulkErrors_count(H03)
H04_errors = bulkErrors_count(H04)
H05_errors = bulkErrors_count(H05)
H06_errors = bulkErrors_count(H06)
H07_errors = bulkErrors_count(H07)
H08_errors = bulkErrors_count(H08)
H09_errors = bulkErrors_count(H09)
H10_errors = bulkErrors_count(H10)
H11_errors = bulkErrors_count(H11)
H12_errors = bulkErrors_count(H12)


# plot errors for strands synthesized per template sequence

import matplotlib.ticker as plticker

plt.rcParams['figure.figsize'] = 8,8
sns.set(font_scale=1)
sns.set_style("white")

numBases = list(range(0,9))

errors_to_plot = [H01_errors, H02_errors, H03_errors, H04_errors, H05_errors, H06_errors,
                 H07_errors, H08_errors, H09_errors, H10_errors, H11_errors, H12_errors]

fig, axes = plt.subplots(nrows=3, ncols=len(errors_to_plot), sharex=True, sharey=True, figsize=(15,7))
axes_list = axes.flatten()
i=0

loc = plticker.MultipleLocator(base=2.0)

for strands in errors_to_plot:

    curr_mismatches = strands['mismatches']
    curr_mismatches_normalized = [entry/strands['total'] for entry in curr_mismatches]

    curr_insertions = strands['insertions']
    curr_insertions_normalized = [entry/strands['total'] for entry in curr_insertions]

    curr_deletions = strands['deletions']
    curr_deletions_normalized = [entry/strands['total'] for entry in curr_deletions]

    axes_list[i].hist(numBases,9, weights = curr_mismatches_normalized, edgecolor='black', linewidth=1,color='#D8D8D8')
    axes_list[i].set_ylim(0,1)
    axes_list[i].set_xlim(0,9)
    axes_list[i].xaxis.set_major_locator(loc)
    sns.despine(ax=axes_list[i])

    axes_list[i+1*len(errors_to_plot)].hist(numBases,9, weights = curr_insertions_normalized,edgecolor='black', linewidth=1,color='#D8D8D8')
    axes_list[i+1*len(errors_to_plot)].set_ylim(0,1)
    axes_list[i+1*len(errors_to_plot)].set_xlim(0,9)
    axes_list[i+1*len(errors_to_plot)].xaxis.set_major_locator(loc)
    sns.despine(ax=axes_list[i+len(errors_to_plot)])

    axes_list[i+2*len(errors_to_plot)].hist(numBases,9, weights = curr_deletions_normalized,edgecolor='black', linewidth=1,color='#D8D8D8')
    axes_list[i+2*len(errors_to_plot)].set_ylim(0,1)
    axes_list[i+2*len(errors_to_plot)].set_xlim(0,9)
    axes_list[i+2*len(errors_to_plot)].xaxis.set_major_locator(loc)
    sns.despine(ax=axes_list[i+2*len(errors_to_plot)])

    i+=1

fig.text(0.5, -0.01, 'Bases of strandC', ha='center')
fig.text(-0.01, 0.5, 'Fraction of all strands', va='center', rotation='vertical')
fig.suptitle('Mismatches, Insertions, Missing (rows) for all H01-H12 (cols)', fontsize=16, y=1)

plt.tight_layout()
plt.savefig("helloworld_individualerrors.pdf")

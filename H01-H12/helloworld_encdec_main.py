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
from Classes.helloworld_Encode_Decode import EncodeDecode

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
data = pd.read_csv("H01-H12/Data/H01-H12_dataTable.csv",sep=",", header=0,
           dtype={"template_ID":str,"match":int,"template":str,"strandC":str,"strandR":str,"strandR_len":int,"strandC_len":int,"template_align":str,"strand_align":str})

data["possiblehit"] = data.apply(EncodeDecode.findpossiblehits,axis=1)

templates_to_decode = ['H01','H02','H03','H04','H05','H06','H07','H08','H09','H10','H11','H12']

for template in templates_to_decode:
    EncodeDecode.decoding(data,template)

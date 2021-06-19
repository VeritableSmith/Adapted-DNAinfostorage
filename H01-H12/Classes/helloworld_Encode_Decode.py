import glob, os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns

class EncodeDecode:

    #function input previous base, desired trit
    @classmethod
    def nextBase(cls, prevBase, wantedTrit):
        transitionTable = {'A0':'G', 'A1':'C', 'A2':'T','C0':'T', 'C1':'G', 'C2':'A','G0':'A', 'G1':'T', 'G2':'C','T0':'C', 'T1':'A', 'T2':'G'}

        return transitionTable[prevBase+str(wantedTrit)]

    @classmethod
    def encodeStr(cls, inputInfo_string):
        inputInfo = [int(s) for s in list(inputInfo_string)]
        outStr = []
        for i in range(len(inputInfo)):
            if(i==0):
                prevBase = 'G'
            else:
                prevBase = outStr[i-1]
            currtrit = inputInfo[i]

            outStr.append(cls.nextBase(prevBase,currtrit))
        return outStr

    @classmethod
    def findpossiblehits(cls, row):
        curr_seq = row['strandC']
        desired_seq = row['template']

        # filter rule - look for strands of a set length and the desired terminal 'C'
        if((len(curr_seq) == len(desired_seq)) and len(curr_seq)>0 and (curr_seq[-1] is desired_seq[-1])):
            return 1
        else:
            return 0

    @classmethod
    def decoding(cls, df, templateID):
        currdf = df[df["template_ID"]==templateID]
        currdf_possiblehits  = currdf[currdf["possiblehit"]==1]
        print(templateID,currdf.iloc[0]["template"])
        print(currdf_possiblehits['strandC'].value_counts(ascending=False).head(5).to_frame(),"\n\n")

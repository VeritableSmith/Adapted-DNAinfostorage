import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set_color_codes("bright")
import glob
from pylab import rcParams

class Decode_ONT:
    def findpossiblehits(row):
        curr_seq = row['strandC']
        desired_seq = row['template']

        # filter rule - look for strands of a set length and the desired terminal 'C'
        if((len(curr_seq) == len(desired_seq)) and len(curr_seq)>0 and (curr_seq[-1] is desired_seq[-1])):
            return 1
        else:
            return 0

    def decodingprobs(df,hrs,num_iterations):
        df_prob_by_hrs = []
        for timepts in hrs:
            num_matches = []
            for i in range(num_iterations):
                mini = df.sample(n=df.loc[df["hr_from_start"]<=timepts].shape[0])
                mini_onlypossibles = mini
                mini_summary = mini_onlypossibles[["strandC","match","possiblehit"]].groupby(['strandC'], as_index=False).agg(['sum']).reset_index()

                newcols = [''.join(t) for t in mini_summary]
                mini_summary.columns = newcols
                mini_summary.reset_index(level=0, inplace=True)

                all_max = mini_summary["possiblehitsum"].max()
                onlymax = mini_summary[mini_summary["possiblehitsum"]==all_max]

                if((onlymax.shape[0] == 1) and (onlymax.iloc[0]["matchsum"])>0):
                    num_matches.append(1)
                else:
                    num_matches.append(0)
            df_prob_by_hrs.append(sum(num_matches)/num_iterations)
        return df_prob_by_hrs

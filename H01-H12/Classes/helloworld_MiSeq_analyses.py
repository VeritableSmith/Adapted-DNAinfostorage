import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pylab import rcParams
import numpy as np
import glob, os
from re import sub

class MiSeq_Analysis:
    @classmethod
    def rle(cls,df):
        run_counts = {'0':[],'1':[],'2':[],'3':[],'4':[],'5':[],'6':[],'7':[]}
        for index, row in df.iterrows():
            curr_rawseq = row['strandR']
            curr_rle = sub(r'(.)\1*', lambda m: ","+str(len(m.group(0)))+"," + m.group(1),curr_rawseq)
            curr_rle_tolist = curr_rle[1:].split(',')
            curr_onlyNumbersString = curr_rle_tolist[::2]
            curr_onlyNumbers = [int(x) for x in curr_onlyNumbersString]
            for i in range(8):
                run_counts[str(i)].append(curr_onlyNumbers[i])
        run_counts_df = pd.DataFrame.from_dict(run_counts)
        return run_counts_df

    # tabulate the extension lengths for each of the 12 types of transitions
    @classmethod
    def countTransitions(cls,df):
        transitions = {'AC':[],'AG':[],'AT':[],'CA':[],'CG':[],'CT':[],
                    'GA':[],'GC':[],'GT':[],'TA':[],'TC':[],'TG':[]}

        for index, row in df.iterrows():
            text = row.strandR
            currstr = sub(r'(.)\1*', lambda m: ","+str(len(m.group(0)))+"," + m.group(1),text)
            currstr_tolist = currstr[1:].split(',')
            onlyNumbersString = currstr_tolist[::2]
            onlyChars = currstr_tolist[1::2]
            onlyNumbers = [int(x) for x in onlyNumbersString]

            prev = 'G'
            for (char,num) in zip(onlyChars,onlyNumbers):
                if(prev != char): #there is a GG case where the initiator ends in G and first base is G
                    currTrans = prev+char
                    transitions[currTrans].append(num)
                prev=char
        return transitions

    #for those sequences that are 9 long, we need to remove the last extension from each sequence
    #so that our length data is consistent
    @classmethod
    def snipTrailing(cls,row):
        template = row["template_ID"]
        if(len(template)>8):
            curr_seqraw = row["strandR"]
            currstr = sub(r'(.)\1*', lambda m: ","+str(len(m.group(0)))+"," + m.group(1),curr_seqraw)
            currstr_tolist = currstr[1:].split(',')
            onlyNumbersString = currstr_tolist[::2]
            onlyNumbers = [int(x) for x in onlyNumbersString]
            del onlyNumbers[-1]
            trimmedlenraw = sum(onlyNumbers)
            trimmedlencompressed = len(onlyNumbers)
        else:
            trimmedlenraw = row["strandR_len"]
            trimmedlencompressed = len(row.strandC)
        return (trimmedlenraw,trimmedlencompressed)

    # tabulates mismatches, insertions, and missing from pre-computed alignments
    @classmethod
    def bulkErrors_count(cls, df):
        #to account for 0s
        matches = [0]*(9)
        mismatches = [0]*(9)
        deletions = [0]*(9)
        insertions = [0]*(9)

        num_analyzed = 0
        for index, row in df.iterrows():
            currMatches, currMismatches, currDeletions, currInsertions = 0,0,0,0
            template = row['template']
            template_aligned = list(row['template_align'])
            query_aligned = list(row['strand_align'])

            if(len(template)>8):
                #find terminal C index
                for i in range(len(template)):
                    if(template_aligned[-i]=='C'):
                        break
                template_aligned = template_aligned[:-i]
                query_aligned = query_aligned[:-i]
            num_analyzed+=1
            for char_template, char_query in zip(template_aligned, query_aligned):
                action = ''
                if(char_template == char_query): #match
                    currMatches+=1
                    action = 'match'
                elif(char_template == '-'):   #insertion
                    currInsertions+=1
                    action = 'insertion'
                elif((char_template != char_query) and char_query != '-'): #mismatch if not deletion
                    currMismatches+=1
                    action = 'mismatch'
                elif(char_query == '-'): #deletion in query
                    currDeletions+=1
                    action = 'deletion'
            matches[currMatches]+=1
            mismatches[currMismatches]+=1
            deletions[currDeletions]+=1
            insertions[currInsertions]+=1

        bulkData = {'matches':matches,'mismatches':mismatches,'deletions':deletions,'insertions':insertions, 'total':num_analyzed}
        return bulkData

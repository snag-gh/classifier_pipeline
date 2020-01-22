#! /usr/bin/env python3

from sys import argv
import pandas as pd
from pprint import pprint as pp
import os

df = pd.read_csv(argv[1], sep = '\t')
df['seg.length'] = df['loc.end'] - df['loc.start'] + 1
threshold_1p = -0.2
threshold_19q = -0.273
chr1_breakpoint = 120425000  #lowest observed end point for chr1 p arm 
#chr19_breakpoint = 29050000  #highest observed starting point for chr19 q arm
chr19_breakpoint = 28400000  #highest observed starting point for chr19 q arm
status = 'not codeleted'

def weighted_mean(df):
    arm_len = df[['seg.length']].sum()[0]
    #pp(arm_len)
    df['frac'] = (df['seg.length'] / arm_len) * df['seg.mean']
    #pp(df)
    wtmean = df[['frac']].sum()[0]
    #pp(wtmean)
    return wtmean


def find_mean(df, chrm, breakpoint):
    #pp(df)
    if chrm == 'chr1':
        seg_gt_brk = df[(df['chrom'] == chrm) & (df['loc.end'] >= breakpoint)]
        #pp(seg_gt_brk)
        partial = 0                
        if not seg_gt_brk.empty:
            pend = seg_gt_brk.loc[seg_gt_brk['loc.end'].idxmin()]
            #pp(pend)
            pend_loc = df.index[(df['loc.start'] == pend['loc.start']) & (df['loc.end'] == pend['loc.end'])][0]   #line number of chr p arm end in the .seg file
            #pp(pend_loc)
            beg = df[df.chrom == chrm].index[0]
            #print(beg)
            p_arm = df.iloc[beg:pend_loc+1].copy()
            frag_ge1m = p_arm.loc[(p_arm['seg.length'] >= 1000000) & (p_arm['seg.mean'] > -0.2)] #Check if any fragment larger than 1Mbp has normal or high copy number
            #pp(p_arm)
            if not frag_ge1m.empty:
                partial = 1
            wtmean = weighted_mean(p_arm)
            return wtmean, partial
        else:
            return None, None
    elif chrm == 'chr19':
        seg_lt_brk = df[(df['chrom'] == chrm) & (df['loc.start'] <= breakpoint)]
        #pp(seg_lt_brk)
        partial = 0
        if not seg_lt_brk.empty:
            qstart = seg_lt_brk.loc[seg_lt_brk['loc.start'].idxmax()]
            #pp(qstart)
            qstart_loc = df.index[(df['loc.start'] == qstart['loc.start']) & (df['loc.end'] == qstart['loc.end'])][0]
            #pp(qstart_loc)
            end = df[df.chrom == chrm].index[-1]
            #print(end)
            q_arm = df.iloc[qstart_loc:end+1].copy()
            #pp(q_arm)
            frag_ge1m = q_arm.loc[(q_arm['seg.length'] >= 1000000) & (q_arm['seg.mean'] > -0.273)]
            if not frag_ge1m.empty:
                partial = 1
            wtmean = weighted_mean(q_arm)
            return wtmean, partial
        else:
            return None, None     
            

if argv[2] == "EPIC":
    chr1_meanp, partial1p = find_mean(df, 'chr1', chr1_breakpoint)
#    print(argv[1], 'chr1', chr1_mean_p, chr1_mean_q, chr1_dist)
    chr19_meanq, partial19q = find_mean(df, 'chr19', chr19_breakpoint)
    if (chr1_meanp <= threshold_1p) & (chr19_meanq <= threshold_19q) & (partial1p == 0) & (partial19q == 0):
        status = 'codeleted'
    print('Sample,codeletion status,chr1p,chr1p_partial_deletion_status,chr19q,chr19q_partial_deletion_status,threshold_1p,threshold_19q')
    print('{},{},{},{},{},{},{},{}'.format(str(os.path.basename(argv[1])), status, chr1_meanp, partial1p, chr19_meanq, partial19q, threshold_1p, threshold_19q))


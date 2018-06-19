"""

* File Name : soc_analysis_lib.py

* Purpose : Lib contains SOC data preprocessing and statistical analysis

* Creation Date : 18-07-2017

* Last Modified : Thu 20 Jul 2017 05:07:38 PM EDT

* Created By : Shijie Shu

"""

import pandas as pd
import numpy as np


def getvarxls(data, var, profid, loc):
    """ Get variable from synthesis xlssheet corresponds to proifle ID
    input: 
        data: {DataFrame} whole data sheet read into DF
        var: {string} variable of interest
        profid: {array} 1-D array of profile ID
        loc: 0 for top value, -1 for last value, or ':' for all values
    output: 
        1-D array (same length of profid) of var
    """
    if isinstance(profid, int):
        l = data.loc[profid][var]
    else:
        l = [data.loc[i:i][var] for i in profid if ~np.isnan(i)]
    mat = []
    for i in l:
        try: # if has multi-layers
            if isinstance(loc, str):
                dum = eval('i.iloc[' + loc + '].values')                
                try: 
                    mat.append(list(dum.astype(float)))
                except ValueError: # is unicode
                    mat.append(list(dum))
            else:
                mat.append(i.iloc[loc])
        except AttributeError:
            if isinstance(loc, str):
                try:
                    mat.append([float(i)])
                except ValueError:
                    mat.append([i])
            else:
                mat.append(i)       
    if isinstance(loc, str):
        return np.array(sum(mat, []))
    else:        
        return np.array(mat)    
 
def subset_by_cols(data, colname, val, tof):
    """ Extract a subset of a pandas data frame if value of the sepcified column matches
    given value
    input: 
        data: {DataFrame} From where extract the subset
        colname: The column name being specified
        val: The given value being used as criteria
        tof: A logical flag to identify weather the given value need to be equaled or excluded
    output: 
        subset: the extracted subset of the data
        sub_profid: the profile id of the subset
    """

    col_idx = [i for i, x in enumerate(data.columns == colname) if x][0]
    if(tof):
        subset = data[data.iloc[:,col_idx] == val]
    else:
        subset = data[data.iloc[:,col_idx] != val]
    # profid = data.index.unique()
    sub_profid = subset.index.unique()

    return subset, sub_profid

def subset_by_ids(data, sel_profid):
    """ Extract a subset of a pandas data frame by profid value 
    input: 
        data: {DataFrame} From where extract the subset, must contain an index column
        sel_profid: The profile ID list being selected
    output: 
        subset: the extracted subset of the data
    """

    lgc = np.ones((len(data.index)), dtype=bool)
    for i in range(0, len(data.index)):
        lgc[i] = data.index[i] in sel_profid
    subset = data[lgc]

    return subset

def choose_best_prof(data, prof_std):
    """ Select the profile matches one standard profile the best
    input: 
        data: {pd.DataFrame} The profile collection, must contain an index column
        prof_std: The standard profile, must have the same length as prof
    output: 
        score: the score describing the "matchness"
    """

    score_sel = 0.
    id_sel = 0
    for i in range(0,len(data.index)):
        score = 0.
        for j in range(0, len(prof_std)):
            score = score + np.absolute(data.iloc[i,j] - prof_std[j])
        if(score > score_sel):
            score_sel = score
            id_sel = data.index[i]
     
    return score_sel, id_sel

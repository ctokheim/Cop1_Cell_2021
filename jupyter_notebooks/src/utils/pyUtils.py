#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : MIT
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 27/01/2020
# Last Modified Date: 10/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>


import os
try:
    import cPickle as pickle
except:
    import pickle
import pandas as pd
import numpy as np
from Biopyutils import Comm
import logging
import statsmodels
import statsmodels.stats.multitest as multi
import statsmodels.api as sm
import sys
sys.path.append("..")
from functools import reduce
from scipy import stats, linalg

def addSymbolOnProtein(df,protien_col):
    """addSymbolOnProtein
    Annotate Protein with their corresponding gene name.
    The <protein_col> will be modified as <protein ID>-<Hugo Symbol>

    Parameters
    ----------
    df : pd.DataFrame
        df is a dataframe with protein ID need to map
    protien_col : str
        protien_col is the column name of protien ID in the df

    Returns:
    ----------
    pd.DataFrame
       The <protein_col> will be modified as <protein ID>-<Hugo Symbol>
    """
    from Tigger.configs.data_configs  import hg38_ANNO_PROTEIN
    ref = pd.read_csv(hg38_ANNO_PROTEIN,sep='\t')[['Entry','EntryName']]
    ref.EntryName = ref.EntryName.map(lambda x: x.split('_')[0])
    ref.rename(columns={'Entry':protien_col,'EntryName':'Hugo Symbol'},inplace=True)
    result = df.merge(ref,on=protien_col)
    result[protien_col] = result[protien_col] +'-'+ result['Hugo Symbol']
    result.drop('Hugo Symbol',axis=1)
    return result


def addGeneName(df):
    """ Add corresponding gene name according to the Symbol name on the input table

    Parameters
    ----------
    df : pd.DataFrame
        Must contains a column named 'Symbol

    Returns
    -------
    pd.DataFrame
        With additional column with exact gene name be added in.
    """
    from Tigger.configs.data_configs import hg38_ANNO_GENE

    ref = pd.read_csv(hg38_ANNO_GENE, sep='\t')[['Symbol', 'name']]
    ref.rename(columns={'name':'Gene Name'},inplace=True)
    result = df.merge(ref,on='Symbol',how='left')
    return result

def quantileNorm(df):
    """ Quantile normalization across samples

    Parameters
    ----------
    df : pd.DataFrame
        indexed by features and samples are in the column

    Returns
    -------
    pd.DataFrame
        A dataframe that has been quantile normalized
    """
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    result = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return result

def getExprsnSig(df, exprsn):
    """ Calculation signature from gene expression matrix.
    (average expression if we don't have weight, correlation if we have weights)
    Parameters
    ----------
    df : list / pd.Series
        Weighted or un-weighted gene list
    exprsn: pd.DataFrame
        A gene expression profile, indexded by Hugo Symbol and columned by samples
    Returns
    -------
    pd.Series
        Signature for samples

    Raises
    ------
    ValueError
        input df is not valid
    """
    exprsn = Comm.idConvert(df=exprsn,species='hg',map_id='Entrez')
    if isinstance(df, list):
        ol_gene = exprsn.index.intersection(df)
        miss_gene = float(len(df) - len(ol_gene))
        if miss_gene/len(df) > .5:
            raise ValueError('%.1f %% input genes are missing.' %
                             ((miss_gene/len(df)) * 100))

        sig = exprsn.loc[ol_gene, :].mean(axis=0)
    elif isinstance(df, pd.Series):
        ol_gene = exprsn.index.intersection(df.index)
        miss_gene = float(df.size - len(ol_gene))
        if miss_gene/df.size > .5:
            raise ValueError('%.1f %% input genes are missing.' %
                             ((miss_gene/df.size) * 100))
        sig = exprsn.loc[ol_gene, :].corrwith(df[ol_gene].T)
    else:
        raise ValueError(
            '%s data type is not available, only accept str or list data type' % type(df))

    return sig

def normTumorNormal(df, normal_samples, tumor_samples, min_normal_require=3):
    """ Normalize input expression profile
    1. Quantile normalization
    2. Subtract mean across normal samples (>min_normal_required) or mean across tumor sample
    3. Return normalize expression profile for tumor sample only
    Parameters
    ----------
    df : pd.DataFrame
        indexed by features and samples are in the column
    normal_samples : list
        A name list of normal samples
    tumor_samples : list
        A name list of tumor samples
    min_normal_require : int, optional
        minimum number of normal samples required to subtract mean of normal instead of mean of tumor samples,
        by default 3
    """
    select_samples = normal_samples + tumor_samples
    norm_df = quantileNorm(df.loc[:, select_samples])
    if len(normal_samples) >= min_normal_require:
        norm_factor = norm_df[normal_samples].mean(axis=1)
    else:
        norm_factor = norm_df[tumor_samples].mean(axis=1)
    # subtract mean across samples
    norm_df = norm_df[tumor_samples].sub(norm_factor, axis='index')
    # divide std across samples
    #norm_df = norm_df.div(norm_df.std(axis=1), axis='index')

    return norm_df

def fdrAdjust(p_value, method='fdr_bh'):
    """ Adjust p value

    Parameters
    ----------
    p_value : pd.Series
        Original p value
    method : str, optional
        Method used for testing and adjustment of pvalues. Default: "fdr_bh"
        Can be either the full name or initial letters. Available methods are:
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)

    Returns
    -------
    pd.Series
        p-values corrected for multiple tests
    """
    adjust_p = multi.multipletests(p_value, method=method)[1]
    return adjust_p

def ppcor(C, scale=True):
    """
    Returns the sample linear partial correlation coefficients between pairs of variables in C, controlling
    for the remaining variables in C.
    Parameters
    ----------
    C : array-like, shape (n, p)
        Array with the different variables. Each column of C is taken as a variable
    scale: Bolean, optional
        When True, standardize data before calculating correlation (which will returen identical results of R ppcor).
        Default: True
    Returns
    -------
    P_corr : array-like, shape (p, p)
        P[i, j] contains the partial correlation of C[:, i] and C[:, j] controlling
        for the remaining variables in C.
    P_p : array-like, shape (p, p)
        P[i, j] contains the p value of partial correlation of C[:, i] and C[:, j] controlling
        for the remaining variables in C.
    """

    C = np.asarray(C)
    if scale is True:
        std = C.std(axis=0)
        mean = C.mean(axis=0)
        C = (C-mean)/std
    p = C.shape[1]
    P_corr = np.zeros((p, p), dtype=np.float)
    P_p = np.zeros((p, p), dtype=np.float)
    for i in range(p):
        P_corr[i, i] = 1
        for j in range(i+1, p):
            idx = np.ones(p, dtype=np.bool)
            idx[i] = False
            idx[j] = False
            beta_i = linalg.lstsq(C[:, idx], C[:, j])[0]
            beta_j = linalg.lstsq(C[:, idx], C[:, i])[0]

            res_j = C[:, j] - C[:, idx].dot(beta_i)
            res_i = C[:, i] - C[:, idx].dot(beta_j)

            corr, p_value = stats.pearsonr(res_i, res_j)
            P_corr[i, j] = corr
            P_corr[j, i] = corr

            P_p[i, j] = p_value
            P_p[j, i] = p_value

    return P_corr, P_p

def ctrlCorNoPurity( x, y, scale=True):
    ''' Calculate pearson correlation between x and y while controling ctrl.

    Parameters
    ----------
    ctrl : pandas.DataFrame
        control variable
    x : pandas.Series
        independent variable
    y : pandas.Series
        independent variable
    scale : Bolean, optional
        When True, standardize data before calculating correlation (which will returen identical results of R ppcor).
        Default: True

    Returns
    -------
    tuple
        (partial correlation, p value of partial correlation)
    '''

    ol_samples = reduce(
        lambda x, y: x.intersection(y),
        [
            x.replace([np.inf, -np.inf], np.nan).dropna().index,
            y.replace([np.inf, -np.inf], np.nan).dropna().index
        ]).unique()

    if len(ol_samples) > 10 and not all(y[ol_samples] == 0) and not all(x[ol_samples] == 0):
        C = pd.concat([
                       x[ol_samples], y[ol_samples]], axis=1)
        result = ppcor(C, scale=scale)
        cor = result[0][-1, -2]
        p = result[1][-1, -2]
    else:
        cor = np.nan
        p = np.nan
    return cor, p

def ctrlCor(ctrl, x, y, scale=True):
    ''' Calculate pearson correlation between x and y while controling ctrl.

    Parameters
    ----------
    ctrl : pandas.DataFrame
        control variable
    x : pandas.Series
        independent variable
    y : pandas.Series
        independent variable
    scale : Bolean, optional
        When True, standardize data before calculating correlation (which will returen identical results of R ppcor).
        Default: True

    Returns
    -------
    tuple
        (partial correlation, p value of partial correlation)
    '''

    ol_samples = reduce(
        lambda x, y: x.intersection(y),
        [
            ctrl.replace([np.inf, -np.inf],
                         np.nan).dropna(how='any', axis=0).index,
            x.replace([np.inf, -np.inf], np.nan).dropna().index,
            y.replace([np.inf, -np.inf], np.nan).dropna().index
        ]).unique()

    if len(ol_samples) > 10 and not all(y[ol_samples] == 0) and not all(x[ol_samples] == 0):
        C = pd.concat([ctrl.loc[ol_samples, :],
                       x[ol_samples], y[ol_samples]], axis=1)
        result = ppcor(C, scale=scale)
        cor = result[0][-1, -2]
        p = result[1][-1, -2]
    else:
        cor = np.nan
        p = np.nan
    return cor, p

def tcgaPick(df, source='tumor', transpose=False):
    ''' Picking samples by where it dissect from.

    Parameters
    ----------
    df : pandas.DataFrame
        A data frame with gene indentifier on rows and samples on columns.
    source : str, optional
        The region to pick (the default is 'tumor')
    transpose : bool, optional
        Whether the samples id is on the index (the default is False)

    Raises
    ------
    KeyError
        if input a invalid data type

    Returns
    -------
    pandas.DataFrame
        A data frame with gene indentifier on rows and picked samples on columns.
    '''

    source_map = {'tumor': '0[0-9]$',
                  'normal': '1[0-9]$'}

    if not source in source_map.keys():
        raise KeyError("""
        {0} is not a valid type of source, only accept following input: {1}
        """.format(source, ','.join(source_map.keys())))

    if transpose:
        return df.loc[df.index.str.contains(source_map[source]), :]
    else:
        return df.loc[:, df.columns.str.contains(source_map[source])]

def TcgaZscore(df, pair_TN=True):
    ''' Ways to normalize tumor expression data on TCGA

    Parameters
    ----------
    df : pandas.DataFrame
        A data frame with gene indentifier on rows and samples on columns.
        NOTICE: Value on the data frame ARE already been log scaled!
    pair_TN : bool, optional
        Tell whether normalize tumor expression profile by paired tumor samples (the default is True)

    Raises
    ------
    ValueError
        cannot find enough paired normal samples

    Returns
    -------
    pandas.DataFrame
        Normlized tumor gene expression profile
    '''

    tumor = tcgaPick(df, source='tumor')

    if pair_TN:
        normal = tcgaPick(df, source='normal')
        if normal.shape[1] > 10:
            norm_factor = normal.mean(axis=1)
        else:
            raise ValueError('Cannot find enough paired normal samples (<10)')
    else:
         norm_factor = tumor.mean(axis=1)

    result = tumor.sub(norm_factor, axis='index')
    result = result.div(result.std(axis=1), axis='index')

    return result

#################### Regression ####################

def  linearReg (df,x,y,**kwargs):
    """ Using logistic regresion to exclude the effect from covariable
    Parameters
    ----------
    df : pd.DataFrame
        A Data frame has independent variable, covariable, and dependent variable on columns
    x : str
        Column name of independent variable
    y : str
        Column name of dependent variable

    Returns
    -------
    tuple
    (z-scores,p-values)
        z score and p value of independent variable
    """
    X = df.drop(y,axis=1)
    for key, value in kwargs.items():  # include intercept into linear regression
        if key == 'intcp' and value == 1:
            X = sm.add_constant(X)
    try:
        result = sm.OLS(df[y],X).fit(disp=False)

    except np.linalg.LinAlgError:
        return (np.nan,np.nan)

    except statsmodels.tools.sm_exceptions.PerfectSeparationError:
        return (np.nan,np.nan)

    return (result.tvalues[x], result.pvalues[x])

def logitReg(df,x,y,**kwargs):
    """ Using logistic regresion to exclude the effect from covariable

    Parameters
    ----------
    df : pd.DataFrame
        A Data frame has independent variable, covariable, and dependent variable on columns
    x : str
        Column name of independent variable
    y : str
        Column name of dependent variable

    Returns
    -------
    tuple
    (z-scores,p-values)
        z score and p value of independent variable
    """
    X = df.drop(y,axis=1)
    #print(X)
    try:
        result = sm.Logit(df[y],X).fit(disp=False)

    except np.linalg.LinAlgError:
        return (np.nan,np.nan)

    except statsmodels.tools.sm_exceptions.PerfectSeparationError:
        return (np.nan,np.nan)

    return (result.tvalues[x], result.pvalues[x])

def olsReg(df, x, y, **kwargs):
    """ Using linear regresion to exclude the effect from covariable

    Parameters
    ----------
    df : pd.DataFrame
        A Data frame has independent variable, covariable, and dependent variable on columns
    x : str
        Column name of independent variable
    y : str
        Column name of dependent variable

    Returns
    -------
    tuple
    (z-scores,p-values)
        z score and p value of independent variable
    """
    X = df.drop(y, axis=1)
    try:
        result = sm.OLS(df[y], X).fit(disp=False)

    except np.linalg.LinAlgError:
        return (np.nan, np.nan)

    except statsmodels.tools.sm_exceptions.PerfectSeparationError:
        return (np.nan, np.nan)

    return (result.tvalues[x], result.pvalues[x])


def reuse_check(fileObj):
    """ Load pickle obj if file exists
    Parameters
    ----------
    fileObj : str
        fileObj is the path for requied pickle obj

    Returns:
    ----------
    None/unserialized_data
    """
    data = None
    if os.path.isfile(fileObj):
        with open(fileObj, 'rb') as handle:
            data = pickle.load(handle)

    return data

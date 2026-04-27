# Find differential regulons from SLE public data
# Sofia Salazar
# --------------------
'''
Usage:

python3 dif_regulons.py --outdir "/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/regulon/results" \
--loompy_file "/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/multi_runs_regulons_auc_trk.loom" \
--metadata "/Users/sofiasalazar/Desktop/LAB/meta-analysis-SLE/combined/data/all_data.csv"
'''

import os
import re
import math
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu as mwu
from statsmodels.stats.multitest import multipletests as multest
from scipy.stats import shapiro
from scipy.stats import ttest_ind
from itertools import compress
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--outdir')
parser.add_argument('--loompy_file')
parser.add_argument('--metadata')

args = parser.parse_args()

# ---------Functions---------

def test_normality(auc):
    all_values = auc.values.flatten()
    stat, p_value = shapiro(all_values)
    if p_value < 0.05:
        return False # "not normal" 
    else: return True # "normal"

def mwu_test(auc1, auc2, directional = False, adjust = True):
    # compute mann whitney u test non parametric
    if directional :
        regulonsPvalueVectorUp = []
        regulonsPvalueVectorDown = []
    else:
        regulonsPvalueVector = []
    
    regulons = []
    for regulon_num in range((len(auc1.columns) - 1)):
        regulons.append(auc1.columns[regulon_num]) # get regulon name

        auc1_vector = auc1.iloc[:,regulon_num]
        auc2_vector = auc2.iloc[:,regulon_num]
        if directional:
            res = mwu(auc1_vector, auc2_vector, method="asymptotic", alternative="greater")
            regulonsPvalueVectorUp.append(res[1])
            res = mwu(auc1_vector, auc2_vector, method="asymptotic", alternative="less")
            regulonsPvalueVectorDown.append(res[1])

        else:
            res = mwu(auc1_vector, auc2_vector, method="asymptotic")
            regulonsPvalueVector.append(res[1])


    if directional:
        if adjust:
            # Adjust p-values using Benjamini/Hochberg method
            regulonsPvalueVectorUp_Adj = multest(pvals = regulonsPvalueVectorUp,method = 'fdr_bh')[1]
            regulonsPvalueVectorDown_Adj = multest(pvals = regulonsPvalueVectorDown,method = 'fdr_bh')[1]
            out = pd.DataFrame({"regulon": regulons,
                            "upPvalues": regulonsPvalueVectorUp,
                            "AdjupPvalues":regulonsPvalueVectorUp_Adj,
                            "downPvalues":regulonsPvalueVectorDown,
                            "AdjdownPvalues": regulonsPvalueVectorDown_Adj})
        else:
            out = pd.DataFrame({"regulon": regulons,
                            "upPvalues": regulonsPvalueVectorUp,
                            "downPvalues":regulonsPvalueVectorDown })
    else:
        if adjust:
            # Adjust p-values using Benjamini/Hochberg method
            regulonsPvalueVector_Adj = multest(pvals = regulonsPvalueVector,method = 'fdr_bh')[1]
            out = pd.DataFrame({"regulon": regulons,
                            "Pvalues":regulonsPvalueVector,
                             "AdjPvalues": regulonsPvalueVector_Adj})
        else:
            out = pd.DataFrame({"regulon": regulons,
                            "Pvalues":regulonsPvalueVector })
    return out

def t_student_test(auc1, auc2, directional=False, adjust=True):
    # Initialize p-value vectors
    if directional:
        regulonsPvalueVectorUp = []
        regulonsPvalueVectorDown = []
    else:
        regulonsPvalueVector = []
    
    regulons = []
    for regulon_num in range((len(auc1.columns))):
        regulons.append(auc1.columns[regulon_num]) # Get regulon name

        auc1_vector = auc1.iloc[:, regulon_num]
        auc2_vector = auc2.iloc[:, regulon_num]
        
        if directional:
            # Perform one-tailed t-test for greater alternative
            res = ttest_ind(auc1_vector, auc2_vector, alternative='greater')
            regulonsPvalueVectorUp.append(res.pvalue)
            
            # Perform one-tailed t-test for less alternative
            res = ttest_ind(auc1_vector, auc2_vector, alternative='less')
            regulonsPvalueVectorDown.append(res.pvalue)
        else:
            # Perform two-tailed t-test
            res = ttest_ind(auc1_vector, auc2_vector)
            regulonsPvalueVector.append(res.pvalue)

    if directional:
        if adjust:
            # Adjust p-values using Benjamini/Hochberg method
            regulonsPvalueVectorUp_Adj = multest(pvals = regulonsPvalueVectorUp,method = 'fdr_bh')[1]
            regulonsPvalueVectorDown_Adj = multest(pvals = regulonsPvalueVectorDown,method = 'fdr_bh')[1]
            out = pd.DataFrame({"regulon": regulons,
                            "upPvalues": regulonsPvalueVectorUp,
                            "AdjupPvalues":regulonsPvalueVectorUp_Adj,
                            "downPvalues":regulonsPvalueVectorDown,
                            "AdjdownPvalues": regulonsPvalueVectorDown_Adj})
        else:
            out = pd.DataFrame({"regulon": regulons,
                            "upPvalues": regulonsPvalueVectorUp,
                            "downPvalues":regulonsPvalueVectorDown })
    else:
        if adjust:
            # Adjust p-values using Benjamini/Hochberg method
            regulonsPvalueVector_Adj = multest(pvals = regulonsPvalueVector,method = 'fdr_bh')[1]
            out = pd.DataFrame({"regulon": regulons,
                            "Pvalues":regulonsPvalueVector,
                             "AdjPvalues": regulonsPvalueVector_Adj})
        else:
            out = pd.DataFrame({"regulon": regulons,
                            "Pvalues":regulonsPvalueVector })
    
    return out

def l2fc(auc1, auc2):
    # compute log2 foldchanges between two auc matrices
    all_fc = []
    regulons = []
    for regulon_num in range((len(auc1.columns) - 1)):
        regulons.append(auc1.columns[regulon_num]) # get regulon name
        auc1_vector = auc1.iloc[:,regulon_num]
        auc2_vector = auc2.iloc[:,regulon_num]
        mean1 = np.mean(auc1_vector) + 1
        mean2 = np.mean(auc2_vector) + 1 

        res = np.log2(mean1 / mean2) # log2FC for a regulon

        all_fc.append(res)
    out = pd.DataFrame({"regulon": regulons,
                        "log2FC": all_fc})
    return out

def get_df(pval_df, lfc_df):
    out = pd.merge(pval_df, lfc_df, on = "regulon", how = "inner")
    return out
# ---------Main---------

print("INFO --- Reading metadata...")
meta = pd.read_csv(args.metadata, sep = ',')
metadf_all = pd.DataFrame(meta).set_axis(['num','ID','disease','study', 'batch'], axis = 1)

metadf = metadf_all[["ID", "disease"]]

print("INFO --- Splitting metadata into SLE and control samples")
meta_sle = metadf[metadf['disease'] == 'SLE']
meta_ctrl = metadf[metadf['disease'] == 'CONTROL']

print(f"INFO --- Metadata SLE shape: {meta_sle.shape}, metadata CONTROL shape: {meta_ctrl.shape}")

print("INFO --- Reading loompy file...")
lf = lp.connect(args.loompy_file, mode='r+', validate=False)
print("INFO --- Retrieving AUC matrix")

auc_mtx_multi = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
auc_mtx_multi.to_csv(os.path.join(args.outdir, "AUC_mtx.csv"), sep=",", index=True)
print("INFO --- Saved AUC matrix")


print(f"INFO --- AUC matrix shape: {auc_mtx_multi.shape}")

regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
    regulons[i] =  list(r[r==1].index.values)

print(f"INFO --- Total regulons found: {len(regulons.keys())}")

tf_targets = pd.DataFrame(lf.ra.Regulons, columns = list(auc_mtx_multi.columns), index= list(lf.ra.Gene))
tf_targets.to_csv(os.path.join(args.outdir, "tf_targets.csv"), sep=",", index=True)
print("INFO --- Saved tf targets data")
print("INFO --- Splitting AUC matrix into SLE and control samples")

auc_ctrl = auc_mtx_multi[auc_mtx_multi.index.isin(meta_ctrl['ID'])]
auc_sle = auc_mtx_multi[auc_mtx_multi.index.isin(meta_sle['ID'])]

print(f"INFO --- AUC matrix SLE shape: {auc_sle.shape}, AUC matrix CONTROL shape: {auc_ctrl.shape}")

# Sort regulons alphabetically
auc_ctrl = auc_ctrl.reindex(sorted(auc_ctrl.columns), axis=1)
auc_sle = auc_sle.reindex(sorted(auc_sle.columns), axis=1)

print("INFO --- Testing for normality of AUC values")

normal_flag = test_normality(auc_mtx_multi)
if normal_flag:
    print("INFO --- Data is normal, will conduct t-student test...")
    # non-directional result
    res_nd = t_student_test(auc_sle, auc_ctrl, directional=False, adjust=True)

else:
    print("INFO --- Data is not normal, will conduct Mann-Whitney U test")
    res_nd = mwu_test(auc_sle, auc_ctrl, directional=False, adjust=True)

n = (res_nd['Pvalues'] < 0.05).sum()
print(f"INFO --- Number of significantly different regulons: {n}")

plt.figure(figsize=(10, 6))
plt.hist(-np.log10(res_nd['AdjPvalues']), color="lightblue", edgecolor="lightblue")
plt.axvline(-np.log10(0.05), color='red', linestyle="dashed", linewidth=1)
plt.title('SLE vs Control differential regulons')
plt.ylabel('Frequency')
plt.xlabel('-log10(p-value)')
plt.savefig(os.path.join(args.outdir, 'histogram_SLE_Ctrl_regulons.png'), dpi=300)

print("INFO --- Computing log2Fold Changes...")
res_fc = l2fc(auc_sle, auc_ctrl)

print("INFO -- Assemblying output...")
res_df = get_df(pval_df=res_nd, lfc_df=res_fc)
print(f"INFO --- Output shape: {res_df.shape}, columns: {res_df.columns}")
print("INFO -- saving output")
res_df.to_csv(os.path.join(args.outdir, "difregs_SLE_Ctrl.csv"), sep=",", index=False)
print("FINISHED!")
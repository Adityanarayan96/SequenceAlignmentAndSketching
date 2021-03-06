from os.path import join
import numpy as np
import importlib
import pandas as pd
from sklearn.metrics import roc_curve,roc_auc_score
import matplotlib.pyplot as plt

PROJ_DIR = '/home/gcgreen2/alignment'
# GT_DIR = join(PROJ_DIR, 'spectral_jaccard_similarity/groundTruths')

# def load_dfs(dataset, pred_dir, gt_dir=GT_DIR): # load the alignment files
#     gt_path = join(gt_dir, dataset+'_daligner_ground_truth.txt')
#     pred_path = join(pred_dir, dataset+'_aln.tsv')
#     gt_df = pd.read_csv(gt_path, sep='\t', header=None, names=['i1','i2','overlap','l1','l2'])
#     pred_df = pd.read_csv(pred_path, sep='\t', header=None, names=['i1','i2','overlap','l1','l2'])
#     return gt_df, pred_df
# def load_dfs(pred_path, gt_path): # load the alignment files
#     pred_df = pd.read_csv(pred_path, sep='\t', header=None, names=['i1','i2','overlap','l1','l2'])
#     gt_df = pd.read_csv(gt_path, sep='\t', header=None, names=['i1','i2','overlap','l1','l2'])
#     return pred_df, gt_df

# def get_overlaps(pred_df, gt_df): 
#     '''returns [pred, gt]'''
#     overlaps = {}
#     for i in range(len(pred_df)):
#         line = pred_df.iloc[i]
#         i1, i2, pred_overlap = line.i1, line.i2, line.overlap
#         if i1 in overlaps:
#             overlaps[i1][i2] = (pred_overlap,0)
#         else:
#             overlaps[i1] = {i2:(pred_overlap,0)}

#     for i in range(len(gt_df)):
#         line = gt_df.iloc[i]
#         i1, i2, gt_overlap = line.i1, line.i2, line.overlap
#         if i1 in overlaps:
#             pred_overlap = overlaps[i1][i2][0] if i2 in overlaps[i1] else 0
#             overlaps[i1][i2] = (pred_overlap, gt_overlap)
#         else:
#             overlaps[i1] = {i2:(0,gt_overlap)}
#     overlaps = np.concatenate([list(overlaps[i1].values()) for i1 in overlaps])
#     return overlaps

def get_overlaps(pred_df, gt_df): 
    '''returns [pred, gt]'''
    n_seq = max(max(pred_df.i1),max(pred_df.i2),max(gt_df.i1),max(gt_df.i2))
    overlaps = np.full((n_seq,n_seq,2), 0)
    for i in range(len(pred_df)):
        line = pred_df.iloc[i]
        i1, i2, pred_overlap = int(line.i1-1), int(line.i2-1), line.overlap
        if i2<i1: continue
        overlaps[i1,i2,0] = pred_overlap

    for i in range(len(gt_df)):
        line = gt_df.iloc[i]
        i1, i2, gt_overlap = int(line.i1-1), int(line.i2-1), line.overlap
        if i2<i1: continue
        overlaps[i1,i2,1] = gt_overlap
    
    overlaps = overlaps.reshape(n_seq**2, 2)
    overlaps = overlaps[np.where(np.logical_or(overlaps[:,0]!=0, overlaps[:,1]!=0))]
    return overlaps

def roc(pred_dfs, gt_df, names, title):
    plt.figure()
    for pred_df,name in zip(pred_dfs,names):
        overlaps = get_overlaps(pred_df, gt_df)
        gt = overlaps[:,1] > 0 
        pred = overlaps[:,0]
        fpr,tpr,thresh = roc_curve(gt, pred)
        auroc = round(roc_auc_score(gt, pred),4)

        plt.plot(
            fpr,
            tpr,
            lw=2,
            label="%s, AUROC = %0.2f" % (name, auroc),
        )
    
    plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(title)
    plt.legend(loc="lower right")
    plt.show()
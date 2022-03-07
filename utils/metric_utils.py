from os.path import join
import numpy as np
import importlib
import pandas as pd
from sklearn.metrics import roc_curve,roc_auc_score
import matplotlib.pyplot as plt

PROJ_DIR = '/home/gcgreen2/alignment'
GT_DIR = join(PROJ_DIR, 'spectral_jaccard_similarity/groundTruths')

def load_dfs(dataset, pred_dir, gt_dir=GT_DIR): # load the alignment files
    gt_path = join(gt_dir, dataset+'_daligner_ground_truth.txt')
    pred_path = join(pred_dir, dataset+'_aln.tsv')
    gt_df = pd.read_csv(gt_path, sep='\t', header=None, names=['i1','i2','overlap','l1','l2'])
    pred_df = pd.read_csv(pred_path, sep='\t', header=None, names=['i1','i2','overlap','l1','l2'])
    return gt_df, pred_df

def get_overlaps(gt_df, pred_df):
    overlaps = {}
    for i in range(len(gt_df)):
        line = gt_df.iloc[i]
        i1, i2, overlap = line.i1, line.i2, line.overlap
        if i1 in overlaps:
            overlaps[i1][i2] = (overlap,0)
        else:
            overlaps[i1] = {i2:(overlap,0)}

    for i in range(len(pred_df)):
        line = pred_df.iloc[i]
        i1, i2, overlap = line.i1, line.i2, line.overlap
        if i1 in overlaps:
            gt_overlap = overlaps[i1][i2][0] if i2 in overlaps[i1] else 0
            overlaps[i1][i2] = (gt_overlap, overlap)
        else:
            overlaps[i1] = {i2:(0,overlap)}

    overlaps = np.concatenate([list(overlaps[i1].values()) for i1 in overlaps])
    return overlaps


def roc(dataset, pred_dir, runs, gt_dir=GT_DIR):
    plt.figure()
    lw=2
    
    for run in runs:
        run_dir = os.path.join(pred_dir, run)
        gt_df, pred_df = load_dfs(dataset, run_dir, gt_dir)
        overlaps = get_overlaps(gt_df, pred_df)
        gt = overlaps[:,0] > 0 
        pred = overlaps[:,1]
        fpr,tpr,thresh = roc_curve(gt, pred)
        auroc = round(roc_auc_score(gt, pred),4)

        plt.plot(
            fpr,
            tpr,
            lw=lw,
            label="%s, AUROC = %0.2f" % (run, auroc),
        )
    
    plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC for {dataset}")
    plt.legend(loc="lower right")
    plt.show()
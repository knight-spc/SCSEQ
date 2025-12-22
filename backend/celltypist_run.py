import argparse
import pandas as pd
import numpy as np
import sys
import os
import scanpy as sc
import celltypist
from celltypist import models

from check_file import check_files, copy_file_with_number

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--work_path', type=str, required=True)
    parser.add_argument('--annotation_model', type=str, required=True)
    parser.add_argument('--umap_data_path', type=str, required=True)
    args = parser.parse_args()
    print("work_path:", args.work_path)
    print("annotation_model:", args.annotation_model)
    print("umap_data_path:", args.umap_data_path)


    # 加载模型
    if f'{args.annotation_model}.pkl' in list(models.models_description().model):
        model = models.Model.load(model = f'{args.annotation_model}.pkl')
    else:
        model = models.Model.load(model = f'{args.work_path}/reference_seurat_data/pkl/{args.annotation_model}.pkl')
    
    # load data
    h5ad_path = os.path.join(args.work_path,'seurat_data.h5ad')
    adata = sc.read_h5ad(h5ad_path)
    predictions_new = celltypist.annotate(adata, model = model, majority_voting = True, p_thres=0.5)
    adata_new = predictions_new.to_adata()

    umap_data = pd.read_csv(args.umap_data_path)
    umap_data["annotation"] = list(adata_new.obs["majority_voting"])
    umap_data["annotation"] = umap_data["annotation"].str.replace(' ', '_')

    umap_data.to_csv(args.umap_data_path, index=False)
    print("注释完成")
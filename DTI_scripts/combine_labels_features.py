#coding=utf-8
#python3.8
#used for combine labels and all features

import numpy as np
import pandas as pd
import pickle

#read lables
labels_df = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/data_labels_all.csv", header=0)
labels_df.head()

#combine pocket features
pocket_physicochemical_features = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/pockets_physicochemical_property.csv", header=0)
pocket_geometric_features = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/pockets_geometric_features.csv", header=0)
pocket_dynamic_features = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/pockets_dynamic_property.csv", header=0)
pocket_features = pd.concat([pocket_physicochemical_features, pocket_geometric_features.iloc[:,2:],pocket_dynamic_features.iloc[:,2:]], axis=1)

#load ligand smiles
with open("/home/databank/zlp/md-40/DTI_prediction/ligands_descriptors.pkl",'rb') as f:
    ligand_fps = pickle.load(f)

ligand_fps_names = ["erg", "ecfp4", "ecfp6", "rdkit2d", "daylight", "pubchem", "maccs"]
all_data = pd.DataFrame()

for index in labels_df.index.values:
    data = pd.DataFrame([labels_df.iloc[index,:].to_dict()], index=[index])

    file = labels_df.loc[index, "File"]
    pocket_name = labels_df.loc[index, "pocket name"]
    lig_smi = labels_df.loc[index, "ligand smiles"]

    pocket_feature = pocket_features[(pocket_features["File"] == file) & (pocket_features["pocket name"]==pocket_name)].iloc[:,2:].copy()
    #change index for concat
    pocket_feature.index = [index]

    #combile pocket features and lables
    data_pocket_feature = pd.concat([data, pocket_feature], axis=1)
    all_data = pd.concat([all_data, data_pocket_feature], axis=0, ignore_index=True)
    
    #add ligand fps
    for fp in ligand_fps_names:
        lig_fp = "_".join([str(i) for i in list(ligand_fps[lig_smi][fp])])
        all_data.loc[index, "ligand_"+fp] = lig_fp
    
    #print(all_data)
all_data.to_csv("all_labels_features.csv", index=False, header=True)

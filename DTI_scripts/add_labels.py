#coding=utf-8
#python3.8
#标记label

import pandas as pd
from pymol import cmd
import numpy as np
from glob import glob
import os
from rdkit import Chem


def convert_sdf2smi(lig):
    mol = Chem.MolFromMolFile(lig)
    canonical_smi = Chem.MolToSmiles(mol)
    return canonical_smi


data_df = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/PDBbind_drugbank_data_2.csv", header=0)
path = "/home/databank/zlp/md-40/"
lig_path = "/home/databank/zlp/md-40/DTI_prediction/ligands/"
f = open("data_labels.csv","w")
f.write("File,pocket name,pdb id,ligand name,ligand smiles,labels"+"\n")

for index in data_df.index.values:
    file = data_df.loc[index, "File"]
    lig_occopied_poc = data_df.loc[index, "occupied poc"]

    pdb_id = data_df.loc[index, 'PDB ID']
    ligand_name = data_df.loc[index, 'lig name']
    lig_sdf = glob(lig_path+pdb_id.lower()+"_[A-Z]_"+ligand_name+".sdf")[0]
    smi = convert_sdf2smi(lig_sdf)

    pse = path+file+"/workspace/PocketStablility.pse"
    cmd.delete("all")
    cmd.load(pse)
    poc_names = sorted([i for i in cmd.get_names("all") if i[0] == "P"])
    for poc in poc_names:
        if poc == lig_occopied_poc:
            label = 1
        else:
            label = 0
        line = [file,poc,pdb_id,ligand_name,smi,str(label)]
        f.write(",".join(line))
        f.write("\n")

f.close()



    
    

    


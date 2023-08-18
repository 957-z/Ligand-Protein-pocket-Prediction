#coding=utf-8
#python3.8
#used for physicochemical property calculation for protein pocket

import numpy as np
import pandas as pd
from pymol import cmd

def get_chain_resi_resn(sele):
    space = {"res_info":[]}
    cmd.iterate(sele,"res_info.append([chain,resi,resn])", space=space)
    res_info = space['res_info']
    #remove duplicates
    res_info = [list(i) for i in set(tuple(_) for _ in res_info)]
    return res_info



#reference: R. Krivák and D. Hoksza. P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. Journal of Cheminformatics. 10, 39,(2018).
#modified artificially
# HYDROPHOBIC = ['PHE','ILE','TRP','GLY','LEU','VAL','MET','ALA', 'PRO']
# HYDROPHILIC = ['ARG','ASN','ASP','GLN', 'GLU','LYS','CYS','TYR', "SER", "THR","HIS"]
HYDROPATHY_INDEX = {'ALA':1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,'CYS':2.5,'GLU':-3.5, 'GLN':-3.5,'GLY':-0.4, 'HIS':-3.2, 'ILE':4.5,'LEU':3.8,'LYS':-3.9, 'MET':1.9, 'PHE':2.8,'PRO':-1.6,'SER':-0.8, 'THR':-0.7,'TRP':-0.9, 'TYR':-1.3, 'VAL':4.2}
ALIPHATIC = ['ALA', 'LEU', 'ILE', 'VAL', 'GLY', 'PRO']
AROMATIC = ['PHE', 'TRP', 'TYR',"HIS"]
SULFUR = ['CYS', 'MET']
HYDROXYL = ['SER', 'THR']
BASIC = ['ARG', 'LYS', 'HIS']
ACIDIC = ['ASP', 'GLU']
AMIDE = ['ASN', 'GLN']
#CHARGE = {'ASP':-1,'GLU':-1,'ARG':1,'HIS':1,'LYS':1}
POLAR = ['ARG','ASN','ASP','GLN','GLU','HIS','LYS','SER','THR','TYR','CYS']
IONIZABLE = ['ASP','GLU','HIS','LYS','ARG','CYS','TYR']
#only sidechains are considered
HB_DONOR = ['ARG', 'LYS', 'TYR']
HB_ACCEPTOR = ['ASP', 'GLU']
HB_DONOR_ACCEPTOR = ['ASN', 'GLN', 'HIS', 'SER', 'THR', 'TYR']

properties = ["hydrophobicity", "aliphatic", "aromatic", "sulfur", "hydroxyl", "basic", "acidic", "amide", "polar", "ionizable", "hb_donor", "hb_acceptor", "hb_donor_acceptor"]

data_df = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/data_labels.csv", header=0)
path = "/home/databank/zlp/md-40/"
output = open("pockets_physicochemical_property.csv","w")
output.write("File,pocket name,"+",".join(properties)+"\n")

for key,val in data_df.groupby("File"):
    pse = path+key+"/workspace/PocketStablility.pse"
    cmd.delete("all")
    cmd.load(pse)
    poc_names = sorted([i for i in cmd.get_names("all") if i[0] == "P"])
    for poc in poc_names:
        cmd.select("poc_surf","byres (polymer.protein within 3 of "+poc+")")
        poc_chain_resi_res = get_chain_resi_resn("poc_surf")
        total_num = len(poc_chain_resi_res)

        for property in properties:
            vars()[property] = 0
            
        for res in poc_chain_resi_res:
            resn = res[-1].upper()
            if resn == "ASH":
                resn = "ASP"
            for property in properties:
                if property == "hydrophobicity":
                    vars()[property] += HYDROPATHY_INDEX[resn]
                else:
                    if resn in vars()[property.upper()]:
                        vars()[property] += 1
        line = [key,poc]
        for property in properties:
            if property == "hydrophobicity":
                line.append(str(np.round(vars()[property]/total_num,2)))
            else:
                line.append(str(np.round((vars()[property])/total_num,2)))
        output.write(",".join(line))
        output.write("\n")

output.close()


            







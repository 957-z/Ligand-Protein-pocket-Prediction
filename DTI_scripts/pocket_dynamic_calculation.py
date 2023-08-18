#coding=utf-8
#python3.8
#used for dynamic property calculation for protein pocket

import numpy as np
import pandas as pd
from pymol import cmd

def get_stability_mean_std(sele):
    space = {"atm_info":[]}
    cmd.iterate(sele,"atm_info.append(b)", space=space)
    stability = space['atm_info']
    
    return np.round(np.mean(stability),2), np.round(np.std(stability),2)

data_df = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/data_labels.csv", header=0)
path = "/home/databank/zlp/md-40/"
output = open("pockets_dynamic_property.csv","w")
output.write("File,pocket name,stability_mean,stability_std,stability_3A_mean,stability_3A_std"+"\n")

for key,val in data_df.groupby("File"):
    pse = path+key+"/workspace/PocketStablility.pse"
    cmd.delete("all")
    cmd.load(pse)
    poc_names = sorted([i for i in cmd.get_names("all") if i[0] == "P"])
    for poc in poc_names:
        cmd.select("poc_surf",poc + " within 3 of polymer.protein")
        stability_mean, stability_std = get_stability_mean_std(poc)
        stability_3A_mean, stability_3A_std = get_stability_mean_std("poc_surf")
        line = [key, poc, str(stability_mean), str(stability_std), str(stability_3A_mean), str(stability_3A_std)]
        output.write(",".join(line))
        output.write("\n")

output.close()    


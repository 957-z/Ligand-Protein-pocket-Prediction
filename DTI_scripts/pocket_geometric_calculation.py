#coding=utf-8
#python3.8

import numpy as np
import pandas as pd
from pymol import cmd

def get_resi_resn(sele):
    space = {"res_info":[]}
    cmd.iterate(sele,"res_info.append([resi,resn])", space=space)
    res_info = space['res_info']
    #remove duplicates
    res_info = [list(i) for i in set(tuple(_) for _ in res_info)]
    return res_info

data_df = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/data_labels.csv", header=0)
path = "/home/databank/zlp/md-40/DTI_prediction/grid_labeled_pse_1.9/"
output = open("pockets_geometric_features.csv","w")
output.write("File,pocket name,volume,depth,sa,enclosure"+"\n")

for key,val in data_df.groupby("File"):
#for key in ["era"]:
    pse = path+key+".pse"
    cmd.delete("all")
    cmd.load(pse)
    poc_names = sorted([i for i in cmd.get_names("all") if len(i) == 3])
    
    for poc in poc_names:
        volume = cmd.count_atoms(poc)
        surface_grids = get_resi_resn(poc+"_surface_grids")
        #surface area is defined as The number of surface grids multiplied by the squared grid spacing (1 A)
        sa = len(surface_grids)
        
        try:
            exposed_grids = get_resi_resn(poc+"_exposed_grids")
            #enclosure is defined as the proportion of the number of surface grids to the number of hull grids
            enclosure = sa/(len(surface_grids)+len(exposed_grids))
            surface_grids_depth = {}
            for i in surface_grids:
                grid1 = "(resi "+str(i[0])+" and resn "+ i[1]+")"
                all_distances = []
                for j in exposed_grids:
                    grid2 = "(resi "+str(j[0])+" and resn "+ j[1]+")"
                    distance = cmd.distance("tmp",grid1, grid2)
                    cmd.delete("tmp")
                    all_distances.append(distance)
                surface_grids_depth["_".join(i)] = min(all_distances)
                #add property label to each grid, and user can color them in pse file accodring to different depths
                cmd.set_atom_property("depth",min(all_distances),grid1,proptype=3)
            depth = max(surface_grids_depth.values())
        
        #in some systems, there is no exposed grids. the depth is define the maximum distance between surface grids.
        except Exception as e:
            print(e)
            all_distances = []
            for i in surface_grids:
                grid1 = "(resi "+str(i[0])+" and resn "+ i[1]+")"
                for j in surface_grids:
                    grid2 = "(resi "+str(j[0])+" and resn "+ j[1]+")"
                    distance = cmd.distance("tmp",grid1, grid2)
                    cmd.delete("tmp")
                    all_distances.append(distance)
            depth = max(all_distances)
            enclosure = 1
        
        line = [key,poc,str(np.round(volume,2)), str(np.round(depth,2)), str(np.round(sa,2)), str(np.round(enclosure,2))]
        output.write(",".join(line))
        output.write("\n")
    
    cmd.save("/home/databank/zlp/md-40/DTI_prediction/grid_labeled_depth/"+key+"_depth.pse")

output.close()



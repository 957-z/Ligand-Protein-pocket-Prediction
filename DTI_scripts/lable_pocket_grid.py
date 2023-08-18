#coding=utf-8
#python3.8
#used for classification to three classes

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

#vdw raidus + 1.6 (the radius of small shpere used to pocket detection)+0.3(for tolerance)
def get_neighbor_protein_atm(sele):
    elements_vdw = {"C":1.70,'N':1.55,'O':1.52,"H":1.20}
    neighbor_protein_num = 0
    for ele in elements_vdw:
        cutoff = elements_vdw[ele]+1.9
        neighbor_protein_num += cmd.count_atoms("(polymer.protein and e. "+ele+") within "+str(cutoff)+" of "+sele)
    return neighbor_protein_num

def color_grids(poc,surface_grids, exposed_grids, inside_grids):
    try:
        cmd.select(poc+"_inside_grids"," or ".join(["(resi "+str(grid[0])+" and resn "+grid[1]+")" for grid in inside_grids]))
        cmd.color("gray90",poc+"_inside_grids")
    except Exception as e:
        print(e)
        print(poc,len(exposed_grids))

    cmd.select(poc+"_surface_grids"," or ".join(["(resi "+str(grid[0])+" and resn "+grid[1]+")" for grid in surface_grids]))
    cmd.color("blue",poc+"_surface_grids")
    
    try:
        cmd.select(poc+"_exposed_grids"," or ".join(["(resi "+str(grid[0])+" and resn "+grid[1]+")" for grid in exposed_grids]))
        cmd.color("red",poc+"_exposed_grids")
    except Exception as e:
        print(e)
        print(poc,len(exposed_grids))

data_df = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/data_labels.csv", header=0)
path = "/home/databank/zlp/md-40/"

for key,val in data_df.groupby("File"):
    pse = path+key+"/workspace/PocketStablility.pse"
    cmd.delete("all")
    cmd.load(pse)
    poc_names = sorted([i for i in cmd.get_names("all") if i[0] == "P"])
    
    for poc in poc_names:
        surface_grids = []
        exposed_grids = []
        poc_grids = get_resi_resn(poc)
        inside_grids = poc_grids[:]
        for grid in poc_grids:
            grid_ = "(resi "+str(grid[0])+" and resn "+grid[1]+")"
            neighbor_grids_num = cmd.count_atoms(poc+" within 1.2 of "+grid_)
            if neighbor_grids_num < 7:
                inside_grids.remove(grid)
                neighbor_protein_num = get_neighbor_protein_atm(grid_)
                if neighbor_protein_num > 0:
                    surface_grids.append(grid)
                else:
                    exposed_grids.append(grid)

        color_grids(poc, surface_grids, exposed_grids, inside_grids)
        #cmd.save("test-1.pse")
        #regroup for exposed_grids (some exposed girds are classified into wrong class)
        # exposed_grids_new = exposed_grids[:]
        # surface_grids_new = surface_grids[:]
        # for exposed_grid in exposed_grids:
        #     exposed_grid_ = "(resi "+str(exposed_grid[0])+" and resn "+exposed_grid[1]+")"
        #     #neighbor_exposed_grids_num = cmd.count_atoms(poc+"_exposed_grids within 2 of "+exposed_grid_)
        #     neighbor_surface_grids_num = cmd.count_atoms(poc+"_surface_grids within 1.2 of "+exposed_grid_)
        #     #neighbor_inside_grids_num = cmd.count_atoms(poc+"_inside_grids within 2 of "+exposed_grid_)
        #     if exposed_grid[0] == 484:
        #         print(neighbor_surface_grids_num)
        #     if neighbor_surface_grids_num >= 3:
        #         print(exposed_grid_)
        #         exposed_grids_new.remove(exposed_grid)
        #         surface_grids_new.append(exposed_grid)
        # color_grids(poc, surface_grids_new, exposed_grids_new, inside_grids)
        # cmd.save("test-2.pse")
        
    cmd.save("/home/databank/zlp/md-40/DTI_prediction/grid_labeled_pse_1.9/"+key+".pse")
    



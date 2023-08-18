#coding=utf-8
#python3.8
#used for calculation fingerprint for ligands

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.rdReducedGraphs import GetErGFingerprint
from descriptastorus.descriptors import rdDescriptors, rdNormalizedDescriptors
from pybiomed_helper import calcPubChemFingerAll
from rdkit.Chem import MACCSkeys

def smiles2erg(s):
    try:
        mol = Chem.MolFromSmiles(s)
        features = np.array(GetErGFingerprint(mol))
    except:
        print('rdkit cannot find this SMILES for ErG: ' + s + 'convert to all 0 features')
        features = np.zeros((315, ))
    return features

def smiles2morgan(s, radius = 2, nBits = 1024):
    try:
        mol = Chem.MolFromSmiles(s)
        features_vec = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
        features = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(features_vec, features)
    except:
        print('rdkit not found this smiles for morgan: ' + s + ' convert to all 0 features')
        features = np.zeros((nBits, ))
    return features

def smiles2rdkit2d(s):    
    try:
        generator = rdNormalizedDescriptors.RDKit2DNormalized()
        features = np.array(generator.process(s)[1:])
        NaNs = np.isnan(features)
        features[NaNs] = 0
    except:
        print('descriptastorus not found this smiles: ' + s + ' convert to all 0 features')
        features = np.zeros((200, ))
    return np.array(features)

def smiles2daylight(s):
	try:
		NumFinger = 2048
		mol = Chem.MolFromSmiles(s)
		bv = FingerprintMols.FingerprintMol(mol)
		temp = tuple(bv.GetOnBits())
		features = np.zeros((NumFinger, ))
		features[np.array(temp)] = 1
	except:
		print('rdkit not found this smiles: ' + s + ' convert to all 0 features')
		features = np.zeros((2048, ))
	return np.array(features)

def smiles2pubchem(s):
	try:
		features = calcPubChemFingerAll(s)
	except:
		print('pubchem fingerprint not working for smiles: ' + s + ' convert to 0 vectors')
		features = np.zeros((881, ))
	return np.array(features)

def smiles2maccs(s):
    try:
        mol = Chem.MolFromSmiles(s)
        fp = MACCSkeys.GenMACCSKeys(mol)
        fp_arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp,fp_arr)
        #print(np.array(fp_arr))
    except:
        print('maccs fingerprint not working for smiles: ' + s + ' convert to 0 vectors')
        fp_arr = np.zeros((166, ))
    #only 1-166 bits represent maccs
    return fp_arr[1:]

from collections import defaultdict
import pickle

data_df = pd.read_csv("/home/databank/zlp/md-40/DTI_prediction/data_labels.csv", header=0)
all_fps = defaultdict(dict)

for key,val in data_df.groupby("ligand smiles"):
      erg = smiles2erg(key)
      all_fps[key]["erg"] = erg
      ecfp4 = smiles2morgan(key, radius=2, nBits=1024)
      all_fps[key]["ecfp4"] = ecfp4
      ecfp6 = smiles2morgan(key, radius=3, nBits=1024)
      all_fps[key]["ecfp6"] = ecfp6
      rdkit2d = smiles2rdkit2d(key)
      all_fps[key]["rdkit2d"] = rdkit2d
      daylight = smiles2daylight(key)
      all_fps[key]["daylight"] = daylight
      pubchem = smiles2pubchem(key)
      all_fps[key]["pubchem"] = pubchem
      maccs = smiles2maccs(key)
      all_fps[key]["maccs"] = maccs

#save to pikle file
with open('ligands_descriptors.pkl', 'wb') as f:
    pickle.dump(all_fps, f)


#using ecfp4 fingerprint
import multiprocessing
from itertools import product
from tqdm import tqdm
from sklearn.model_selection import KFold
from collections import defaultdict
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, recall_score, f1_score,roc_auc_score, average_precision_score
from imblearn.over_sampling import RandomOverSampler, SMOTE
from xgboost import XGBClassifier
from sklearn import preprocessing
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore")

def oversampling(features, df): 
    X = df[features]
    y = df["labels"]
    sampling = SMOTE(random_state=111)
    X_, y_ = sampling.fit_resample(X, y)
    return X_, y_

def get_k_folds(data, k):
    kf = KFold(n_splits=k, random_state=None, shuffle=False)
    kfolds = {}
    for i, (train_index, test_index) in enumerate(kf.split(data)):
        kfolds[i] = (train_index, test_index)
    return kfolds

def train_val_split(train_index, val_index, df_train_val):
    df_train_systems = np.unique(df_train_val["File"])[train_index].tolist()
    df_train = df_train_val[df_train_val["File"].str.contains("|".join(df_train_systems))]
    
    df_val_systems = np.unique(df_train_val["File"])[val_index].tolist()
    df_val = df_train_val[df_train_val["File"].str.contains("|".join(df_val_systems))]
    #print(np.unique(df_train["File"]), np.unique(df_val["File"]))
    
    return df_train, df_val

def train_eval_xgboost_model(parameter):
    features = ['hydrophobicity', 'aliphatic', 'aromatic', 'sulfur',
        'hydroxyl', 'basic', 'acidic', 'amide', 'polar', 'ionizable',
        'hb_donor', 'hb_acceptor', 'hb_donor_acceptor', 'volume', 'depth', 'sa',
        'enclosure', 'stability_mean', 'stability_std'] + list(ligand_ecfp4_fp.columns)
    
    max_depth, n_estimators, learning_rate, subsample, gamma, min_child_weight = parameter
    kfolds = get_k_folds(np.unique(df_train_val["File"]), 10)
    
    
    metrics_val = defaultdict(list)
    for k,index in kfolds.items():
        #print(index)
        train_index = index[0]
        val_index = index[1]
        df_train, df_val = train_val_split(train_index, val_index, df_train_val)
        X_train, y_train = oversampling(features, df_train)
        X_val, y_val = oversampling(features, df_val)
        
        model = XGBClassifier(max_depth=max_depth, n_estimators=n_estimators, learning_rate=learning_rate, subsample=subsample, 
                            gamma=gamma, min_child_weight=min_child_weight, 
                            objective='binary:logistic', seed=0, eval_metric="logloss",tree_method='gpu_hist', gpu_id=0)
        model.fit(X_train, y_train)
        
        y_val_pred_proba = model.predict_proba(X_val)
        ap = np.round(average_precision_score(y_val, y_val_pred_proba[:,1]),2)
        roc_auc = np.round(roc_auc_score(y_val, y_val_pred_proba[:,1]),2)
        
        y_val_pred = model.predict(X_val)
        f1 = np.round(f1_score(y_val, y_val_pred),2)
        recall = np.round(recall_score(y_val, y_val_pred),2)
        
        metrics_val["ap"].append(ap)
        metrics_val["roc_auc"].append(roc_auc)
        metrics_val["f1"].append(f1)
        metrics_val["recall"].append(recall)
    #print(np.mean(pd.DataFrame(metrics_val), axis=0))
    matrics_val_mean = [str(i) for i in parameter]+np.round(np.mean(pd.DataFrame(metrics_val), axis=0), 2).astype(str).tolist()
    print("done!")
    return matrics_val_mean
    
def write_output(line):
    output.write(",".join(line))
    output.write("\n")

max_depth = np.arange(2,9,2)
n_estimators = np.arange(20,101,10)
learning_rate = [0.1]
subsample = np.arange(0.4,1.05,0.1)
gamma = np.arange(0,1.1,0.2)
min_child_weight = [1,3,5,7]
# max_depth = [3]
# n_estimators = [5]
# learning_rate = [0.1]
# subsample = [0.5]
# gamma = [0.1]
# min_child_weight = [3,6]

paras = []
for max_depth, n_estimators, learning_rate, subsample, gamma, min_child_weight in product(max_depth, n_estimators, learning_rate, subsample, gamma, min_child_weight):
    paras.append((max_depth, n_estimators, learning_rate, subsample, gamma, min_child_weight))

data = pd.read_csv("../all_labels_features.csv", header=0)
ligand_ecfp4_fp = data["ligand_ecfp4"].str.split("_", expand=True).astype(float)
ligand_ecfp4_fp.columns = ["ligand_ecfp4_"+ str(i) for i in range(ligand_ecfp4_fp.shape[1])]

df = pd.concat([data, ligand_ecfp4_fp], axis=1)

df_train_val = df.drop(index=df[(df["File"] == "1kvl") | (df["File"] == "1z95") | (df["File"] == "fkbp5")].index.values)  
output = open("xgboost_ecfp4_kfolds_gridsearch_results-gpu.csv","a+")
output.write("max_depth, n_estimators, learning_rate, subsample, gamma, min_child_weight,ap,roc_auc,f1,recall"+"\n")

#pool = multiprocessing.Pool(16)
import time
complete_num = 0
print(time.asctime())
for para in paras:
    metrics_mean = train_eval_xgboost_model(para)
    write_output(metrics_mean)
    complete_num += 1
    if complete_num % 100 == 0:
        print(complete_num,time.asctime())
#     pool.apply_async(func=train_eval_xgboost_model, args=(para,), callback=write_output)

# pool.close()
# pool.join()
output.close()

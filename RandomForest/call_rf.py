import joblib
import numpy as np
import pandas as pd
import time 
import argparse
import subprocess

import sys
import os
sys.path.append(os.path.join(os.getcwd(), '../utils'))
from features import *
from datautil import hard_filter_keeporg

SOM_SNV_FEATURES=len(som_selected_features) + 1
SOM_INDEL_FEATURES = SOM_SNV_FEATURES + 3

type_to_label = {"SNV": 0, "Deletion": 1, "Insertion": 2, "Complex": 3, "MNV":3}
varLabel_to_label = {
    "Germline":0, "StrongLOH":1, "LikelyLOH":2, "StrongSomatic":3,
    "LikelySomatic":4, "AFDiff":5, "SampleSpecific":6
}
fe2i = som_features_to_index
fvc_sf = som_selected_features
def format_snv_data_item(jri, fisher):
    if not fisher:
        print("not support to train if you not run rabbitvar without --fiser!!")
        exit(-1)
    data = list()
    # key is chrom:pos like "chr1:131022:A:T"
    key = jri[2] + ":" + jri[3] + ":" + jri[5] + ":" + jri[6] #TODO: case sensitive
    for sf in fvc_sf:
        data.append(jri[fe2i[sf]])
    data.append(varLabel_to_label[jri[fe2i["VarLabel"]]])#varlabel
    if len(data) != SOM_SNV_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        exit(-1)
    return key, data

def format_indel_data_item(jri, fisher):
    if not fisher:
        print("not support to train if you not run rabbitvar without --fiser!!")
        exit(-1)
    data = list()
    # key is chrom:pos like "chr1:131022:A:T"
    key = jri[2] + ":" + jri[3] + ":" + jri[5] + ":" + jri[6] #TODO: case sensitive
    data.append(len(jri[5]))
    data.append(len(jri[6]))
    data.append(type_to_label[jri[fe2i["VarType"]]])
    for sf in fvc_sf:
        data.append(float(jri[fe2i[sf]]))
    data.append(varLabel_to_label[jri[fe2i["VarLabel"]]])#varlabel
    #data.append(jri[fe2i["VarLabel"]])
    if len(data) != SOM_INDEL_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        exit(-1)
    return key, data

def my_predict(clf, data, scale):
    proba = clf.predict_proba(data)
    return np.asarray([1 if x > scale else 0 for x in proba[:,1]])

def call_rf(args):
    in_file = args.in_file
    vartype = args.var_type
    cr = list()
    raw = list()
    cr_start = time.time()
    with open(in_file, 'r') as f:
        for var in f:
            if var[0] == '#':
                continue
            items = var.strip().split('\t')
            if(len(items) == 61): #36 or 38(fisher)
                if vartype == "SNV" and items[fe2i["VarType"]] != "SNV":
                    continue
                elif vartype == "INDEL" and items[fe2i["VarType"]] == "SNV":
                    continue
                elif items[fe2i["VarLabel"]] == "Germline":
                    continue
                if vartype == 'SNV':
                    key, data = format_snv_data_item(items, True)
                elif vartype == 'INDEL':
                    key, data = format_indel_data_item(items, True)
                cr.append(data)
                raw.append(var)
            else:
                print("wrong data format!!", len(items))
    cr_end = time.time()
    print("time of load data to cr: {} s".format(cr_end - cr_start) )
    
    cr_start = time.time()
    clf = joblib.load(args.model)
    cr_end = time.time()
    print("time of load model: {} s".format(cr_end - cr_start) )

    cr_start = time.time()
    #data = np.asfarray(cr)
    data = pd.DataFrame(cr)
    data.columns = ["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel"]
    cr_end = time.time()
    print("time of np.asaray: {} s".format(cr_end - cr_start) )
    
    #--hardfilter--#
    print("before hard filter: {} data".format(len(data)))
    data = hard_filter_keeporg(data)
    print("after hard filter: {} data".format(len(data)))

    cr_start = time.time()
    #pred = clf.predict(data)
    scale = args.scale
    pred = my_predict(clf, data[data.columns[:-1]], scale)
    data['pred'] = pred
    cr_end = time.time()
    print("length {} - {}".format(len(data), len(pred)))
    print("time of pred: {} s".format(cr_end - cr_start) )

    cr_start = time.time()
    with open(args.out_file, 'w') as f:
        for i in range(len(pred)):
            if pred[i] == 1 and data['hard_flag'][i] == 0:
                f.write(str(i) + ":"+str(raw[i]))
    cr_end = time.time()
    print("time of write: {} s".format(cr_end - cr_start) )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "train your network")
    parser.add_argument('--in_file', help = "RabbitVar intermidiate file(with fisher test)", type=str, required = True)
    parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--scale', help = 'scale for random forest predict_proba to filter, 0.5 perform the same result as function predict. (default > 0.2)', type = float, default = 0.2)
    #parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--nthreads', help = "number of thread", type=int, default=20)
    parser.add_argument('--model', help = "random forest model", type=str, required = True)
    parser.add_argument('--out_file', help = "file to output", type=str, required = True)
    args = parser.parse_args()
    call_rf(args)

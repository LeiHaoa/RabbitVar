#!/usr/bin/env python3
import sys
from multiprocessing import Process, Queue

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
import os

som_features = [
    "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt",
    "Var1Depth", "Var1AltDepth", "Var1RefFwdReads", "Var1RefRevReads", "Var1AltFwdReads", "Var1AltRevReads", "Var1Genotype", "Var1AF",
    "Var1Bias", "Var1PMean", "Var1PStd", "Var1QMean", "Var1QStd", "Var1MQ", "Var1Sig_Noise", "Var1HiAF", "Var1ExtraAF", "Var1NM", "Var1Pvalue", "Var1Oddr",
    "Var2Depth", "Var2AltDepth", "Var2RefFwdReads", "Var2RefRevReads", "Var2AltFwdReads", "Var2AltRevReads", "Var2Genotype", "Var2AF",
    "Var2Bias", "Var2PMean", "Var2PStd", "Var2QMean", "Var2QStd", "Var2MQ", "Var2Sig_Noise", "Var2HiAF", "Var2ExtraAF", "Var2NM", "Var2Pvalue", "Var2Oddr",
    "shift3", "MSI", "MSI_NT", "5pFlankSeq", "3pFlankSeq", "Seg", "VarLabel", "VarType",
    "Duprate1", "SV_info1", "Duprate2", "SV_info2", "Pvalue", "Oddratio"
    ]
som_features_to_index = {
    "Sample":0, "Gene":1, "Chr":2, "Start":3, "End":4, "Ref":5, "Alt":6,
    "Var1Depth":7, "Var1AltDepth":8, "Var1RefFwdReads":9, "Var1RefRevReads":10, "Var1AltFwdReads":11, "Var1AltRevReads":12, "Var1Genotype":13, "Var1AF":14,
    "Var1Bias":15, "Var1PMean":16, "Var1PStd":17, "Var1QMean":18, "Var1QStd":19, "Var1MQ":20, "Var1Sig_Noise":21, "Var1HiAF":22, "Var1ExtraAF":23, "Var1NM":24, "Var1Pvalue":25, "Var1Oddr":26,
    "Var2Depth":27, "Var2AltDepth":28, "Var2RefFwdReads":29, "Var2RefRevReads":30, "Var2AltFwdReads":31, "Var2AltRevReads":32, "Var2Genotype":33, "Var2AF":34,
    "Var2Bias":35, "Var2PMean":36, "Var2PStd":37, "Var2QMean":38, "Var2QStd":39, "Var2MQ":40, "Var2Sig_Noise":41, "Var2HiAF":42, "Var2ExtraAF":43, "Var2NM":44, "Var2Pvalue":45, "Var2Oddr":46,
    "shift3":47, "MSI":48, "MSI_NT":49, "5pFlankSeq":50, "3pFlankSeq":51, "Seg":52, "VarLabel":53, "VarType":54,
    "Duprate1":55, "SV_info1":56, "Duprate2":57, "SV_info2":58, "Pvalue":59, "Oddratio":60
    }
som_selected_features = [
    "Var1Depth", "Var1AltDepth", "Var1RefFwdReads", "Var1RefRevReads", "Var1AltFwdReads", "Var1AltRevReads", "Var1AF",
    "Var1PMean", "Var1PStd", "Var1QMean", "Var1QStd", "Var1MQ", "Var1Sig_Noise", "Var1HiAF", "Var1ExtraAF", "Var1NM", "Var1Pvalue", "Var1Oddr",
    "Var2Depth", "Var2AltDepth", "Var2RefFwdReads", "Var2RefRevReads", "Var2AltFwdReads", "Var2AltRevReads", "Var2AF",
    "Var2PMean", "Var2PStd", "Var2QMean", "Var2QStd", "Var2MQ", "Var2Sig_Noise", "Var2HiAF", "Var2ExtraAF", "Var2NM", "Var2Pvalue", "Var2Oddr",
    "shift3", "MSI", "MSI_NT",
    "Pvalue", "Oddratio"
    ]
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
    data.append( str(varLabel_to_label[jri[fe2i["VarLabel"]]]) )#varlabel
    #data.append(jri[fe2i["VarLabel"]])
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
    data.append( str(len(jri[5])) )
    data.append( str(len(jri[6])) )
    data.append( str(type_to_label[jri[fe2i["VarType"]]]) )
    for sf in fvc_sf:
        data.append(jri[fe2i[sf]])
    data.append( str(varLabel_to_label[jri[fe2i["VarLabel"]]]) )#varlabel
    #data.append(jri[fe2i["VarLabel"]])
    if len(data) != SOM_INDEL_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        exit(-1)
    return key, data

def get_label(truth_vars, key):
    if key in truth_vars:
        return 1
    else :
        return 0
    return 0

import sys

if len(sys.argv) != 5:
    print('intput parameter number wrong: '.format(len(sys.argv)))
    exit(-1)

fastvc_file = sys.argv[1]
truth_file = sys.argv[2]
vartype = sys.argv[3]
tsv_file = sys.argv[4]

truth_vars = set()
with open(truth_file, 'r') as f:
    for var in f:
        if var[0] == '#':
            continue
        items = var.split('\t')
        chrom, pos, id, ref, alt, _, filter = items[:7]         
        if len(chrom) <= 6:
            site = chrom + ":" + pos + ":" + ref + ":" + alt
            truth_vars.add(site)
cr = list()
count = 0
with open(fastvc_file, 'r') as f, open(tsv_file, 'a') as tsv:
    for var in f:
        if var[0] == '#':
            continue
        items = var.strip().split('\t')
        if(len(items) == 61): #36 or 38(fisher)
            if vartype == "SNV" and items[fe2i["VarType"]] != "SNV":
                continue
            elif vartype == "INDEL" and items[fe2i["VarType"]] == "SNV":
                continue
            if vartype == "SNV":
                key, data = format_snv_data_item(items, True)
            elif vartype == "INDEL":
                key, data = format_indel_data_item(items, True)
            label = get_label(truth_vars, key)
            raw = ','.join([*data, str(label)])
            tsv.write(raw+'\n')
            count += 1
        else:
            print("wrong data format!!", len(items))

print("totally write %d data!" % count)

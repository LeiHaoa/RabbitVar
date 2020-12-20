#!/usr/bin/env python3
import sys
from multiprocessing import Process, Queue

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
import os

features_to_index = {
    "RefAllel":5, "VarAllel": 6, "TotalPosCov":7, "positioncoverage" : 8,"refForwardcoverage" : 9,"refReversecoverage" : 10,
    "varsCountOnForward" : 11, "VarsCountOnReverse" : 12, "genotype" : 13, "frequency" : 14, "strandbiasflag" : 15, 
    "meanPosition" : 16, "pstd" : 17, "meanQuality" : 18, "qstd" : 19, "pvalue" : 20, "ratio" : 21, "mapq" : 22, 
    "qratio" : 23, "higreq" : 24, "extrafreq" : 25, "shift3" : 26, "msi" : 27, "msint" : 28, "nm" : 29, "hicnt" : 30, 
    "hicov" : 31, "leftSequence" : 32, "rightSequence" : 33, "region" : 34, "varType" : 35, "duprate" : 36
}

All_features = ["RefAllel", "VarAllel", "TotalPosCov", "positioncoverage","refForwardcoverage","refReversecoverage",
    "varsCountOnForward", "VarsCountOnReverse","genotype","frequency","strandbiasflag","meanPosition","pstd","meanQuality",
    "qstd", "pvalue","ratio","mapq","qratio","higreq", "extrafreq","shift3","msi","msint","nm","hicnt","hicov",
    "leftSequence","rightSequence","region","varType","duprate"
]
selected_features = [
    "TotalPosCov", "positioncoverage","refForwardcoverage","refReversecoverage","varsCountOnForward", "VarsCountOnReverse",
    "frequency","meanPosition","pstd","meanQuality","pvalue","ratio","mapq","qratio","higreq", "shift3","msi","msint",
    "nm","hicnt","hicov",
]

type_to_label = {"SNV": 0, "Deletion": 1, "Insertion": 2, "Complex": 3}
fe2i = features_to_index
fvc_sf = selected_features

GERM_SNV_FEATURES=21
GERM_INDEL_FEATURES=24

def format_snv_data_item(jri, fisher):
    if not fisher:
        print("not support to train if you not run rabbitvar without --fiser!!")
        exit(-1)
    data = list()
    # key is chrom:pos like "chr1:131022:A:T"
    key = jri[2] + ":" + jri[3] + ":" + jri[5] + ":" + jri[6] #TODO: case sensitive
    for sf in fvc_sf:
        data.append(jri[fe2i[sf]])
    if len(data) != GERM_SNV_FEATURES:
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
    if len(data) != GERM_INDEL_FEATURES:
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
        if len(chrom) < 6:
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

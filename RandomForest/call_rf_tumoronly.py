import joblib
import numpy as np
import pandas as pd
import time 
import argparse

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

type_to_label = {"SNV": 0, "Deletion": 1, "Insertion": 2, "Complex": 3, "MNV":3}
varLabel_to_label = {
    "Germline":0, "StrongLOH":1, "LikelyLOH":2, "StrongSomatic":3,
    "LikelySomatic":4, "AFDiff":5, "SampleSpecific":6
}
GERM_SNV_FEATURES=21
GERM_INDEL_FEATURES=24

fe2i = features_to_index
fvc_sf = selected_features

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
    data.append(len(jri[5]))
    data.append(len(jri[6]))
    data.append(type_to_label[jri[fe2i["varType"]]])
    for sf in fvc_sf:
        data.append(jri[fe2i[sf]])
    #data.append(jri[fe2i["VarLabel"]])
    if len(data) != GERM_INDEL_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        exit(-1)
    return key, data

#out_file = "/home/old_home/haoz/workspace/FastVC/detection_result/somatic/train_set_10_21/FDSynthetic.notloose.txt"
def call_rf(args):
    #out_file = "/home/old_home/haoz/workspace/FastVC/detection_result/NC_DATASET/NC_DATA_1.txt"
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
            if(len(items) == 38): 
                if vartype == "SNV" and items[fe2i["varType"]] != "SNV":
                    continue
                elif vartype == "INDEL" and items[fe2i["varType"]] == "SNV":
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
    #clf = joblib.load("/home/old_home/haoz/workspace/FastVC/RandomForest/models/ramdom_forest_model_somatic_test0d02.pkl")
    clf = joblib.load(args.model)
    cr_end = time.time()
    print("time of load model: {} s".format(cr_end - cr_start) )

    cr_start = time.time()
    data = np.asfarray(cr)
    cr_end = time.time()
    print("time of np.asaray: {} s".format(cr_end - cr_start) )
    
    cr_start = time.time()
    pred = clf.predict(data)
    cr_end = time.time()
    print("time of pred: {} s".format(cr_end - cr_start) )

    cr_start = time.time()
    with open(args.out_file, 'w') as f:
        for i in range(len(pred)):
            if pred[i] == 1:
                f.write(str(raw[i]))
    cr_end = time.time()
    print("time of write: {} s".format(cr_end - cr_start) )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "train your network")
    parser.add_argument('--in_file', help = "RabbitVar intermidiate file(with fisher test)", type=str, required = True)
    parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--nthreads', help = "number of thread", type=int, default=20)
    parser.add_argument('--model', help = "random forest model", type=str, required = True)
    parser.add_argument('--out_file', help = "file to output", type=str, required = True)
    args = parser.parse_args()
    call_rf(args)

import sys
from multiprocessing import Process, Queue

import joblib
import numpy as np
import pandas as pd
import sklearn.datasets
import sklearn.metrics
import sklearn.model_selection
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
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
    data.append(len(jri[5]))
    data.append(len(jri[6]))
    data.append(type_to_label[jri[fe2i["varType"]]])
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

def train_rf(args):
    #----prepare label list----#
    #truth_file = "/home/old_home/haoz/workspace/data/FD/Truth/FDtruth_Data_1.snv.vcf"
    truth_file = args.truth_file
    vartype = args.var_type

    #--- prepare source datafram from file ---#

    #fastvc_file = "/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/FD_DATA_1.txt"
    fastvc_file = args.train_data
    cr = list()

    if args.tsv == None:
        truth_vars = set()
        with open(truth_file, 'r') as f:
            for var in f:
                if var[0] == '#':
                    continue
                items = var.split('\t')
                chrom, pos, id, ref, alt, _, filter = items[:7]         
                #if len(chrom) < 6 and filter == "PASS" and (len(ref) > 1 or len(alt) > 1) :
                if len(chrom) < 6:
                    site = chrom + ":" + pos + ":" + ref + ":" + alt
                    truth_vars.add(site)

        print("length of truth: ", len(truth_vars))
        with open(fastvc_file, 'r') as f:
            for var in f:
                if var[0] == '#':
                    continue
                items = var.strip().split('\t')
                if(len(items) == 38): #36 or 38(fisher)
                    if vartype == "SNV" and items[fe2i["varType"]] != "SNV":
                        continue
                    elif vartype == "INDEL" and items[fe2i["varType"]] == "SNV":
                        continue
                    if vartype == "SNV":
                        key, data = format_snv_data_item(items, True)
                    elif vartype == "INDEL":
                        key, data = format_indel_data_item(items, True)
                    label = get_label(truth_vars, key)
                    cr.append([*data, label])
                else:
                    print("wrong data format!!", len(items))
                    exit(0)
    else:
        for var in f:
            if var[0] == '#':
                continue
            cr.append(var.strip().split(','))
                

    print("length of cr: ", len(cr))
    if vartype == "SNV":
        data = pd.DataFrame(cr, columns=[*selected_features, "label"])
        data[data.columns] = data[data.columns].apply(pd.to_numeric)

        #--- data prepare ----#
        #data.rename(columns={0:'input',1:'label'},inplace=True)
        data['is_train'] = np.random.uniform(0, 1, len(data)) <= .9
        train, test = data[data['is_train'] == True], data[data['is_train'] == False]
        
        #--select 1/100 faslse data form fasle
        false = train[ train['label'] == 0 ]
        false = false.sample(frac = 0.02)
        truth = train[ train['label'] == 1]
        print("faslse number: {}, truth number: {}".format(len(false), len(truth)) )
        aug_data = shuffle(pd.concat([truth, false, truth, truth, truth, truth], axis = 0))
        #aug_data = shuffle(pd.concat([truth, false], axis = 0))

        #aug_data = train
        train_set = aug_data[aug_data.columns[:-2]].to_numpy()
        #y, _ = pd.factorize(aug_data["label"])
        y = aug_data["label"]#.to_numpy()
        print("type: ", type(train_set), type(y))
        print(len(train_set), len(y))
    elif vartype == "INDEL":
        data = pd.DataFrame(cr, columns=["RefLength", "AltLength", "VarType", *selected_features, "label"])
        data[data.columns] = data[data.columns].apply(pd.to_numeric)

        #--- data prepare ----#
        #data.rename(columns={0:'input',1:'label'},inplace=True)
        data['is_train'] = np.random.uniform(0, 1, len(data)) <= .9
        train, test = data[data['is_train'] == True], data[data['is_train'] == False]
        
        #--select 1/100 faslse data form fasle
        false = train[ train['label'] == 0 ]
        false = false.sample(frac = 0.02)
        truth = train[ train['label'] == 1]
        print("faslse number: {}, truth number: {}".format(len(false), len(truth)) )
        aug_data = shuffle(pd.concat([truth, false, truth, truth, truth, truth], axis = 0))
        #aug_data = shuffle(pd.concat([truth, false], axis = 0))
        #aug_data = train
        train_set = aug_data[aug_data.columns[:-2]].to_numpy()
        print(train_set[:2])
        y = aug_data["label"]
        print("type: ", type(train_set), type(y))
        print(len(train_set), len(y))

    print("data prepare done, start fiting ...")
    if args.pretrained_model == None:
        clf = RandomForestClassifier(n_jobs=args.nthreads, max_depth=12, min_samples_leaf=50, 
                                     class_weight = 'balanced',
                                     n_estimators=50, max_features=None, verbose=1)
    else:
        print("loading model: ", args.pretrained_model)
        clf = joblib.load(args.pretrained_model)

    clf.fit(train_set, y)
    
    joblib.dump(clf, args.out_model, compress=9)
    print("store model done, now testing...")
    #--- test
    #test_set = np.asarray(list(map(lambda x: np.asarray(x), test["input"])))
    test_set = test[test.columns[:-2]].to_numpy()
    #test_set_reduction = pca.transform(test_set)
    ytest, _ = pd.factorize(test["label"])
    ypred = clf.predict(test_set)
    print("Accuracy score", sklearn.metrics.accuracy_score(ytest, ypred))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "train your network")
    parser.add_argument('--train_data', help = "RabbitVar intermidiate file(with fisher test)", type=str, required = True)
    parser.add_argument('--truth_file', help = "truth file / the ground truth(.vcf)", type=str, required = True)
    parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--tsv', help = "tsv file, which contained data and label", type=str, required = False)
    #parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--nthreads', help = "number of thread", type=int, default=20)
    parser.add_argument('--pretrained_model', help = "pretrained model", type=str, required = False)
    parser.add_argument('--out_model', help = "out model name (just for experiments)", type=str, required = True)
    args = parser.parse_args()
    train_rf(args)

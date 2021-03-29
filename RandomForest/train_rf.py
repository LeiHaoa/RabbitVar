import sys
import os
from multiprocessing import Process, Queue
sys.path.append(os.path.join(os.getcwd(), '../utils'))

import joblib
import numpy as np
import pandas as pd
import sklearn.datasets
import sklearn.metrics
import sklearn.model_selection
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
import argparse

from features import *
from datautil import get_data_fromcsv, hard_filter

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
    data.append(len(jri[5]))
    data.append(len(jri[6]))
    data.append(type_to_label[jri[fe2i["VarType"]]])
    for sf in fvc_sf:
        data.append(jri[fe2i[sf]])
    data.append(varLabel_to_label[jri[fe2i["VarLabel"]]])#varlabel
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

def train_rf_old(args):
    #----prepare label list----#
    #truth_file = "/home/old_home/haoz/workspace/data/FD/Truth/FDtruth_Data_1.snv.vcf"
    truth_file = args.truth_file
    vartype = args.var_type

    #--- prepare source datafram from file ---#

    #fastvc_file = "/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/FD_DATA_1.txt"
    fastvc_file = args.train_data
    #fastvc_file = "/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/demo.txt"
    cr = list()
    if args.tsv == None:  #use file and truth file 
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
                    cr.append([*data, label])
                else:
                    print("wrong data format!!", len(items))
    else: #use tsv file 
        '''
        with open(args.tsv, 'r') as f:
            for var in f:
                if var[0] == '#':
                    continue
                cr.append(var.strip().split(','))
        '''
        #cr = np.fromfile(args.tsv) much slower than pd.read_csv
        data = pd.read_csv(args.tsv, header=None)
        print("length of data: ", len(data))

    if vartype == "SNV":
        if args.tsv == None:
            data = pd.DataFrame(cr, columns=[*som_selected_features, "VarLabel", "label"])
            data[data.columns] = data[data.columns].apply(pd.to_numeric)
            #data['VarLabel'] = data['VarLabel'].astype('category')
            #data["label"] = data["label"]
        else:
            data.columns = [*som_selected_features, "VarLabel", "label"]

        #--- data prepare ----#
        #data.rename(columns={0:'input',1:'label'},inplace=True)
        data['is_train'] = np.random.uniform(0, 1, len(data)) <= .9
        train, test = data[data['is_train'] == True], data[data['is_train'] == False]
        
        #--select 1/100 faslse data form fasle
        false = train[ train['label'] == 0 ]
        #false = false.sample(frac = 0.2)
        truth = train[ train['label'] == 1]
        print("faslse number: {}, truth number: {}".format(len(false), len(truth)) )
        aug_data = shuffle(pd.concat([truth, false, truth, truth, truth, truth], axis = 0))
        # aug_data = shuffle(pd.concat([truth, false], axis = 0))
        # aug_data = train

        train_set = aug_data[aug_data.columns[:-2]].to_numpy()
        #y, _ = pd.factorize(aug_data["label"])
        y = aug_data["label"]#.to_numpy()
        print("type: ", type(train_set), type(y))
        print(len(train_set), len(y))
    elif vartype == "INDEL":
        if args.tsv == None:
            data = pd.DataFrame(cr, columns=["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel", "label"])
            data[data.columns[:-2]] = data[data.columns[:-2]].apply(pd.to_numeric)
            #data["label"] = data["label"].astype('category')
            data[data.columns] = data[data.columns].apply(pd.to_numeric)
        else:
            data.columns = ["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel", "label"]

        #--- data prepare ----#
        #data.rename(columns={0:'input',1:'label'},inplace=True)
        data['is_train'] = np.random.uniform(0, 1, len(data)) <= .9
        train, test = data[data['is_train'] == True], data[data['is_train'] == False]

        #--select frac faslse data form fasle
        false = train[ train['label'] == 0 ]
        truth = train[ train['label'] == 1 ]
        print("faslse number: {}, truth number: {}".format(len(false), len(truth)*5 ) )
        aug_data = shuffle(pd.concat([truth, false, truth, truth, truth, truth], axis = 0))
        #aug_data = train
        train_set = aug_data[aug_data.columns[:-2]].to_numpy()
        y = aug_data["label"]
        print("type: ", type(train_set), type(y))
        print(len(train_set), len(y))

    print("data prepare done, start fiting ...")
    if args.pretrained_model == None:
        clf = RandomForestClassifier(n_jobs=args.nthreads, max_depth=20, min_samples_leaf=50, 
                                     n_estimators=150, max_features=None, verbose=1)
    else:
        print("model {} exists!".format(args.pretrained_model))
        exit(-1)
        if not os.path.exists(args.pretrained_model):
            print("model {} not exists, create new one!".format(args.pretrained_model))
            clf = RandomForestClassifier(n_jobs=args.nthreads, max_depth=20, min_samples_leaf=50, 
                                         n_estimators=150, max_features=None, verbose=1)
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

def train_rf(args):
    #----prepare label list----#
    truth_file = args.truth_file
    vartype = args.var_type
    #--- prepare source datafram from file ---#
    fastvc_file = args.train_data
    cr = list()

    if vartype == "SNV":
        if args.tsv == None:
            data = pd.DataFrame(cr, columns=[*som_selected_features, "VarLabel", "label"])
            data[data.columns] = data[data.columns].apply(pd.to_numeric)
            #data['VarLabel'] = data['VarLabel'].astype('category')
            #data["label"] = data["label"]
        else:
            data = get_data_fromcsv(args.tsv, [*som_selected_features, "VarLabel", "label"], vtype = 'SNV')

        data = hard_filter(data)

        #--- data prepare ----#
        false = data[ data['label'] == 0 ]
        truth = data[ data['label'] == 1 ]
        print("faslse number: {}, truth number: {}".format(len(false), len(truth)) )
        aug_data = shuffle(pd.concat([truth, false, truth, truth], axis = 0))
        y = aug_data["label"]#.to_numpy()
        print(len(aug_data), len(y))
    elif vartype == "INDEL":
        if args.tsv == None:
            data = pd.DataFrame(cr, columns=["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel", "label"])
            data[data.columns[:-2]] = data[data.columns[:-2]].apply(pd.to_numeric)
            #data["label"] = data["label"].astype('category')
            data[data.columns] = data[data.columns].apply(pd.to_numeric)
        else:
            data = get_data_fromcsv(args.tsv, ["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel", "label"], vtype = 'SNV')

        data = hard_filter(data)
        #--- data prepare ----#
        false = data[ data['label'] == 0 ]
        truth = data[ data['label'] == 1 ]
        print("faslse number: {}, truth number: {}".format(len(false), len(truth)) )
        aug_data = shuffle(pd.concat([truth, false, truth, truth], axis = 0))
        y = aug_data["label"]#.to_numpy()
        print(len(aug_data), len(y))

    print("data prepare done, start fiting ...")
    if args.pretrained_model == None:
        clf = RandomForestClassifier(n_jobs=args.nthreads, max_depth=20, min_samples_leaf=50, 
                                     n_estimators=150, max_features=None)
    else:
        print("model {} exists!".format(args.pretrained_model))
        exit(-1)
        if not os.path.exists(args.pretrained_model):
            print("model {} not exists, create new one!".format(args.pretrained_model))
            clf = RandomForestClassifier(n_jobs=args.nthreads, max_depth=20, min_samples_leaf=50, 
                                         n_estimators=150, max_features=None)
        else:
            print("loading model: ", args.pretrained_model)
            clf = joblib.load(args.pretrained_model)

    clf.fit(aug_data, y)
    
    joblib.dump(clf, args.out_model, compress=9)
    print("store model done")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "train your network")
    parser.add_argument('--train_data', help = "RabbitVar intermidiate file(with fisher test)", type=str, required = False)
    parser.add_argument('--truth_file', help = "truth file / the ground truth(.vcf)", type=str, required = False)
    parser.add_argument('--tsv', help = "tsv file, which contained data and label", type=str, required = False)
    parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    #parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--nthreads', help = "number of thread, default -1", type=int, default=-1)
    parser.add_argument('--pretrained_model', help = "pretrained model", type=str, required = False)
    parser.add_argument('--out_model', help = "out model name (just for experiments)", type=str, required = True)
    args = parser.parse_args()
    train_rf(args)

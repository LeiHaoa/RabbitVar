import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../utils'))
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
import features
import datautil

#---------------xgboost
import xgboost as xgb
from sklearn.metrics import mean_squared_error

def train_rf(args):
    #----prepare label list----#
    truth_file = args.truth_file
    vartype = args.var_type
    #--- prepare source datafram from file ---#
    fastvc_file = args.train_data

    if vartype == "SNV":
        data = get_data_fromcsv(args.tsv, columns=[*features.som_rf_snv_input_features, 'Label'], vtype = 'SNV')
        print("before hard filter: {} data".format(len(data)))
        #data = hard_filter(data)
        print("after hard filter: {} data".format(len(data)))
        test_data = get_data_fromcsv(args.tsv + ".test", columns=[*features.som_rf_snv_input_features, "Label"], vtype = "SNV")
        #test_data = hard_filter(test_data)
        tmp_feature = features.som_rf_snv_input_features
    elif vartype == "INDEL":
        data = get_data_fromcsv(args.tsv, columns=[*features.som_rf_indel_input_features, 'Label'], vtype = 'INDEL')
        print("before hard filter: {} data".format(len(data)))
        data = hard_filter(data)
        print("after hard filter: {} data".format(len(data)))
        #-----------test data -------------#
        #test_data = get_data_fromcsv(args.tsv + ".test", columns=[*features.som_rf_indel_input_features, "Label"], vtype = "INDEL")
        test_data = get_data_fromcsv("./data_INDEL_HCC1395_ALL.tsv.test", columns=[*features.som_rf_indel_input_features, "Label"], vtype = "INDEL")
        test_data = hard_filter(test_data)
        tmp_feature = features.som_rf_indel_input_features
    else:
        print(f'undefied variant type: {vartype}')
        exit(-1)

    print(data.head())
    print("data prepare done, start fiting ...")

    #clf = RandomForestClassifier(n_jobs=args.nthreads, max_depth=len(features.som_rf_indel_input_features), min_samples_leaf=50, 
    #                             n_estimators=150, max_features=None, class_weight = "balanced")
    #-----xgboost------#
    xg_w_scale = sum(data['Label'] == 0) / sum(data['Label'] == 1)
    print("using xgboost pos weight: ", xg_w_scale)
    # snv xgboost param
    # clf = xgb.XGBClassifier(objective ='binary:logistic', scale_pos_weight = xg_w_scale, 
    #                         learning_rate = 0.1, max_depth = 10, min_child_weight = 17,
    #                         n_jobs = -1, n_estimators = 1300,
    #                         tree_method = 'gpu_hist', gpu_id = 0
    #                         )
    # indel xgboost param
    clf = xgb.XGBClassifier(objective ='binary:logistic', scale_pos_weight = xg_w_scale, 
                            learning_rate = 0.01, max_depth = 25, min_child_weight = 18,
                            n_jobs = -1, n_estimators = 1200,
                            tree_method = 'gpu_hist', gpu_id = 1
                            )

    clf.fit(data[tmp_feature], data['Label'])
    
    joblib.dump(clf, args.out_model, compress=9)
    print("store model done, now testing...")
    test_set = test_data[tmp_feature]
    ytest = test_data['Label']
    ypred = clf.predict(test_set)
    print("Accuracy score", sklearn.metrics.accuracy_score(ytest, ypred))

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

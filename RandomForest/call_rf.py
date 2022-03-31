import joblib
import numpy as np
import pandas as pd
import time 
import argparse
import subprocess

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../utils'))
from features import *
from datautil import hard_filter
import datautil
from vcf_writer import *

all_inform = False

if all_inform:
    SOM_SNV_FEATURES=len(som_selected_features) + 1 + (5 + 2)
else:
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
        data.append(float(jri[fe2i[sf]]))
    data.append(varLabel_to_label[jri[fe2i["VarLabel"]]])#varlabel
    if len(data) != SOM_SNV_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        exit(-1)
    return key, data

def get_gt(af):
    GTFREQ = 0.2
    FREQ  = 0.01
    af = float(af)
    if 1 - af < GTFREQ:
        gt = "3"
    else:
        if af <= 0.5:
            gt = "2"
        else:
            if af > FREQ:
                gt = "1"
            else:
                gt = "0"
    return gt

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

    if all_inform:
        t_vb = jri[fe2i['Var1Bias']].split(";") #var1bias
        var1sbr, var1sba = t_vb if len(t_vb) == 2 else [t_vb[0], t_vb[0]]
        data.append(var1sbr)
        data.append(var1sba)
        t_vb = jri[fe2i['Var2Bias']].split(";") #var2bias
        var2sbr, var2sba = t_vb if len(t_vb) == 2 else [t_vb[0], t_vb[0]] 
        data.append(var2sbr)
        data.append(var2sba)
        data.append( "1" if jri[fe2i['Var1Genotype']] == jri[fe2i['Var2Genotype']] else "0") # SameGenotype
        data.append(get_gt(jri[fe2i["Var1AF"]])) #gt
        data.append(get_gt(jri[fe2i["Var2AF"]])) #gtm

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
    return np.asarray([1 if x >= scale else 0 for x in proba[:,1]])

def depth_normalization(data, tdepth, ndepth):
    #to_be_scale = ["Var1AltDepth", "Var1RefFwdReads", "Var1RefRevReads", "Var1AltFwdReads", "Var1AltRevReads", "Var2AltDepth", "Var2RefFwdReads", "Var2RefRevReads", "Var2AltFwdReads", "Var2AltRevReads"]
    data['Var1AltFreq']    = data['Var1AltDepth'] / data['Var1Depth']
    data['Var1RefFwdFreq'] = data['Var1RefFwdReads'] / (data['Var1Depth'] - data['Var1AltDepth'])
    data['Var1AltFwdFreq'] = data['Var1AltFwdReads'] / data['Var1AltDepth']
    data['Var2AltFreq']    = data['Var2AltDepth'] / data['Var2Depth']
    data['Var2RefFwdFreq'] = data['Var2RefFwdReads'] / (data['Var2Depth'] - data['Var2AltDepth'])
    data['Var2AltFwdFreq'] = data['Var2AltFwdReads'] / data['Var2AltDepth']
    data = data.replace([np.inf, -np.inf], np.nan).fillna(-1)
    data['Var1DepthFreq'] = data['Var1Depth'] / tdepth
    data['Var2DepthFreq'] = data['Var2Depth'] / ndepth
    return data 

def rf_filter(args):
  cr = list()
  raw = list()
  in_file = args.in_file
  scale = float(args.scale)
  time_step0 = time.time()
  #cr = pd.read_csv(in_file, delimiter = '\t', header = None, engine = 'c', skipinitialspace = True, na_filter = False)
  #cr.columns = [*som_features, 'None'] #TODO: i should change the code of c++ to avoid the None colum
  cr = datautil.get_data_fromtxt(in_file, args.var_type)
  cr['VarLabel'] = cr['VarLabel'].map(varLabel_to_label)
  cr['VarType'] = cr['VarType'].map(type_to_label)
  cr['RefLength'] = cr['Ref'].str.len()
  cr['AltLength'] = cr['Alt'].str.len()

  time_step1 = time.time()
  print("load time: ", time_step1 - time_step0)
  #hard filter
  print('before hardf:', len(cr))
  cr = hard_filter(cr)
  print('after hardf:', len(cr))
  time_step2 = time.time()
  print("hard filter time: ", time_step2 - time_step1)

  #snv data process 
  if args.var_type == 'SNV':
    time_start = time.time()
    snvs = cr[cr['VarType'] == 0]
    if scale == 0:
      print('Just use hard filter')
      return snvs
    inputs = snvs[[*som_selected_features, "VarLabel"]].to_numpy()
    clf = joblib.load(args.models)
    clf.verbose = False
    snv_pred = my_predict(clf, inputs, scale)
    snv_result = snvs.loc[snv_pred == 1]
    time_end = time.time()
    print("time filter snv: {} s".format(time_end - time_start))
    return snv_result
  #indel data process
  elif args.var_type == 'INDEL':
    time_start = time.time()
    indels = cr[cr['VarType'] != 0]
    if scale == 0:
      print('Just use hard filter')
      return indels

    #iii = ["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel"]
    #inputs = indels[["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel"]]
    #depth normal lization
    # indels = depth_normalization(indels, args.tdepth, args.ndepth)
    #tmp = indels[["Var1DepthFreq", "Var1AltFreq", "Var1RefFwdFreq", "Var1AltFwdFreq",  "Var1AF","Var2DepthFreq", "Var2AltFreq", "Var2RefFwdFreq", "Var2AltFwdFreq",  "Var2AF"]]
    #tmp.to_csv("./result_vcf.csv", sep="\t", index='False', encoding='ascii', float_format='%10.4f')
    inputs = indels[som_rf_input_features]
    print("inputs size: ", len(inputs))
    #clf = joblib.load(args.indel_model)
    clf = joblib.load(args.model)
    clf.verbose = False
    indel_pred = my_predict(clf, inputs, scale)
    indel_result = indels.loc[indel_pred == 1]
    time_end = time.time()
    print("time filter indel: {} s, indel size: {}".format(time_end - time_start, len(indel_result)))
    return indel_result
  else: 
    print("ERROR: Unsupported variant type:", args.var_type)
    exit(-1)

def call_rf(args):
  vcf = rf_filter(args)
  cr_start = time.time()
  vcf['Chr'] = vcf['Chr'].astype(str)
  vcf = vcf.sort_values(by=['Chr', 'Start'])
  print(vcf.head())
  vcf_file = args.out_file
  with open(vcf_file, 'w') as f:
    write_header(f)
    tmp = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", 'Sample']) #TODO what if user not specified a sample name???
    tmp += '\n'
    f.write(tmp)
    for i, record in vcf.iterrows():
      f.write(format_record(record) + '\n')
  cr_end = time.time()
  print("time of write: {} s".format(cr_end - cr_start) )


def call_rf_old(args):
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
    #clf = joblib.load("/home/old_home/haoz/workspace/FastVC/RandomForest/models/ramdom_forest_model_somatic_test0d02.pkl")
    clf = joblib.load(args.model)
    cr_end = time.time()
    print("time of load model: {} s".format(cr_end - cr_start) )

    cr_start = time.time()
    #data = np.asfarray(cr)
    data = pd.DataFrame(cr)
    if vartype == 'INDEL':
        if all_inform:
            data.columns = ["RefLength", "AltLength", "VarType", "Var1SBR", "Var1SBA", "Var2SBR", "Var2SBA", "SameGenotype", "GT", "GTM", *som_selected_features, "VarLabel"]
        else:
            data.columns = ["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel"]
    else:
        data.columns = [*som_selected_features, "VarLabel"]
    cr_end = time.time()
    print("time of np.asaray: {} s".format(cr_end - cr_start) )
    
    #--hardfilter--#
    print("before hard filter: {} data".format(len(data)))
    data = hard_filter_keeporg(data)
    print("after hard filter: {} data".format(len(data)))
    ###-----------test for not depth data-------------------#
    #data = data[["RefLength", "AltLength", "VarType", *som_rf_input_features, "VarLabel", "hard_flag"]]
    ###
    cr_start = time.time()
    scale = args.scale
    pred = my_predict(clf, data[data.columns[:-1]], scale) #remove hard flag
    data['pred'] = pred
    cr_end = time.time()
    print("length {} - {}".format(len(data), len(pred)))
    print("time of pred: {} s".format(cr_end - cr_start) )

    cr_start = time.time()
    with open(args.out_file, 'w') as f:
        for i in range(len(pred)):
            if pred[i] == 1 and data['hard_flag'][i] == 0:
                #if data['hard_flag'][i] == 0:
                f.write(str(raw[i]))
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
    parser.add_argument('--tdepth', help = "tdepth", type=int, required = True)
    parser.add_argument('--ndepth', help = "ndepth", type=int, required = True)
    args = parser.parse_args()
    call_rf(args)

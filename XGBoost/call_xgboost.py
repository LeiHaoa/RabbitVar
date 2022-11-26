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

#---------------xgboost
import xgboost as xgb
from sklearn.metrics import mean_squared_error

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

def my_predict_uniform(clf, data, scale):
    tmp = list(data.columns)
    tmp.remove('Var1AF')
    data = data[tmp]
    proba = clf.predict_proba(data)[:,1]
    scale = float(scale)
    return proba, np.asarray([1 if x >= scale else 0 for x in proba])

def my_predict(clf, data, scale, depth):
    start_t = time.time()
    af = data['Var1AF'].to_numpy()
    #tmp = list(data.columns)
    #tmp.remove('Var1AF')
    #data = data[tmp]

    proba = clf.predict_proba(data)[:,1]
    res = []
    assert(len(af) == len(proba))
    #sp = float(scale.split(':')[0])
    sp, lowaf_scale, highaf_scale = [float(x) for x in scale.split(':')]
    #lowaf_scale = 0.9 + min(0.09, 0.09 * depth / 300)
    print(f"sp: {sp}, lowaf_scale: {lowaf_scale}, highaf_scale: {highaf_scale}")
    for i in range(len(proba)):
      #if (af[i] < sp and proba[i] > lowaf_scale) or (af[i] >= sp and proba[i] > highaf_scale):
      if (af[i] < sp and proba[i] > (1 - math.pow(0.1, 3 - af[i]*10)) ) or (af[i] >= sp and proba[i] > highaf_scale):
        res.append(1)
      else:
        res.append(0)
    end_t = time.time()
    print(f'predict time: {end_t - start_t}s')

    return proba, np.asarray(res)


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
  scale = args.scale
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

  tdepth = args.tdepth
  ndepth = args.ndepth
  if os.path.exists(args.in_file+'.info'):
    with open(args.in_file + ".info", "r") as f:
      tdepth = round(float(f.readline().strip()))
      ndepth = round(float(f.readline().strip()))  
      print("find .info file!, tdepth: {}, ndepth: {}".format(tdepth, ndepth))

  #snv data process 
  if args.var_type == 'SNV':
    time_start = time.time()
    snvs = cr[cr['VarType'] == 0]
    if scale == '0':
      print('Just use hard filter')
      return snvs
    #inputs = snvs[[*som_rf_snv_input_features, "Var1AF"]]
    inputs = snvs[som_rf_snv_input_features]
    clf = joblib.load(args.model)
    clf.verbose = False
    if len(scale.split(':')) == 3:
      snv_proba, snv_pred = my_predict(clf, inputs, scale, tdepth)
    elif len(scale.split(':')) == 1:
      snv_proba, snv_pred = my_predict_uniform(clf, inputs, scale)
    else:
      print(f"error: unsupported scale format: {scale}")
    snvs['pred'] = snv_proba
    snv_result = snvs.loc[snv_pred == 1]
    time_end = time.time()
    print("time filter snv: {} s".format(time_end - time_start))
    return snv_result
  #indel data process
  elif args.var_type == 'INDEL':
    time_start = time.time()
    indels = cr[cr['VarType'] != 0]
    if scale == '0':
      print('Just use hard filter')
      return indels

    #depth normal lization
    #indels = depth_normalization(indels, tdepth, ndepth)
    tmp_feature = som_rf_indel_input_features
    inputs = indels[tmp_feature]
    print("inputs size: ", len(inputs))
    clf = joblib.load(args.model)
    load_time = time.time()
    print(f'load model time: {load_time - time_start}')
    clf.verbose = False
    if len(scale.split(':')) == 3:
        indel_proba, indel_pred = my_predict(clf, inputs, scale, tdepth)
    elif len(scale.split(':')) == 1:
        indel_proba, indel_pred = my_predict_uniform(clf, inputs, scale)
    else:
        print(f"error: unsupported scale format: {scale}")
    indels['pred'] = indel_proba
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
  vcf_file = args.out_file
  just_hard_filter = True if args.scale == '0' else False
  with open(vcf_file, 'w') as f:
    write_header(f)
    tmp = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", 'Sample']) #TODO what if user not specified a sample name???
    tmp += '\n'
    f.write(tmp)
    for i, record in vcf.iterrows():
      f.write(format_record(record, just_hard_filter) + '\n')
  cr_end = time.time()
  print("time of write: {} s".format(cr_end - cr_start) )


def call_rf_keep_all(args):
  cr = list()
  raw = list()
  in_file = args.in_file
  scale = args.scale
  time_step0 = time.time()
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

  tdepth = args.tdepth
  ndepth = args.ndepth
  if os.path.exists(args.in_file+'.info'):
    with open(args.in_file + ".info", "r") as f:
      tdepth = round(float(f.readline().strip()))
      ndepth = round(float(f.readline().strip()))  
      print("find .info file!, tdepth: {}, ndepth: {}".format(tdepth, ndepth))

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
    #snv_result = snvs.loc[snv_pred == 1]
    snvs['pred'] = snv_pred
    time_end = time.time()
    print("time filter snv: {} s".format(time_end - time_start))
    vcf = snvs
  #indel data process
  elif args.var_type == 'INDEL':
    time_start = time.time()
    indels = cr[cr['VarType'] != 0]
    if scale == 0:
      print('Just use hard filter')
      return indels

    #depth normal lization
    #indels = depth_normalization(indels, tdepth, ndepth)
    inputs = indels[som_rf_indel_input_features]
    print("inputs size: ", len(inputs))
    clf = joblib.load(args.model)
    clf.verbose = False
    #indel_pred = my_predict(clf, inputs, scale)
    indel_pred = clf.predict_proba(inputs) 
    #indel_result = indels.loc[indel_pred == 1]
    indels['pred'] = indel_pred[:,1]
    time_end = time.time()
    print("time filter indel: {} s, indel size: {}".format(time_end - time_start, len(indels)))
    vcf = indels
  else: 
    print("ERROR: Unsupported variant type:", args.var_type)
    exit(-1)

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
      f.write(format_record_keep(record, just_hard_filter) + '\n')
  cr_end = time.time()
  print("time of write: {} s".format(cr_end - cr_start) )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "train your network")
    parser.add_argument('--in_file', help = "RabbitVar intermidiate file(with fisher test)", type=str, required = True)
    parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--scale', help = 'scale for xgboost predict_proba to filter, 0.5 perform the same result as function predict, formula: a:b-c means when VAF < a, chose variants whose predict prob > b and when VAF >= a, chose variants whose predict prob > c. (default 0.2:0.9-0.25)', type = str, default = "0.2:0.9-0.25")
    #parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--nthreads', help = "number of thread", type=int, default=20)
    parser.add_argument('--model', help = "random forest model", type=str, required = True)
    parser.add_argument('--out_file', help = "file to output", type=str, required = True)
    parser.add_argument('--tdepth', help = "tdepth", type=int, required = True)
    parser.add_argument('--ndepth', help = "ndepth", type=int, required = True)
    args = parser.parse_args()
    call_rf(args)
    #call_rf_keep_all(args)

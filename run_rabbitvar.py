#!/usr/bin/env python3
import time 
import os
import joblib
import numpy as np
import pandas as pd
import argparse
import subprocess
#import pybedtools
import math
import re
from cmdparser import *
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), './utils'))
from vcf_writer import *
import features
from datautil import *

fe2i = features.som_features_to_index
fvc_sf = features.som_selected_features

SOM_SNV_FEATURES=len(fvc_sf) + 1
SOM_INDEL_FEATURES = SOM_SNV_FEATURES + 3

type_to_label = {"SNV": 0, "Deletion": 1, "Insertion": 2, "Complex": 3, "MNV":3}
varLabel_to_label = {
    "Germline":0, "StrongLOH":1, "LikelyLOH":2, "StrongSomatic":3,
    "LikelySomatic":4, "AFDiff":5, "SampleSpecific":6
}

def get_socks():
  for line in str(subprocess.check_output('lscpu')).split('\\n'):
    contex = line.split()
    if contex[0].strip() == 'NUMA':
      return int(contex[-1])
  print("can not find NUMA struct, use default: 1")
  return 1

def prepare_cmd(BIN, bed_file_name, out_file_name, param):
  cmd = list()
  if 'bed' in param:
    param['bed'] = bed_file_name
  param['out'] = out_file_name
  cmd.append(BIN)
  for p, v in param.items():
    cmd.append("--" + p)
    cmd.append(v)
  
  return cmd

def bed_save_to(regions, path):
  with open(path, 'w') as f:
    for line in regions:
      f.write(str(line))

def split_bed(beds, parts, workspace):
  splited_info = []
  stride = len(beds) // parts
  for i in range(parts - 1):
    bed_path = os.path.join(workspace, "part{}.bed".format(i))
    out_path = os.path.join(workspace, "out{}.txt".format(i))
    bed_save_to(beds[(i)*stride:(i+1)*stride], bed_path)
    splited_info.append((bed_path, out_path))
  bed_save_to(beds[(parts - 1) * stride:], os.path.join(workspace, "part{}.bed".format(parts - 1)))
  splited_info.append((os.path.join(workspace, "part{}.bed".format(parts - 1)),
                       os.path.join(workspace, "out{}.txt".format(parts - 1)))
                     )
  print("split info: ", splited_info)
  return splited_info

def run_rabbitvar(BIN, workspace, param):
  socks = get_socks()
  splited_info = list()
  if not os.path.exists(workspace):
    os.mkdir(workspace)
  if socks > 1:
    #re-distribure bed file
    #beds = pybedtools.example_bedtool(param['bed'])
    #sorted(beds)
    beds = list()
    with open(param['bed'], 'r') as f:
      for line in f:
        beds.append(line)
    #split_info:[(bed1, out1), (bed2, out2), ...]
    splited_info = split_bed(beds, socks, workspace) #split file according to numa
    return splited_info
    ps = list()
    for i in range(int(socks)):
      cmd = prepare_cmd(BIN, splited_info[i][0], splited_info[i][1], param)
      print(cmd)
      try:
        ps.append(subprocess.Popen(cmd, stderr=subprocess.STDOUT))
      except subprocess.CalledProcessError as e:
        print("an error occured while runing cmd {}, error: {}".format(' '.join(cmd), e.output))
    for p in ps:
      p.wait()
  else:
    bed = param['bed'] if 'bed' in param else None
    out = param['out']
    cmd = prepare_cmd(BIN, bed, out, param)
    #print(cmd)
    splited_info.append((bed, out))
    subprocess.Popen(cmd, stderr=subprocess.STDOUT)

  print("All process run over!")
  return splited_info

def rf_filter(param, in_file):
  cr = list()
  raw = list()
  cr = pd.read_csv(in_file, delimiter = '\t', header = None, engine = 'c', skipinitialspace = True)
  cr.columns = [*features.som_features, 'None'] #TODO: i should change the code of c++ to avoid the None colum
  cr['VarLabel'] = cr['VarLabel'].map(varLabel_to_label)
  cr['VarType'] = cr['VarType'].map(type_to_label)
  cr['RefLength'] = cr['Ref'].str.len()
  cr['AltLength'] = cr['Alt'].str.len()
  #hard filter
  cr = hard_filter(cr)

  #snv data process 
  time_start = time.time()
  snvs = cr[cr['VarType'] == 0]
  inputs = snvs[[*som_selected_features, "VarLabel"]].to_numpy()
  clf = joblib.load(param['snvmod'])
  clf.verbose = False
  scale = float(param['snvscale'])
  snv_pred = my_predict(clf, inputs, scale)
  snv_result = snvs.loc[snv_pred == 1]
  time_end = time.time()
  print("time filter snv: {} s".format(time_end - time_start))
  
  #indel data process
  time_start = time.time()
  indels = cr[cr['VarType'] != 0]
  iii = ["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel"]
  inputs = indels[["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel"]].to_numpy()
  #clf = joblib.load(args.indel_model)
  scale = float(param['indelscale'])
  clf = joblib.load(param['indelmod'])
  clf.verbose = False
  indel_pred = my_predict(clf, inputs, scale)
  indel_result = indels.loc[indel_pred == 1]
  time_end = time.time()
  print("time filter indel: {} s".format(time_end - time_start))

  return snv_result, indel_result

def my_predict(clf, data, scale):
  proba = clf.predict_proba(data)
  return np.asarray([1 if x > scale else 0 for x in proba[:,1]])

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "rabbitvar")
  detector_parser = parser.add_argument_group("detector_parser")
  detectorParam(detector_parser)
  rabbitvar_parser = parser.add_argument_group("rabbitvar_parser")
  rabbitvarParam(rabbitvar_parser)
  filter_parser = parser.add_argument_group("filter_parser")
  filterParam(filter_parser)
  
  args = parser.parse_args()
  detector_param = {}
  for x in detector_parser._group_actions:
    k, v = x.dest, getattr(args, x.dest, None)
    if v:
      detector_param[k] = '' if isinstance(v, bool) else str(v)
  print("[INFO] detector paramters: ", detector_param)
  splited_info = run_rabbitvar(args.BIN, args.workspace, detector_param)
  filter_param = {}
  for x in filter_parser._group_actions:
    k, v = x.dest, getattr(args, x.dest, None)
    #print('filter parser:', k, v)
    if v:
      filter_param[k] = '' if isinstance(v, bool) else str(v)

  vcf_file = args.vcf
  vcf_list = []
  for detector_out in splited_info:
    vcf_list.extend(rf_filter(filter_param, detector_out[1]))
  vcf = pd.concat(vcf_list)
 
  #sort and write vcf file
  vcf = vcf.sort_values(['Chr', 'Start'])
  with open(vcf_file, 'w') as f:
    write_header(f)
    tmp = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", args.Name]) #TODO what if user not specified a sample name???
    tmp += '\n'
    f.write(tmp)
    for i, record in vcf.iterrows():
      f.write(format_record(record) + '\n')

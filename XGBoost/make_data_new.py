#!/usr/bin/env python3
import joblib
import numpy as np
import pandas as pd
import time 
import argparse
import subprocess

import sys
import os
sys.path.append('/home/user_home/haoz/workspace/RabbitVar/utils/')
from features import *
import features
from datautil import hard_filter
import datautil
from vcf_writer import *

type_to_label = {"SNV": 0, "Deletion": 1, "Insertion": 2, "Complex": 3, "MNV":3}
varLabel_to_label = {
    "Germline":0, "StrongLOH":1, "LikelyLOH":2, "StrongSomatic":3,
    "LikelySomatic":4, "AFDiff":5, "SampleSpecific":6
}

def get_truth(truth_file):
  truth_vars = set()
  with open(truth_file, 'r') as f:
    for var in f:
      if var[0] == '#': continue
      items = var.split('\t')
      chrom, pos, id, ref, alt, _, filter = items[:7]         
      #if len(chrom) < 6 and filter == "PASS" and (len(ref) > 1 or len(alt) > 1) :
      if (filter.find('PASS') == -1) and (filter.strip() != '.') :
        #print('pass:', filter.strip())
        continue
      #if len(ref) == 1 and len(alt) == 1: continue
      #if len(ref) != 1 or len(alt) != 1: continue
      if len(chrom) < 6:
        site = chrom + ":" + pos + ":" + ref.upper() + ":" + alt.upper()
        truth_vars.add(site)
    return truth_vars

def get_label(truth_vars, key):
  if key in truth_vars:
    return 1
  else :
    return 0
  return 0

def depth_normalization(data, tdepth, ndepth):
    #data['Var1AltFreq']    = data['Var1AltDepth'] / data['Var1Depth']
    data['Var1RefFwdFreq'] = data['Var1RefFwdReads'] / (data['Var1Depth'] - data['Var1AltDepth'])
    data['Var1AltFwdFreq'] = data['Var1AltFwdReads'] / data['Var1AltDepth']
    #data['Var2AltFreq']    = data['Var2AltDepth'] / data['Var2Depth']
    data['Var2RefFwdFreq'] = data['Var2RefFwdReads'] / (data['Var2Depth'] - data['Var2AltDepth'])
    data['Var2AltFwdFreq'] = data['Var2AltFwdReads'] / data['Var2AltDepth']
    data = data.replace([np.inf, -np.inf], np.nan).fillna(-1)
    data['Var1DepthFreq'] = data['Var1Depth'] / tdepth
    data['Var2DepthFreq'] = data['Var2Depth'] / ndepth
    return data 

truth_vars = set()
def get_bed_sorted(bed):
  regions = dict()
  with open(bed, 'r') as f:
    for line in f:
      chr, start, end = line.split('\t')[:3]
      if str(chr) in regions:
        regions[str(chr)].append([int(start), int(end)])
      else:
        regions[str(chr)] = list()
        regions[str(chr)].append([int(start), int(end)])
  for k, v in regions.items():
    regions[k] = sorted(v)
  return regions

def is_in_highconf(chr, start, end, regions):
  regs = regions[chr]
  #print(len(regs))
  #print('processin: ', chr + ":" + str(start) + ":" + str(end))
  for rstart, rend in regs:
    if end < rstart:
      break
    if start > rstart and start < rend:
      #print('in region:', start, rstart, rend)
      return True
  return False

def filter_by_region(data, highconf_regions):
  data['in_highconf'] = data.apply(lambda x: is_in_highconf(x['Chr'], x['Start'], x['End'], highconf_regions), axis = 1)
  return data[data['in_highconf'] == True]

def rf_reader(args):
  cr = list()
  raw = list()
  in_file = args.in_file
  time_step0 = time.time()
  #cr = pd.read_csv(in_file, delimiter = '\t', header = None, engine = 'c', skipinitialspace = True, na_filter = False)
  cr = datautil.get_data_fromtxt(in_file, args.var_type)
  #cr.columns = [*som_features, 'None'] #TODO: i should change the code of c++ to avoid the None colum
  cr['VarLabel'] = cr['VarLabel'].map(varLabel_to_label)
  feature_in = features.som_rf_snv_input_features
  if args.var_type == "INDEL":
    cr['VarType'] = cr['VarType'].map(type_to_label)
    cr['RefLength'] = cr['Ref'].str.len()
    cr['AltLength'] = cr['Alt'].str.len()
    feature_in = features.som_rf_indel_input_features
  time_step1 = time.time()
  print("load time: ", time_step1 - time_step0)
  #print('start of load:', cr.columns, len(cr))
  #hard filter
  cr = hard_filter(cr)
  print('after hard filter:', len(cr))
  #--------filter by confidnece region-----------#
  if args.highconf != "":
    print('filter by region')
    bed = args.highconf
    highconf_regions = get_bed_sorted(bed)
    cr = filter_by_region(cr, highconf_regions)
    print('after filter by region:', cr.columns, len(cr))
  #calc depth
  tdepth = args.tdepth
  ndepth = args.ndepth
  if os.path.exists(args.in_file+'.info'):
    with open(args.in_file + ".info", "r") as f:
      tdepth = round(float(f.readline().strip()))
      ndepth = round(float(f.readline().strip()))
      print("find .info file!, tdepth: {}, ndepth: {}".format(tdepth, ndepth))
  #depth normalization
  #cr = depth_normalization(cr, tdepth, ndepth)
  #print('after depth normalization:', cr.columns, len(cr))
  #--------get truth set-----------#
  global truth_vars
  truth_vars = get_truth(args.truth_file)
  print('length of truth vars is:', len(truth_vars))
  def rabv_get_label(x):
    return get_label(truth_vars, str(x['Chr']) + ':' + str(x['Start']) + ':' + str(x['Ref']).upper() + ':' + str(x['Alt']).upper())

  cr['Label'] = cr.apply(rabv_get_label, axis=1)
  print("truth data: {}, false data: {}, total data: {}".format(sum(cr['Label'] == 1), sum(cr['Label'] == 0), len(cr)))
  #------- split data to train and test -------#
  test_set = cr[cr['Chr'] == 'chr2']
  train_set  = cr[cr['Chr'] != 'chr2']
  train_set[[*feature_in, 'Label']].to_csv(args.tsv, mode='a', header=False, index = False)
  test_set[[*feature_in, 'Label']].to_csv(args.tsv+".test", mode='a', header=False, index = False)
  #cr[[*feature_in, 'Label']].to_csv(args.tsv, mode='a', header=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "train your network")
    parser.add_argument('--in_file', help = "RabbitVar intermidiate file(with fisher test)", type=str, required = True)
    parser.add_argument('--truth_file', help = "truth file / the ground truth(.vcf)", type=str, required = False)
    parser.add_argument('--var_type', help = "var type you want to train(SNV/INDEL)", type=str, required = True)
    parser.add_argument('--tsv', help = "tsv file, which contained data and label", type=str, required = False)
    parser.add_argument('--tdepth', help = "tdepth", type=int, required = False, default=0)
    parser.add_argument('--ndepth', help = "ndepth", type=int, required = False, default=0)
    parser.add_argument('--highconf', help = "high confidence regions", type = str, required = False, default="")
    args = parser.parse_args()
    rf_reader(args)

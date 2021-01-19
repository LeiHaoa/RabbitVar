import os
import joblib
import numpy as np
import pandas as pd
import time 
import argparse
import subprocess
import pybedtools
import math
import re
from cmdparser import *

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

def write_header(fout):
  header = """\
##fileformat=VCFv4.3
##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant Type: SNV Insertion Deletion Complex">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="Chr End Position">
##INFO=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand Bias Info">
##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Reference depth by strand">
##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Variant depth by strand">
##INFO=<ID=PMEAN,Number=1,Type=Float,Description="Mean position in reads">
##INFO=<ID=PSTD,Number=1,Type=Float,Description="Position STD in reads">
##INFO=<ID=QUAL,Number=1,Type=Float,Description="Mean quality score in reads">
##INFO=<ID=QSTD,Number=1,Type=Float,Description="Quality score STD in reads">
##INFO=<ID=SBF,Number=1,Type=Float,Description="Strand Bias Fisher p-value">
##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description="Strand Bias Odds ratio">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mean Mapping Quality">
##INFO=<ID=SN,Number=1,Type=Float,Description="Signal to noise">
##INFO=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency using only high quality bases">
##INFO=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted AF for indels due to local realignment">
##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="No. of bases to be shifted to 3 prime for deletions due to alternative alignment">
##INFO=<ID=MSI,Number=1,Type=Float,Description="MicroSatellite. > 1 indicates MSI">
##INFO=<ID=MSILEN,Number=1,Type=Float,Description="MicroSatellite unit length in bp">
##INFO=<ID=NM,Number=1,Type=Float,Description="Mean mismatches in reads">
##INFO=<ID=LSEQ,Number=1,Type=String,Description="5' flanking seq">
##INFO=<ID=RSEQ,Number=1,Type=String,Description="3' flanking seq">
##INFO=<ID=GDAMP,Number=1,Type=Integer,Description="No. of amplicons supporting variant">
##INFO=<ID=TLAMP,Number=1,Type=Integer,Description="Total of amplicons covering variant">
##INFO=<ID=NCAMP,Number=1,Type=Integer,Description="No. of amplicons don't work">
##INFO=<ID=AMPFLAG,Number=1,Type=Integer,Description="Top variant in amplicons don't match">
##INFO=<ID=HICNT,Number=1,Type=Integer,Description="High quality variant reads">
##INFO=<ID=HICOV,Number=1,Type=Integer,Description="High quality total reads">
##INFO=<ID=SPLITREAD,Number=1,Type=Integer,Description="No. of split reads supporting SV">
##INFO=<ID=SPANPAIR,Number=1,Type=Integer,Description="No. of pairs supporting SV">
##INFO=<ID=DUPRATE,Number=1,Type=Float,Description="Duplication rate in fraction">
##FILTER=<ID=Bias,Description="Strand Bias">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">
"""
  fout.write(header)

def get_socks():
  for line in str(subprocess.check_output('lscpu')).split('\\n'):
    contex = line.split()
    print(contex)
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
  
  print('----------------------\n', cmd, '\n------------------------\n')
  return cmd

def bed_save_to(regions, path):
  with open(path, 'w') as f:
    for line in regions:
      f.write(str(line))

def split_bed(beds, parts, workspace):
  splited_info = []
  stride = len(beds) // parts
  print(stride)
  for i in range(parts - 1):
    print(i)
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
  print(socks)
  splited_info = list()
  if socks > 1:
    #re-distribure bed file
    #beds = pybedtools.example_bedtool(param['bed'])
    #sorted(beds)
    #print(type(beds))
    beds = list()
    with open(param['bed'], 'r') as f:
      for line in f:
        beds.append(f)
    #split_info:[(bed1, out1), (bed2, out2), ...]
    splited_info = split_bed(beds, socks, workspace) #split file according to numa
    ps = list()
    for i in range(int(socks)):
      cmd = prepare_cmd(BIN, splited_info[i][0], splited_info[i][1], param)
      #ps.append(subprocess.Popen(cmd, stderr=subprocess.DEVNULL))
      #for p in ps:
        #p.wait()
  else:
    bed = param['bed'] if 'bed' in param else None
    out = param['out']
    cmd = prepare_cmd(BIN, bed, out, param)
    print(cmd)
    splited_info.append((bed, out))
    #subprocess.Popen(cmd, stderr=subprocess.DEVNULL)

  print("now, all process run over!")
  return splited_info

def rf_filter(args, in_file):
  cr = list()
  raw = list()
  cr = pd.read_csv(in_file, delimiter = '\t', header = None, engine = 'c', skipinitialspace = True)
  cr.columns = [*som_features, 'None'] #TODO: i should change the code of c++ to avoid the None colum
  cr['VarLabel'] = cr['VarLabel'].map(varLabel_to_label)
  cr['VarType'] = cr['VarType'].map(type_to_label)
  #snv data process 
  time_start = time.time()
  snvs = cr[cr['VarType'] == 0]
  inputs = snvs[[*som_selected_features, "VarLabel"]].to_numpy()
  #clf = joblib.load(args.snv_model)
  #scale = args.scale
  clf = joblib.load("/home/old_home/haoz/workspace/FastVC/RandomForest/models/som_snv_0108.pkl")
  scale = 0.2
  snv_pred = my_predict(clf, inputs, scale)
  snv_result = snvs.loc[snv_pred == 1]
  time_end = time.time()
  print("time processing snv: {} s".format(time_end - time_start))
  
  #indel data process
  time_start = time.time()
  indels = cr[cr['VarType'] != 0]
  indels['RefLength'] = cr['Ref'].map(len)
  indels['AltLength'] = cr['Alt'].map(len)
  inputs = indels[["RefLength", "AltLength", "VarType", *som_selected_features, "VarLabel"]].to_numpy()
  #clf = joblib.load(args.indel_model)
  #scale = args.scale
  clf = joblib.load("/home/old_home/haoz/workspace/FastVC/RandomForest/models/som_indel_0108.pkl")
  scale = 0.2
  indel_pred = my_predict(clf, inputs, scale)
  indel_result = indels.loc[indel_pred == 1]
  time_end = time.time()
  print("time processing indel: {} s".format(time_end - time_start))

  return snv_result, indel_result

def my_predict(clf, data, scale):
    proba = clf.predict_proba(data)
    return np.asarray([1 if x > scale else 0 for x in proba[:,1]])

def format_record(record):
  try:
    [sample, gene, chrt, start, end, ref, alt, dp1, vd1, rfwd1, rrev1, vfwd1, \
     vrev1, gt1, af1, bias1, pmean1, pstd1, qual1, qstd1, mapq1, sn1, hiaf1, \
     adjaf1, nm1, sbf1, oddratio1, dp2, vd2, rfwd2, rrev2, vfwd2, vrev2, \
     gt2, af2, bias2, pmean2, pstd2, qual2, qstd2, mapq2, sn2, hiaf2, \
     adjaf2, nm2, sbf2, oddratio2, shift3, msi, msilen, lseq, rseq, seg, \
     status, vtype, sv1, duprate1, sv2, duprate2, pvalue, oddratio]  = record[:61]
  except ValueError:
    print("invalide record: \n", record, "\n record length ---> ", len(record))
    exit(-1)
  
  rd1 = rfwd1 + rrev1
  rd2 = rfwd2 + rrev2
  if vtype == "" : 
    vtype = "REF"
  GTFREQ = 0.2
  FREQ = 0.02
  gt  = "1/1" if (1-af1 < GTFREQ) else ("1/0" if af1 >= 0.5 else ("0/1" if af1 >= FREQ else "0/0"))
  gtm = "1/1" if (1-af2 < GTFREQ) else ("1/0" if af2 >= 0.5 else ("0/1" if af2 >= FREQ else "0/0"))
  bias1 = bias1.replace(';',',')
  bias2 = bias2.replace(';',',')
  if bias1 == "0": bias1 = "0,0" 
  if bias2 == "0": bias2 = "0,0" 
  mapq1 = f'{mapq1:.0}'
  mapq2 = f'{mapq2:.0}'
  qual = int(math.log(vd1)/math.log(2) * qual1) if vd1 > vd2  else int(math.log(vd2)/math.log(2) * qual2)

  pinfo1 = "\t".join([chrt, str(start), ".", ref, alt, str(qual)])
  filters = "PASS"
  sample_nowhitespace = re.sub(r'\s', '_', sample)

  pinfo2_1 = "STATUS={};SAMPLE={};TYPE={};DP={};VD={};AF={:.6f};SHIFT3={};MSI={};MSILEN={};SSF={};SOR={};LSEQ={};RSEQ={}".format(status, sample_nowhitespace, vtype, dp1, vd1, af1, shift3, msi, msilen, pvalue, oddratio, lseq, rseq)
  pinfo2_2 = "{}:{}:{}:{}:{}:{}:{:.6f}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{:.6f}:{:.6f}:{}".format(gt, dp1, vd1, str(vfwd1) + "," + str(vrev1), str(rfwd1)+","+str(rrev1), str(rd1)+","+str(vd1), af1, bias1, pmean1, pstd1, qual1, qstd1, sbf1, oddratio1, mapq1, sn1, hiaf1, adjaf1, nm1) 
  pinfo2_3 = "{}:{}:{}:{}:{}:{}:{:.6f}:{}:{}:{}:{}:{}:{}:{}:{}:{}:{:.6f}:{:.6f}:{}".format(gtm, dp2, vd2, str(vfwd2) + "," + str(vrev2), str(rfwd2)+","+str(rrev2), str(rd2)+","+str(vd2), af2, bias2, pmean2, pstd2, qual2, qstd2, sbf2, oddratio2, mapq2, sn2, hiaf2, adjaf2, nm2)
  pinfo2 = "\t".join([pinfo2_1, "GT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM", pinfo2_2, pinfo2_3])
  return "\t".join([pinfo1, filters, pinfo2])

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "rabbitvar")
  detector_parser = parser.add_argument_group("detector_parser")
  detectorParam(detector_parser)
  rabbitvar_parser = parser.add_argument_group("rabbitvar_parser")
  rabbitvarParam(rabbitvar_parser)
  
  args = parser.parse_args()
  detector_param = {}
  for x in detector_parser._group_actions:
    k, v = x.dest, getattr(args, x.dest, None)
    if v:
      detector_param[k] = v
  print(detector_param)
  splited_info = run_rabbitvar(args.BIN, args.workspace, detector_param)
  exit(0)
  #if keep_intermident:
  #  interm_file = 1
  #else:
  #  for itf in os.path.join(tmpdir): # multi thread maybe better
  #    vcf.append(rf_filter(itf))

  #tmp_file = "/home/old_home/haoz/workspace/FastVC/detection_result/FD_DATASET/demo.txt"
  #vcf_file = "./tmpresult.vcf"
  vcf_file = args.vcf
  vcf = list()
  for detector_out in splited_info:
    snvs, indels = rf_filter(None, detector_out[1])
    print(type(snvs),type(indels))
    vcf.append(snvs)
    vcf.append(indels)
    print(len(vcf), len(vcf[0]))
    #write vcf file
    with open(vcf_file, 'w') as f:
      write_header(f)
      tmp = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", args.Name]) #TODO what if user not specified a sample name???
      tmp += '\n'
      f.write(tmp)
      for pddata in vcf:
        for i, record in pddata.iterrows():
          f.write(format_record(record) + '\n')

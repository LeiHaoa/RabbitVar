#! /usr/bin/python3

from features import *
import pandas as pd

def get_gt(af):
    GTFREQ = 0.2
    FREQ  = 0.02
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

def format_indel_data_item(jri, fisher = True):
    #[FIXME] fisher should always be true, otherwish the map is wrong
    data = list()
    # key is chrom:pos like "chr1:131022:A:AT"
    key = jri[2] + ":" + jri[3] + ":" + jri[5] + ":" + jri[6]
    data.append(len(jri[fe2i["Ref"]])) #refallele len
    data.append(len(jri[fe2i["Alt"]])) #varallel len
    for sf in fvc_sf:
        data.append(jri[fe2i[sf]])
     
    if fisher:
        data.extend(varLabel_to_label[jri[fe2i["VarLabel"]]])#varlabel
        data.extend(indels_label[jri[fe2i["VarType" ]]]) #vartype
    else:
        print("not support to train if you not run rabbitvar without --fiser!!")
        exit(-1)
        data.extend(varLabel_to_label[jri[fe2i["VarLabel"]]])#varlabel
        data.extend(indels_label[jri[fe2i["VarType"]]])

    if len(data) != SOM_INDEL_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        exit(-1)
    return key, data

def format_snv_data_item(jri, fisher = True):
    if not fisher:
        print("not support to train if you not run rabbitvar without --fiser!!")
        exit(-1)
    data = list()
    # key is chrom:pos like "chr1:131022:A:T"
    key = jri[2] + ":" + jri[3] + ":" + jri[5] + ":" + jri[6] #TODO: case sensitive
    for sf in fvc_sf:
        data.append(jri[fe2i[sf]])
    data.extend(varLabel_to_label[jri[fe2i["VarLabel"]]])#varlabel
    if len(data) != SOM_SNV_FEATURES:
        print("fvc data length error: \n", len(data), data, " ori\n", jri)
        exit(-1)
    return key, data

def get_indel_data(fvc_result_path):
    #--- read fastvc result file and format ---#
    fastvc_indel_dict = dict()
    index = -1
    with open(fvc_result_path, 'r') as f:
        for line in f:
            index += 1
            items = line.strip().split("\t")
            if items[fe2i['VarType']] not in indels_label:
                continue
            if len(items) == 55:
                k, d = format_indel_data_item(items, False)
                fastvc_indel_dict[k] = [d, index]
            elif len(items) == 61 :
                k, d = format_indel_data_item(items, True)
                fastvc_indel_dict[k] = [d, index]
            else:
                print("your train file should be 55 or 61 items! but you have {} items".format(len(items)))
    print("get fastvc indels data done: ", len(fastvc_indel_dict))
    return fastvc_indel_dict

def get_snv_data(fvc_result_path):
    fastvc_snv_dict = dict()
    index = -1
    with open(fvc_result_path, 'r') as f:
        for line in f:
            index += 1
            items = line.strip().split("\t")
            if items[fe2i['VarType']] != snv_label:
                continue
            if len(items) == 55:
                k, d = format_snv_data_item(items, False)
                fastvc_snv_dict[k] = [d, index]
            elif len(items) == 61 :
                k, d = format_snv_data_item(items, True)
                fastvc_snv_dict[k] = [d, index]
    print("get fastvc SNV data done: ", len(fastvc_snv_dict))
    return fastvc_snv_dict

def get_data_fromtxt_manual(fvc_result_path, vtype = 'SNV'):
    if vtype.upper() == 'SNV':
        return get_snv_data(fvc_result_path)
    elif vtype.upper() == 'INDEL':
        return get_indel_data(fvc_result_path)
    else:
        print("unrecognized variant type: {} !".format(vtype))
        exit(-1)

def get_data_fromtxt(fvc_result_path, vtype = 'ALL'):
  dtype_dictionary = {'Sample':'str','Gene':'str','Chr':'str','Start':'int64','End':'int64','Ref':'str',
                      'Alt':'str','Var1Depth':'int32','Var1AltDepth':'int32','Var1RefFwdReads':'int32',
                      'Var1RefRevReads':'int32','Var1AltFwdReads':'int32','Var1AltRevReads':'int32',
                      'Var1Genotype':'str','Var1AF':'float32','Var1Bias':'str','Var1PMean':'float32',
                      'Var1PStd':'int8','Var1QMean':'float32','Var1QStd':'int8','Var1MQ':'float32',
                      'Var1Sig_Noise':'float32','Var1HiAF':'float32','Var1ExtraAF':'float32','Var1NM':'float32',
                      'Var1Pvalue':'float32','Var1Oddr':'float32','Var2Depth':'int32','Var2AltDepth':'int32',
                      'Var2RefFwdReads':'int32','Var2RefRevReads':'int32','Var2AltFwdReads':'int32',
                      'Var2AltRevReads':'int32','Var2Genotype':'object','Var2AF':'float32','Var2Bias':'object',
                      'Var2PMean':'float32','Var2PStd':'int8','Var2QMean':'float32','Var2QStd':'int8',
                      'Var2MQ':'float32','Var2Sig_Noise':'float32','Var2HiAF':'float32','Var2ExtraAF':'float32',
                      'Var2NM':'float32','Var2Pvalue':'float32','Var2Oddr':'float32','shift3':'int32',
                      'MSI':'float32','MSI_NT':'int32','5pFlankSeq':'object','3pFlankSeq':'object',
                      'Seg':'object','VarLabel':'object','VarType':'object','Duprate1':'float32',
                      'SV_info1':'int32','Duprate2':'float32','SV_info2':'int32',
                      'Pvalue':'float32','Oddratio':'float32','None':'float32'}
  cr = pd.read_csv(fvc_result_path, delimiter = '\t', names = [*som_features, 'None'], dtype = dtype_dictionary, header = None, engine = 'pyarrow', na_values=[''])
  cr.columns = [*som_features, 'None'] #TODO: i should change the code of c++ to avoid the None colum
  if vtype == 'INDEL':
      return cr[cr['VarType'] != 'SNV']
  elif vtype == 'SNV':
      return cr[cr['VarType'] == 'SNV']
  elif vtype == 'ALL':
      return cr
  else:
      print('ERROR: unsupported variant type: ', vtype)
      exit(-1)

def get_data_fromcsv(data_path, header, vtype = 'SNV'):
  data = pd.read_csv(data_path, header=None)
  if vtype.upper() == 'INDEL':
    data.columns = header
    #data = data[data['VarLabel'] != 0] #vartype != germline
  elif vtype.upper() == 'SNV':
    data.columns = header
    #data = data[data['VarLabel'] != 0] #vartype != germline
  return data 

def hard_filter(data):
  hard_filtered_data = data[~((data['Var1AF'] < 0.01)
                              | (data['Var1QMean'] <= 20)
                              | (data['Var1NM'] >= 6)
                              | (data['Var1NM'] < 0)
                              | ((data['Var1Depth'] == 0) | (data['Var2Depth'] == 0))
                              | ((data['VarLabel'] != 3) & (data['VarLabel'] != 4) & (data['VarLabel'] != 5)))]
  return hard_filtered_data

def hard_filter_keeporg(data):
  #hard flag: 1 means it should filter by hard filter, 0 means it should keep
  data['hard_flag'] = ((data['Var1AF'] < 0.01)
                              | (data['Var1QMean'] <= 20)
                              | (data['Var1NM'] >= 6)
                              | (data['Var1NM'] < 0)
                              | ((data['Var1Depth'] == 0) | (data['Var2Depth'] == 0))
                              | ((data['VarLabel'] != 3) & (data['VarLabel'] != 4) & (data['VarLabel'] != 5)))
  return data

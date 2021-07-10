#! /usr/bin/python3
from features import *
import pandas as pd

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

def get_data_fromtxt(fvc_result_path, vtype = 'SNV'):
    if vtype.upper() == 'SNV':
        return get_snv_data(fvc_result_path)
    elif vtype.upper() == 'INDEL':
        return get_indel_data(fvc_result_path)
    else:
        print("unrecognized variant type: {} !".format(vtype))
        exit(-1)

def get_data_fromcsv(data_path, header, vtype = 'SNV'):
  print(data_path, header, len(header))
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
                              | ((data['VarLabel'] != 3) & (data['VarLabel'] != 4)))]
  return hard_filtered_data

def hard_filter_keeporg(data):
  #hard flag: 1 means it should filter by hard filter, 0 means it should keep
  data['hard_flag'] = ((data['Var1AF'] < 0.01)
                              | (data['Var1QMean'] <= 20)
                              | (data['Var1NM'] >= 6)
                              | (data['Var1NM'] < 0)
                              | ((data['VarLabel'] != 3) & (data['VarLabel'] != 4)))
  print(data['hard_flag'])
  return data

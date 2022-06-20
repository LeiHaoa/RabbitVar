#!/usr/bin/env python3
import pandas as pd
import sys
import os
sys.path.append('/home/user_home/haoz/workspace/RabbitVar/utils')
import features
import numpy as np
from sklearn.utils import shuffle

VTYPE="SNV"

def get_data_fromcsv(data_path, columns, vtype = "SNV"):
  data = pd.read_csv(data_path, header=None, names=columns, engine = 'pyarrow')
  return data 

tsv = "./data_SNV_HCC1395_WGS_Spike.tsv"
data1 = get_data_fromcsv(tsv, [*features.som_rf_snv_input_features, 'Label'], VTYPE)
truth_data = data1[data1['Label'] == 1]
false_data = data1[data1['Label'] == 0]
print(f"tsv dataset1: {tsv}, truth: {len(truth_data)}, false: {len(false_data)}")

'''
truth_data['is_selected'] = np.random.uniform(0, 1, len(truth_data)) <= (29480 / 161449)
false_data['is_selected'] = np.random.uniform(0, 1, len(false_data)) <= (2901586 / 3978302)
print(sum(truth_data['is_selected'] == 1), sum(false_data['is_selected'] == 1))
truth_data = truth_data[truth_data['is_selected'] == 1]
false_data = false_data[false_data['is_selected'] == 1]
'''

#---------------all HCC data--------------#
tsv = "./data_SNV_DREAMC.tsv"
data2 = get_data_fromcsv(tsv, [*features.som_rf_snv_input_features, 'Label'], VTYPE)
truth_data = data2[data2['Label'] == 1]
false_data = data2[data2['Label'] == 0]
print(f"tsv dataset1: {tsv}, truth: {len(truth_data)}, false: {len(false_data)}")

#aug_data = shuffle(pd.concat([HCC1395_data, truth_data, false_data], axis = 0))
aug_data = shuffle(pd.concat([data1, data2], axis = 0))
aug_data = aug_data[[*features.som_rf_snv_input_features, 'Label']]
aug_data.to_csv("./data_SNV_SPICKEIN_DREAM.tsv",header = False, index = False)

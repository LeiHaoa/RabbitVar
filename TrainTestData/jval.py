
import sys
import os
sys.path.append('/home/user_home/haoz/workspace/RabbitVar/utils')

import joblib
import numpy as np
import pandas as pd
import sklearn.datasets
import sklearn.metrics
import sklearn.model_selection
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import shuffle
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import f1_score
import argparse

from features import *
from datautil import get_data_fromcsv, hard_filter
import features
import datautil


test_data = get_data_fromcsv("./data_INDEL_HCC1395_ALL.tsv.test", columns=[*features.som_rf_indel_input_features, "Label"], vtype = "INDEL")
#print(test_data.head())
test_data = hard_filter(test_data)
clf = joblib.load("./models/train_INDEL_HCC1395_ALL.gpu.pkl")
test_set = test_data[features.som_rf_indel_input_features]
ytest = test_data['Label']
ypred = clf.predict(test_set)

res = f1_score(ytest, ypred, average=None)
print(res)
#print("Accuracy score", sklearn.metrics.accuracy_score(ytest, ypred))

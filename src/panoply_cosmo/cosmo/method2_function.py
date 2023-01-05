#!/usr/bin/env python
# coding: utf-8

# -*- coding: utf-8 -*-
"""2B_modelAB.ipynb
"""
# Commented out IPython magic to ensure Python compatibility.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib inline
import seaborn as sns

sns.set()
import os

pd.set_option('display.max_rows', 100)

from sklearn.model_selection import StratifiedKFold
from sklearn import linear_model

from sklearn.linear_model import LinearRegression
from sklearn.linear_model import LogisticRegressionCV

from sklearn import svm
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

from sklearn.neighbors import KNeighborsRegressor
from sklearn.metrics import r2_score

import random
import math

from sklearn.decomposition import PCA

from sklearn.linear_model import LassoCV

from sklearn import preprocessing

#import pandas_profiling

from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import GenericUnivariateSelect, SelectFdr

import gc

from sklearn.feature_selection import SelectFromModel
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import f_classif

from sklearn.svm import SVC
from sklearn.svm import LinearSVC

from sklearn.neighbors import KNeighborsClassifier

from sklearn.model_selection import LeaveOneOut

from sklearn.ensemble import VotingClassifier

from sklearn.preprocessing import RobustScaler

from random import seed
from random import randint

import warnings

import argparse
import sys

warnings.filterwarnings('ignore')

rand_times = 10


def ensemble(model, X, Y, Y_pred):
    # loo = LeaveOneOut()
    seed(1)
    for i in range(rand_times):
        print("Repeat %d" % (i))
        s = randint(0, 10 ** 6)
        loo = StratifiedKFold(n_splits=5, random_state=s, shuffle=True)
        split_num = loo.get_n_splits(X)

        split_idx = 0
        for train_index, test_index in loo.split(X, Y):
            X_train, X_test = X[train_index], X[test_index]
            Y_train, Y_test = Y[train_index], Y[test_index]

            model.fit(X_train, Y_train)

            # X_test_pred  = model.predict_proba(X_test)[:,1]
            # for j in range(len(test_index)):
            #     Y_pred[test_index[j],i] = X_test_pred[j]
            Y_pred[test_index, i] = model.predict_proba(X_test)[:, 1]

            # print(Y_pred[:,i][test_index])

    return Y_pred

def run(d1_file:str, d2_file:str, sample_file:str, sample_labels, out_dir="./", prefix="test"):

    print("Processing file: %s and %s\n" % (d1_file, d2_file) )

    test_d1 = pd.read_csv(d1_file, sep='\t')
    test_d1 = test_d1.T.fillna(0)

    test_d2 = pd.read_csv(d2_file, sep='\t')
    test_d2 = test_d2.T.fillna(0)

    test_cli = pd.read_csv(sample_file, sep='\t')

    # to be used features/genes, remove highly correlated ones
    corr_matrix = test_d1.loc[:, test_d1.columns[test_d1.nunique() > 1].tolist()].corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))
    use_pro = [column for column in upper.columns if any(upper[column] <= 0.9)]

    corr_matrix = test_d2.loc[:, test_d2.columns[test_d2.nunique() > 1].tolist()].corr().abs()
    upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))
    use_rna = [column for column in upper.columns if any(upper[column] <= 0.9)]

    del corr_matrix
    del upper
    gc.collect()

    clf1 = LogisticRegressionCV(cv=5, penalty='l1', class_weight='balanced',
                                solver='liblinear',
                                max_iter=5000)
    clf3 = Pipeline([
        ('feature_selection', GenericUnivariateSelect(f_classif, 'fdr', param=0.1)),
        ('classification', SVC(gamma='auto', C=1.3, kernel="rbf",
                               probability=True, class_weight='balanced',
                               max_iter=1000))
    ])
    clf6 = LogisticRegressionCV(cv=5, penalty='l2', class_weight='balanced')
    clf8 = Pipeline([
        ('feature_selection', GenericUnivariateSelect(f_classif, 'fdr', param=0.1)),
        ('classification', KNeighborsClassifier(n_neighbors=20,
                                                weights='distance',
                                                ))
    ])
    model = VotingClassifier(estimators=[
        ('lr', clf1),
        #         ('lgbm', clf2),
        #             ('svc', clf3),
        #         ('lsvn', clf4),
        #         ('knn', clf5),
        ('lr_l2', clf6),
        #         ('nsc', clf7),
        # ('knn_fdr', clf8),
    ],
        voting='soft',
        n_jobs=1)

    #labels = [l for l in test_cli.columns[test_cli.nunique() > 1] if l != "sample" and l != "age"]

    ##Protein Build Model A
    pro_labels = sample_labels # ["gender", "msi"]
    test_d1_data = test_d1.merge(test_cli, left_index=True, right_on="sample").set_index("sample")

    pred_pro = {}
    X = test_d1_data[use_pro].values
    scaler = StandardScaler()
    scaler.fit(X)
    X = scaler.transform(X)

    les = {}
    for l in pro_labels:
        le = preprocessing.LabelEncoder()
        le.fit(test_d1_data[l].values)
        Y = le.transform(test_d1_data[l].values)
        les[l] = le

        pred_res = np.zeros((X.shape[0], rand_times))
        ensemble(model, X, Y, pred_res)
        pred_pro[l] = pred_res

    del X, Y, scaler
    gc.collect()

    modelA_dist = pd.DataFrame(index=test_d1_data.index)
    for l in pro_labels:
        modelA_dist[l + "_dist"] = np.absolute(np.sum(pred_pro[l], axis=1) / rand_times
                                               - les[l].transform(test_cli[l]))
        modelA_dist[l + "_prob"] = np.sum(pred_pro[l], axis=1) / rand_times

        ##RNA Build Model B

    rna_labels = sample_labels #["gender", "msi"
                  #               , "stage", "colon_rectum"
                  #]
    test_d2_data = test_d2.merge(test_cli, left_index=True, right_on="sample").set_index("sample")

    pred_rna = {}
    X = test_d2_data[use_rna].values
    scaler = StandardScaler()
    scaler.fit(X)
    X = scaler.transform(X)

    les = {}
    for l in rna_labels:
        le = preprocessing.LabelEncoder()
        le.fit(test_d2_data[l].values)
        Y = le.transform(test_d2_data[l].values)
        les[l] = le

        pred_res = np.zeros((X.shape[0], rand_times))
        ensemble(model, X, Y, pred_res)
        pred_rna[l] = pred_res

    del X, Y, scaler
    gc.collect()

    modelB_dist = pd.DataFrame(index=test_d2_data.index)
    for l in rna_labels:
        modelB_dist[l + "_prob"] = np.sum(pred_rna[l], axis=1) / rand_times
        modelB_dist[l + "_dist"] = np.absolute(np.sum(pred_rna[l], axis=1) / rand_times -
                                               les[l].transform(test_cli[l]))

    modelA_dist.to_csv(out_dir + "/" + str(prefix) + "_ModelA_result.csv")
    modelB_dist.to_csv(out_dir + "/" + str(prefix) + "_ModelB_result.csv")

def main():
    parser = argparse.ArgumentParser(description='COSMO')
    parser.add_argument('-d1', '--dataset1', default=None, type=str, required=True,
                        help="Quantification data at gene level, for example, proteomics data")
    parser.add_argument('-d2', '--dataset2', default=None, type=str,
                        help="Quantification data at gene level, for example, RNA-Seq data")
    parser.add_argument('-s', '--sample_file', default=None, type=str, required=True,
                       help="RNA expression data")
    parser.add_argument('-l', '--sample_label', default=None, type=str,
                        help="Sample class labels for prediction")
    parser.add_argument('-o', '--out_dir', default="./", type=str,
                        help="Output directory")
    parser.add_argument('-prefix', '--prefix', default="test", type=str,
                        help="The prefix of output file(s)")

    args = parser.parse_args(sys.argv[1:len(sys.argv)])

    d1_file = args.dataset1
    d2_file = args.dataset2
    sample_file = args.sample_file
    sample_label = args.sample_label.split(",")
    out_dir = args.out_dir
    prefix = args.prefix

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    run(d1_file=d1_file, d2_file=d2_file, sample_file=sample_file, sample_labels=sample_label, out_dir=out_dir)

if __name__=="__main__":
    main()


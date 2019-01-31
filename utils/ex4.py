import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.tree import *
from sklearn.metrics import confusion_matrix
import math

import catboost

import lightgbm as lgb
import csv

import xgboost
from sklearn.model_selection import  GridSearchCV

import datetime
import pandas as pd
from dateutil.parser import parse
from sklearn.ensemble import RandomForestClassifier

from sklearn.decomposition import PCA

from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report
import random

import os

def dt2epoch(value):
    if type(value) is not str:
        return value

    d = parse(value)
    epoch = (d - datetime.datetime(1970, 1, 1)).total_seconds()
    return epoch


data = pd.read_csv('~/mnt/ssd/kaggle-talkingdata2/competition_files/train_sample.csv')
parse_dates = ['attributed_time', 'click_time']
# train['attributed_time']=pd.to_numeric(train['attributed_time'])
# train['click_time']=pd.to_numeric(train['click_time'])

for cur in parse_dates:
    data[cur] = data[cur].apply(dt2epoch)
print(data.info())
print(data.shape)

data['is_attributed'].value_counts()


def down_sample(data, sample_rate, random_seed):
    '''
    This function returns a down-sampled version of the dataset. It will
    down sample the 0 class to have sample_rate * (number of 1-class instances).
    '''

    zero_counter = data['is_attributed'].value_counts()[0]
    one_counter = data['is_attributed'].value_counts()[1]
    np.random.seed(random_seed)
    n_downsamples = int(float(sample_rate) * one_counter)
    print("n_downsamples: {}/{}. sample_rate {}".format(n_downsamples, zero_counter, sample_rate))
    choice = np.random.choice(zero_counter, min(n_downsamples, zero_counter), replace=False)
    #   print("choice: {}".format(choice))
    downsampling_0 = data.loc[data["is_attributed"] == 0].iloc[choice, :]
    #   print(downsampling_0.info())
    downsampling_1 = data.loc[data["is_attributed"] == 1]
    #   print(downsampling_1.info())
    reduced_ds = pd.concat((downsampling_0, downsampling_1))
    #   print(reduced_ds.info())
    return reduced_ds


def bagging_on_samples(data, sample_rate, random_seed, classifier, n):
    '''
    This function returns a list of n trained classifiers, each trained
    on a down sampled version of the data. To do that, you should call
    down_sample n times, each time with a different seed.
    The argument classifier is a sklearn classifier, already initialized with
    the relevant arguments.
    '''
    fitted_models = []
    for i in np.arange(n):
        downsampling = down_sample(data, sample_rate, random_seed + i)

        fitted_models.append(classifier.__class__(**classifier.get_params()).fit(downsampling.drop(["is_attributed"], axis=1),
                                            downsampling["is_attributed"]))

    return fitted_models


def prepare_data1(data):
    data_copy = data.copy()  # .iloc[:1000,:]

    '''
    This function adds a few more advanced features to the dataset.
    '''

    tss = data_copy.dropna().groupby(by=['ip'])['attributed_time'].agg("min")

    df_ip = data_copy.groupby(by=['ip'])['click_time'].apply(lambda x: sorted(list(x)))

    #   print(data.iloc[0].loc['click_time'])
    for i, row in data_copy.iterrows():
        l_click_time_head = df_ip[row['ip']][:df_ip[row['ip']].index(row['click_time'])]
        data_copy.loc[i, "nclh"] = len(l_click_time_head) - np.searchsorted(l_click_time_head,
                                                                            row['click_time'] - 1000 * 60 * 60,
                                                                            side='left')
        data_copy.loc[i, "ncld"] = len(l_click_time_head) - np.searchsorted(l_click_time_head,
                                                                            row['click_time'] - 1000 * 60 * 60 * 24,
                                                                            side='left')

        # if i % 1000 == 0:
        #     print(i)

    data_copy['tss'] = data.apply(lambda row: int(row['attributed_time'] > (tss[row['ip']])) \
        if row['ip'] in tss and row['attributed_time'] != pd.NaT else 0, axis=1)

    return data_copy[["ip", "app", "device", "os", "channel", "click_time", "tss" ,"ncld", "nclh", "is_attributed"]]



def prepare_data2(data):
    '''
    This function adds an advanced LDA-based app groups as another feature.
    '''

    # <Your code goes here>


def prepare_data3(data):
    '''
    This function returns a one-hot encoded version of the data.
    '''

    return pd.get_dummies(data)


def grid_search(classifier, arguments, data, n_fold):
    '''
    Returns the an instance of the classifier, initialized with the best configuration.
    n_fold is the number of cross validation folds to use. You may use GridSearchCV.
    '''

    clf = GridSearchCV(classifier, arguments, cv=n_fold)
    clf.fit(data.drop(["is_attributed"], axis=1), data["is_attributed"])
    return classifier.__class__(**clf.best_params_)


n_fold=5
even_data = down_sample(data, 1, 1)
data_preprocessed_1 = prepare_data1(even_data)
data_preprocessed_3 = prepare_data3(data_preprocessed_1)


arguments = {"learning_rate" : [0.01, 0.1, 0.2], "max_depth": [1,3,5,10,15,20], "min_child_weight": [1,5,10,15,20]}
xg_clf=grid_search(xgboost.XGBClassifier(), arguments, data_preprocessed_3, n_fold)

arguments = {"learning_rate" : [0.01, 0.1, 0.2], "depth": [1,3,5],"cat_features" : [[0,1,2,3,4]]}
cat_clf=grid_search(catboost.CatBoostClassifier(), arguments, data_preprocessed_1, n_fold)

arguments = {"learning_rate" : [0.01, 0.1, 0.2], "max_depth": [1,3,5], "min_data_in_leaf": [10,20,30], "num_leaf":[5,10,20], "categorical_feature" : [[0,1,2,3,4]]}
lgb_clf=grid_search(lgb.LGBMClassifier(), arguments, data_preprocessed_1, n_fold)



downsampled=down_sample(data, 5, 1)
downsampled_p1=prepare_data1(downsampled)
downsampled_p3=prepare_data3(downsampled_p1)

choice=np.random.choice(downsampled_p1.shape[0],int(downsampled_p1.shape[0]*0.2),replace=False)
train=downsampled_p1.loc[~np.isin(np.arange(downsampled_p1.shape[0]),choice)]
X_train=train.drop("is_attributed", axis=1)
y_train=train["is_attributed"]

test=downsampled_p1.iloc[choice]
X_test=test.drop("is_attributed", axis=1)
y_test=test["is_attributed"]

train_p3=prepare_data3(train)
X_train_p3=prepare_data3(X_train)
y_train_p3=y_train

X_test_p3=prepare_data3(X_test)
y_test_p3=y_test



forest=RandomForestClassifier(n_estimators=100).fit(downsampled_p3.drop("is_attributed", axis=1),downsampled_p3["is_attributed"])
importances = forest.feature_importances_
indices = np.argsort(importances)[::-1]

# Print the feature ranking
print("Feature ranking:")

for i in np.arange(X_train.shape[1]):
    print("{}. feature {} ({}): {}".format(i + 1, indices[i], X_train.columns.values[indices[i]] , importances[indices[i]]))



def predict_with_bagging(classifiers, data):
    print(data.columns)
    predictions = np.reshape(classifiers[0].predict(data), (1, data.shape[0]))
    for classifier in classifiers[1:]:
        predictions = np.r_[predictions, np.reshape(classifier.predict_proba(data)[:,1], (1, data.shape[0]))]

    #   print(predictions)
    return np.mean(predictions, axis=0)



df_summary=pd.DataFrame()
for n in [2,4,5]:  #
    for sample_rate in [1,2,5]:  # ,2,5
        for classifier_name, classifier in [("xg", xg_clf) , ("lgb", lgb_clf), ("cat", cat_clf)]: # categorical_feature=[0,1,2,3,4]

            random_seed = int(random.random() * 1000)
            print(random_seed)

            if classifier_name == "xg":
                classifiers = bagging_on_samples(train_p3, sample_rate, random_seed, classifier.__class__(**classifier.get_params()), n)
            elif classifier_name == "lgb":
                classifiers = bagging_on_samples(train, sample_rate, random_seed, classifier.__class__(**classifier.get_params()), n)
            else:
                classifiers = bagging_on_samples(train, sample_rate, random_seed, classifier.__class__(**classifier.get_params()), n)

            if classifier_name=="xg":
                prediction = predict_with_bagging(classifiers, X_test_p3)
            elif classifier_name=="lgb":
                prediction = predict_with_bagging(classifiers, X_test)
            else:
                prediction = predict_with_bagging(classifiers, X_test)

            prediction = np.round(prediction)
            # print([x for x in y_test if x != 0])
            # print([x for x in prediction if x != 0])
            accuracy=accuracy_score(y_test, prediction)
            print("{} accuracy".format(classifier_name), accuracy)
            # print(classification_report(y_test, prediction))
            df_summary=df_summary.append({"name": classifier_name, "n" : n, "sample_rate": sample_rate, "accuracy" : accuracy}, ignore_index=True)

print df_summary





classifier_name, classifier = "xg", xg_clf
pca = PCA(n_components=2)
X=downsampled_p3.values
y=downsampled_p3["is_attributed"]
pca.fit_transform(X.astype(np.float64))
fig = plt.figure(1, figsize=(20, 20))
ax = fig.add_subplot(111)
ax.scatter(X[:, 0], X[:, 1], c=y, cmap='jet')
plt.savefig(os.path.join("/home/hag007/Desktop/pca.png"))



classifier_name, classifier = "xg", xg_clf

pca = PCA(n_components=2)
downsampled=down_sample(data, 5, 1)
downsampled_p1=prepare_data1(downsampled)
downsampled_p3=prepare_data3(downsampled_p1)
X=downsampled_p3.copy()
pca.fit_transform(X.astype(np.float64))

choice=np.random.choice(X.shape[0],int(X.shape[0]*0.2),replace=False)
train=X.loc[~np.isin(np.arange(X.shape[0]),choice)]
X_train=train.drop("is_attributed", axis=1)
y_train=train["is_attributed"]

test=X.iloc[choice]
X_test=test.drop("is_attributed", axis=1)
y_test=test["is_attributed"]

train_p3=prepare_data3(train)
X_train_p3=prepare_data3(X_train)
y_train_p3=y_train
X_test_p3=prepare_data3(X_test)
y_test_p3=y_test

clf=classifier.__class__(**classifier.get_params()).fit(X_train_p3,y_train_p3)
prediction = clf.predict(X_test_p3)
accuracy = accuracy_score(y_test_p3, prediction)
print("{} accuracy".format(classifier_name), accuracy)



#####################

classifier_name, classifier = "xg", xg_clf

pca = PCA(n_components=3)
downsampled=down_sample(data, 5, 1)
downsampled_p1=prepare_data1(downsampled)
downsampled_p3=prepare_data3(downsampled_p1)
X=downsampled_p3.copy()
pca.fit_transform(X.astype(np.float64))

choice=np.random.choice(X.shape[0],int(X.shape[0]*0.2),replace=False)
train=X.loc[~np.isin(np.arange(X.shape[0]),choice)]
X_train=train.drop("is_attributed", axis=1)
y_train=train["is_attributed"]

test=X.iloc[choice]
X_test=test.drop("is_attributed", axis=1)
y_test=test["is_attributed"]

train_p3=prepare_data3(train)
X_train_p3=prepare_data3(X_train)
y_train_p3=y_train
X_test_p3=prepare_data3(X_test)
y_test_p3=y_test

clf=classifier.__class__(**classifier.get_params()).fit(X_train_p3,y_train_p3)
prediction = clf.predict(X_test_p3)
accuracy = accuracy_score(y_test_p3, prediction)
print("{} accuracy".format(classifier_name), accuracy)

#####################
classifier_name, classifier = "xg", xg_clf

pca = PCA(n_components=4)
downsampled=down_sample(data, 5, 1)
downsampled_p1=prepare_data1(downsampled)
downsampled_p3=prepare_data3(downsampled_p1)
X=downsampled_p3.copy()
pca.fit_transform(X.astype(np.float64))

choice=np.random.choice(X.shape[0],int(X.shape[0]*0.2),replace=False)
train=X.loc[~np.isin(np.arange(X.shape[0]),choice)]
X_train=train.drop("is_attributed", axis=1)
y_train=train["is_attributed"]

test=X.iloc[choice]
X_test=test.drop("is_attributed", axis=1)
y_test=test["is_attributed"]

train_p3=prepare_data3(train)
X_train_p3=prepare_data3(X_train)
y_train_p3=y_train
X_test_p3=prepare_data3(X_test)
y_test_p3=y_test

clf=classifier.__class__(**classifier.get_params()).fit(X_train_p3,y_train_p3)
prediction = clf.predict(X_test_p3)
accuracy = accuracy_score(y_test_p3, prediction)
print("{} accuracy".format(classifier_name), accuracy)
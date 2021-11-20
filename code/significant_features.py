import pandas as pd
import numpy as np
import os
import matplotlib.pylab as plt
from scipy import interp
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve,auc
import matplotlib.patches as patches
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn import svm

k = 6
celltype = 'cd8'
classif_levels = 'HDvsAll'
import string
filename = "~/Dropbox/MLGROUP/tcrcov/JP/kmers_" + celltype +"/su_" + celltype + "_"+ str(k) + "mers_cleaned.csv"
su = pd.read_csv(filename)

from sklearn.utils import resample
df_majority = su[su.d_cond!="HD"]
df_minority = su[su.d_cond=="HD"]
NEG_COUNT = len(su[su.d_cond!="HD"])

df_minority_upsampled = resample(df_minority, replace=True, n_samples=NEG_COUNT, random_state=430)

su = pd.concat([df_majority, df_minority_upsampled])
print(su.d_cond.value_counts())
print(su.shape)

su = su[su.d_cond.notnull()]
su["d_cond"].replace({"HD": 0, "mild": 1,"moderate":1, "severe":1}, inplace=True)
su.d_cond.value_counts()

new_index = np.arange(0,500,1)
su.index = new_index

labels = np.array(su['d_cond'])
features = su.drop('d_cond', axis = 1)
features = features.drop('Patient_ID', axis = 1)

X = features
y = labels

# Random Forest
### Methods: finding overlap between 1000 top kmers in 500 different random forest models and plotting their mean significance in order
model = RandomForestClassifier()
cv = RepeatedStratifiedKFold(n_splits=5,n_repeats=100,random_state=random_state)
folds = [(train,test) for train, test in cv.split(features, labels)]

fivehundredfold = []
for train, test in tqdm(folds, total=len(folds)):
    model.fit(features.loc[train,:],labels[train])
    rf_top_1000 = pd.Series(model.feature_importances_, index=X.columns).nlargest(10000)
    df = pd.DataFrame(rf_top_1000)
    df.reset_index(inplace=True)
    sig_kmers = list(df['index'])
    fivehundredfold.append(sig_kmers)

import functools
RF_SIG_FEATURES = list(functools.reduce(set.intersection,map(set,fivehundredfold)))

import string
filename = "/Users/violakyoungalee/Dropbox/MLGROUP/tcrcov/Code/VennDiagrams/" + celltype.upper() + "_"+ str(k) +"mer_" + classif_levels + "/"+"RF_SIG_FEATURES.txt"

textfile = open(filename, "w")
for element in RF_SIG_FEATURES:
    textfile.write(element + "\n")
textfile.close()

# Naive Bayes

from sklearn.inspection import permutation_importance
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=20)

model = GaussianNB()
model.fit(X_train,y_train)
imps = permutation_importance(model, X_test, y_test)
print(imps)

fivehundredfold = []
for number in range(500):
    model = GaussianNB()
    model.fit(X,y)
    nb_top_1000 = pd.Series(model.feature_importances_, index=X.columns).nlargest(1000)
    df = pd.DataFrame(nb_top_1000)
    df.reset_index(inplace=True)
    sig_kmers = list(df['index'])
    fivehundredfold.append(sig_kmers)
    
NB_SIG_FEATURES = list(functools.reduce(set.intersection,map(set,fivehundredfold)))
NB_SIG_FEATURES

textfile = open("NB_SIG_FEATURES.txt", "w")
for element in NB_SIG_FEATURES:
    textfile.write(element + "\n")
textfile.close()

import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from tqdm.notebook import tqdm
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.metrics import auc
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import BernoulliNB
from sklearn import svm
from sklearn.utils import resample
import plotly
import plotly.graph_objects as go

k = 6
celltype = 'cd8'
import string
filename = PATH_TO_FILES + celltype +"/su_" + celltype + "_"+ str(k) + "mers_cleaned.csv"
su = pd.read_csv(filename)

# subset for severity
HDvsCOVID = ['HD','moderate'] # modify this line for desired levels of classification
su = su[su.d_cond.isin(HDvsCOVID)]
su["d_cond"].replace({"HD": 0, "moderate":1}, inplace=True) # modify this line for desired levels of classification
su.d_cond.value_counts()

# split train and test sets
train, test = train_test_split(su, test_size=0.2)

# resample to balance the dataset AFTER train-test split
train_majority = train[train.d_cond!=0]
train_minority = train[train.d_cond==0]
NEG_COUNT = len(train[train.d_cond!=0])
train_minority_upsampled = resample(train_minority, replace=True, n_samples=NEG_COUNT, random_state=1)
train_balanced = pd.concat([train_majority, train_minority_upsampled])
print(train_balanced.d_cond.value_counts())
print(train_balanced.shape)

test_majority = test[test.d_cond!=0]
test_minority = test[test.d_cond==0]
NEG_COUNT_test = len(test[test.d_cond!=0])
test_minority_upsampled = resample(test_minority, replace=True, n_samples=NEG_COUNT_test, random_state=1)
test_balanced = pd.concat([test_majority, test_minority_upsampled])
print(test_balanced.d_cond.value_counts())
print(test_balanced.shape)

y_train = np.array(train_balanced['d_cond']) # y_train
train_features = train_balanced.drop('d_cond', axis = 1)
train_features = train_features.drop('Patient_ID', axis = 1)
X_train = np.array(train_features) # X_train

y_test = np.array(test_balanced['d_cond']) # y_test
test_features = test_balanced.drop('d_cond', axis = 1)
test_features = test_features.drop('Patient_ID', axis = 1)
X_test = np.array(test_features) # X_test

# Split the train set again into a train set and a validation set.
from sklearn.model_selection import RepeatedStratifiedKFold
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=100, random_state=random_state)
folds = [(train,test) for train, test in cv.split(X_train, y_train)]


metrics = ['auc', 'fpr', 'tpr', 'thresholds','classif_type']
results = {
    'train': {m:[] for m in metrics},
    'val'  : {m:[] for m in metrics},
    'test' : {m:[] for m in metrics}
}

results2 = {
    'train': {m:[] for m in metrics},
    'val'  : {m:[] for m in metrics},
    'test' : {m:[] for m in metrics}
}

results3 = {
    'train': {m:[] for m in metrics},
    'val'  : {m:[] for m in metrics},
    'test' : {m:[] for m in metrics}
}

results4 = {
    'train': {m:[] for m in metrics},
    'val'  : {m:[] for m in metrics},
    'test' : {m:[] for m in metrics}
}

results5 = {
    'train': {m:[] for m in metrics},
    'val'  : {m:[] for m in metrics},
    'test' : {m:[] for m in metrics}
}

# ML model training
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import BernoulliNB
from sklearn.linear_model import LogisticRegression
from sklearn import svm

#random_state=1

# all models are evaluated for final performance on the test set that was intially split before upsampling.

testset = X_test
testlabels = y_test

knn = KNeighborsClassifier(n_neighbors=3)
rf = RandomForestClassifier(n_estimators=100,random_state=random_state)
gb = GradientBoostingClassifier()
bnb = BernoulliNB(alpha = 1.0, binarize = 0.0, fit_prior = True, class_prior = None)
svc = svm.SVC(kernel='poly',C=20, degree=5,probability=True,random_state=random_state)
cv = RepeatedStratifiedKFold(n_splits=2,n_repeats=100,random_state=random_state)

for train, test in tqdm(folds, total=len(folds)):
#     dtrain = xgb.DMatrix(X_train[train,:], label=y_train[train])
#     dval   = xgb.DMatrix(X_train[test,:], label=y_train[test])
    trainset = X_train[train,:]
    valset = X_train[test,:]
    trainlabels = y_train[train]
    vallabels = y_train[test]

    # Random Forest Loop
    model  = rf.fit(trainset,y_train[train])
    sets = [trainset, valset, testset]
    labelsets = [trainlabels, vallabels, testlabels]
    for i,ds in enumerate(results.keys()):
        y_preds              = model.predict_proba(sets[i])
        labels               = labelsets[i]
        fpr, tpr, thresholds = roc_curve(labels, y_preds[:,1])
        results[ds]['fpr'].append(fpr)
        results[ds]['tpr'].append(tpr)
        results[ds]['thresholds'].append(thresholds)
        results[ds]['auc'].append(roc_auc_score(labels, y_preds[:,1]))
        results[ds]['classif_type'].append("Random Forest")

    # KNN loop
    model  = knn.fit(trainset,y_train[train])
    sets = [trainset, valset, testset]
    labelsets = [trainlabels, vallabels, testlabels]
    for i,ds in enumerate(results2.keys()):
        y_preds              = model.predict_proba(sets[i])
        labels               = labelsets[i]
        fpr, tpr, thresholds = roc_curve(labels, y_preds[:,1])
        results2[ds]['fpr'].append(fpr)
        results2[ds]['tpr'].append(tpr)
        results2[ds]['thresholds'].append(thresholds)
        results2[ds]['auc'].append(roc_auc_score(labels, y_preds[:,1]))
        results2[ds]['classif_type'].append("K Nearest Neighbors")

    # Gradient Boosting Loop
    model  = gb.fit(trainset,y_train[train])
    sets = [trainset, valset, testset]
    labelsets = [trainlabels, vallabels, testlabels]
    for i,ds in enumerate(results3.keys()):
        y_preds              = model.predict_proba(sets[i])
        labels               = labelsets[i]
        fpr, tpr, thresholds = roc_curve(labels, y_preds[:,1])
        results3[ds]['fpr'].append(fpr)
        results3[ds]['tpr'].append(tpr)
        results3[ds]['thresholds'].append(thresholds)
        results3[ds]['auc'].append(roc_auc_score(labels, y_preds[:,1]))
        results3[ds]['classif_type'].append("Gradient Boosting")

    # Bernoulli Naive Bayes Loop
    model  = bnb.fit(trainset,y_train[train])
    sets = [trainset, valset, testset]
    labelsets = [trainlabels, vallabels, testlabels]
    for i,ds in enumerate(results4.keys()):
        y_preds              = model.predict_proba(sets[i])
        labels               = labelsets[i]
        fpr, tpr, thresholds = roc_curve(labels, y_preds[:,1])
        results4[ds]['fpr'].append(fpr)
        results4[ds]['tpr'].append(tpr)
        results4[ds]['thresholds'].append(thresholds)
        results4[ds]['auc'].append(roc_auc_score(labels, y_preds[:,1]))
        results4[ds]['classif_type'].append("Naive Bayes")

    # SVM Loop
    model  = svc.fit(trainset,y_train[train])
    sets = [trainset, valset, testset]
    labelsets = [trainlabels, vallabels, testlabels]
    for i,ds in enumerate(results5.keys()):
        y_preds              = model.predict_proba(sets[i])
        labels               = labelsets[i]
        fpr, tpr, thresholds = roc_curve(labels, y_preds[:,1])
        results5[ds]['fpr'].append(fpr)
        results5[ds]['tpr'].append(tpr)
        results5[ds]['thresholds'].append(thresholds)
        results5[ds]['auc'].append(roc_auc_score(labels, y_preds[:,1]))
        results5[ds]['classif_type'].append("Support Vector Machine")

# save the intermediate data files into csv
df1 = pd.DataFrame(results)
df2 = pd.DataFrame(results2)
df3 = pd.DataFrame(results3)
df4 = pd.DataFrame(results4)
df5 = pd.DataFrame(results5)
df_total = pd.concat([df1, df2, df3, df4, df5])
file_name = celltype.upper() + "_"+ classif_levels + "_" + str(k) + "mer.csv"
df_total.to_csv(file_name)

# generate the final AUROC curve figures

import importlib
importlib.reload(go)

kind = 'test'
fpr_mean    = np.linspace(0, 1, 100)
interp_tprs = []
interp_tprs2 = []
interp_tprs3 = []
interp_tprs4 = []
interp_tprs5 = []

c_grid      = 'rgba(189, 195, 199, 0.5)'
# RF colors
c_fill      = 'rgba(52, 152, 219, 0.1)'
c_line      = 'rgba(52, 152, 219, 0.1)'
c_line_main = 'rgba(0, 0, 119, 1)'
# KNN colors
c_fill2      = 'rgba(255, 255, 25, 0.1)'
c_line2      = 'rgba(255, 255, 25, 0.1)'
c_line_main2 = 'rgba(200, 200, 0, 1.0)'
# Gradient Boosting colors
c_fill3      = 'rgba(144, 144, 200, 0.1)'
c_line3     = 'rgba(144, 144, 200, 0.1)'
c_line_main3 = 'rgba(104, 0, 163, 1)'
# Naive Bayes colors
c_fill4     = 'rgba(190, 40, 40, 0.1)'
c_line4     = 'rgba(190, 40, 40, 0.1)'
c_line_main4 = 'rgba(255, 0, 0, 1.0)'
# SVM colors
c_fill5     = 'rgba(90, 190, 90, 0.1)'
c_line5     = 'rgba(90, 190, 90, 0.1)'
c_line_main5 = 'rgba(0, 108, 0, 1)'

for i in range(100):
    fpr           = results[kind]['fpr'][i]
    tpr           = results[kind]['tpr'][i]
    interp_tpr    = np.interp(fpr_mean, fpr, tpr)
    interp_tpr[0] = 0.0
    interp_tprs.append(interp_tpr)

tpr_mean     = np.mean(interp_tprs, axis=0)
tpr_mean[-1] = 1.0
tpr_std      = 1*np.std(interp_tprs, axis=0)
tpr_upper    = np.clip(tpr_mean+tpr_std, 0, 1)
tpr_lower    = tpr_mean-tpr_std
auc          = np.mean(results[kind]['auc'])


for i in range(100):
    fpr2           = results2[kind]['fpr'][i]
    tpr2           = results2[kind]['tpr'][i]
    interp_tpr2    = np.interp(fpr_mean, fpr2, tpr2)
    interp_tpr2[0] = 0.0
    interp_tprs2.append(interp_tpr2)

tpr_mean2     = np.mean(interp_tprs2, axis=0)
tpr_mean2[-1] = 1.0
tpr_std2      = 1*np.std(interp_tprs2, axis=0)
tpr_upper2    = np.clip(tpr_mean2+tpr_std2, 0, 1)
tpr_lower2    = tpr_mean2-tpr_std2
auc2          = np.mean(results2[kind]['auc'])

for i in range(100):
    fpr3           = results3[kind]['fpr'][i]
    tpr3          = results3[kind]['tpr'][i]
    interp_tpr3   = np.interp(fpr_mean, fpr3, tpr3)
    interp_tpr3[0] = 0.0
    interp_tprs3.append(interp_tpr3)

tpr_mean3     = np.mean(interp_tprs3, axis=0)
tpr_mean3[-1] = 1.0
tpr_std3      = 1*np.std(interp_tprs3, axis=0)
tpr_upper3    = np.clip(tpr_mean3+tpr_std3, 0, 1)
tpr_lower3    = tpr_mean3-tpr_std3
auc3          = np.mean(results3[kind]['auc'])

for i in range(100):
    fpr4           = results4[kind]['fpr'][i]
    tpr4          = results4[kind]['tpr'][i]
    interp_tpr4   = np.interp(fpr_mean, fpr4, tpr4)
    interp_tpr4[0] = 0.0
    interp_tprs4.append(interp_tpr4)

tpr_mean4     = np.mean(interp_tprs4, axis=0)
tpr_mean4[-1] = 1.0
tpr_std4      = 1*np.std(interp_tprs4, axis=0)
tpr_upper4    = np.clip(tpr_mean4+tpr_std4, 0, 1)
tpr_lower4    = tpr_mean4-tpr_std4
auc4          = np.mean(results4[kind]['auc'])

for i in range(100):
    fpr5          = results5[kind]['fpr'][i]
    tpr5          = results5[kind]['tpr'][i]
    interp_tpr5   = np.interp(fpr_mean, fpr5, tpr5)
    interp_tpr5[0] = 0.0
    interp_tprs5.append(interp_tpr5)

tpr_mean5     = np.mean(interp_tprs5, axis=0)
tpr_mean5[-1] = 1.0
tpr_std5     = 1*np.std(interp_tprs5, axis=0)
tpr_upper5    = np.clip(tpr_mean5+tpr_std5, 0, 1)
tpr_lower5    = tpr_mean5-tpr_std5
auc5          = np.mean(results5[kind]['auc'])

fig = go.Figure()

# adding RF
fig.add_trace(go.Scatter(
        x          = fpr_mean,
        y          = tpr_upper,
        line       = dict(color=c_line, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'upper'))
fig.add_trace(     go.Scatter(
        x          = fpr_mean,
        y          = tpr_lower,
        fill       = 'tonexty',
        fillcolor  = c_fill,
        line       = dict(color=c_line, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'lower'),)

# adding KNN
fig.add_trace(go.Scatter(
        x          = fpr_mean,
        y          = tpr_upper2,
        line       = dict(color=c_line2, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'upper'))
fig.add_trace(     go.Scatter(
        x          = fpr_mean,
        y          = tpr_lower2,
        fill       = 'tonexty',
        fillcolor  = c_fill2,
        line       = dict(color=c_line2, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'lower'),)

# adding Gradient Boosting
fig.add_trace(go.Scatter(
        x          = fpr_mean,
        y          = tpr_upper3,
        line       = dict(color=c_line3, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'upper'))
fig.add_trace(     go.Scatter(
        x          = fpr_mean,
        y          = tpr_lower3,
        fill       = 'tonexty',
        fillcolor  = c_fill3,
        line       = dict(color=c_line3, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'lower'),)

# adding Naive Bayes
fig.add_trace(go.Scatter(
        x          = fpr_mean,
        y          = tpr_upper4,
        line       = dict(color=c_line4, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'upper'))
fig.add_trace(     go.Scatter(
        x          = fpr_mean,
        y          = tpr_lower4,
        fill       = 'tonexty',
        fillcolor  = c_fill4,
        line       = dict(color=c_line4, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'lower'),)

# adding SVM
fig.add_trace(go.Scatter(
        x          = fpr_mean,
        y          = tpr_upper5,
        line       = dict(color=c_line5, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'upper'))
fig.add_trace(     go.Scatter(
        x          = fpr_mean,
        y          = tpr_lower5,
        fill       = 'tonexty',
        fillcolor  = c_fill5,
        line       = dict(color=c_line5, width=1),
        hoverinfo  = "skip",
        showlegend = False,
        name       = 'lower'),)

# adding the models in reverse order of performance so that the legend shows up in order of performance
# things added later get added to the top
fig.add_trace(
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_mean2,
        line       = dict(color=c_line_main2, width=2),
        hoverinfo  = "skip",
        showlegend = True,
        name       = f'KNN AUROC: {auc2:.2f}'
))

fig.add_trace(
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_mean4,
        line       = dict(color=c_line_main4, width=2),
        hoverinfo  = "skip",
        showlegend = True,
        name       = f'BNB AUROC: {auc4:.2f}'
))

fig.add_trace(
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_mean,
        line       = dict(color=c_line_main, width=2),
        hoverinfo  = "skip",
        showlegend = True,
        name       = f'Random Forest AUROC: {auc:.2f}'
))

fig.add_trace(
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_mean5,
        line       = dict(color=c_line_main5, width=2),
        hoverinfo  = "skip",
        showlegend = True,
        name       = f'SVM AUROC: {auc5:.2f}'
))

fig.add_trace(
    go.Scatter(
        x          = fpr_mean,
        y          = tpr_mean3,
        line       = dict(color=c_line_main3, width=2),
        hoverinfo  = "skip",
        showlegend = True,
        name       = f'Gradient Boosting AUROC: {auc3:.2f}'
))


fig.add_shape(
    type ='line',
    line =dict(dash='dash'),
    x0=0, x1=1, y0=0, y1=1
)

fig.update_layout(
    title="CD8 HD vs. Moderate 6mer",
    template    = 'plotly_white',
    title_x     = 0.5,
    xaxis_title = "1 - Specificity",
    yaxis_title = "Sensitivity",
    width       = 800,
    height      = 800,
    legend      = dict(
        yanchor="bottom",
        xanchor="right",
        x=0.95,
        y=0.01,
    ),
    font=dict(
        size=23,
        color="Black"
    )
)

fig.update_yaxes(
    range       = [0, 1],
    gridcolor   = c_grid,
    scaleanchor = "x",
    scaleratio  = 1,
    linecolor   = 'black')

fig.update_xaxes(
    range       = [0, 1],
    gridcolor   = c_grid,
    constrain   = 'domain',
    linecolor   = 'black')

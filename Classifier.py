import sys
import numpy as np
import pandas as pd
import Feature_extraction
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score
train100 = np.array(pd.read_csv('./pssm-100.csv', header=None));
train38 = np.array(pd.read_csv('./pssm-38.csv', header=None));
train = np.hstack((train100, train38));
#binding function and other functions
def classifier1(test):
    target=list()
    for i in range(0,4278):
        target.append(0);
    for i in range(0,7306):
        target.append(1);
    target=np.array(target);
    clf = KNeighborsClassifier(n_neighbors=5);
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=None);
    predict=list();
    for trainIndex, testIndex in cv.split(train,target):
        trainSet=train[trainIndex];
        trainTarget=target[trainIndex];
        clf.fit(trainSet,trainTarget);
        predict.append(clf.predict(test)[0]);
    return predict;
#ATP binding function and other binding functions
def classifier2(test):
    trainATP=train[:4278];
    target=list()
    for i in range(0,2318):
        target.append(0);
    for i in range(2318,4278):
        target.append(1);
    target=np.array(target);
    clf = KNeighborsClassifier();
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=None);
    predict=list();
    for trainIndex, testIndex in cv.split(trainATP,target):
        trainSet=trainATP[trainIndex];
        trainTarget=target[trainIndex];
        clf.fit(trainSet, trainTarget);
        predict.append(clf.predict(test)[0]);
    return predict;
#heme binding,zinc ion binding,GTP binding,ADP binding
def classifier3(test):
    trainOther=train[2318:4278];
    target=list()
    # 2318+536=2854
    for i in range(2318, 2854):
        target.append(0)
    # #2854+562=3416
    for i in range(2854, 3416):
        target.append(1)
    # #3416+476=3892
    for i in range(3416, 3892):
        target.append(2)
    # #3892+386=4278
    for i in range(3892, 4278):
        target.append(3)
    target=np.array(target);
    clf = GradientBoostingClassifier(n_estimators=200);
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=None);
    predict=list();
    for trainIndex, testIndex in cv.split(trainOther,target):
        trainSet=trainOther[trainIndex];
        trainTarget=target[trainIndex];
        clf.fit(trainSet, trainTarget);
        predict.append(clf.predict(test)[0]);
    return predict;
if __name__ == '__main__':
    FileName=sys.argv[1];
    test=np.array([Feature_extraction.get138Fea(FileName)]);
    judge1=sum(classifier1(test));
    if judge1==0:
        judge2 = sum(classifier2(test));
        if judge2==0:
            print("ATP binding function");
        elif judge2==10:
            judge3 = sum(classifier3(test));
            if judge3==0:
                print("heme binding function");
            elif judge3==10:
                print("zinc ion binding function");
            elif judge3==20:
                print("GTP binding function");
            elif judge3==30:
                print("ADP binding function");
            else:
                print("four combined funtions");
        else:
            print("binding functions");
    else:
        print("other functions");
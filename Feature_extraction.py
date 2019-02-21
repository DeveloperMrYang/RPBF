import numpy as np
import math
def getPSSM(FileName):
    try:
        with open(FileName,"r") as f:
            l=f.readline();           #remove the first line.
            ls=f.readlines();
        tmp=ls[2:ls.index('\n')];
        pssm=list();
        for line in tmp:
            pssm.append([int(x) for x in line.split()[2:22]]);
        return np.array(pssm);
    except IOError:
        print("File not found");
        exit(-1);
    except ValueError:
        print("The format of file is error");
        exit(-2);
def get380Fea(res,mean):
    row=res.shape[0];
    if row<20:
        print("The length of a protein sequence should be greater than 20");
        exit(-3)
    fea400=np.zeros((20,20));
    for m in range(0,20):
        for n in range(0,20):
            D1=0;
            D2=0;
            for i in range(0,row):
                D1+=(res[i][m]-mean[m])**2;
                D2+=(res[i][n]-mean[n])**2;
            D1=math.sqrt(D1/row);
            D2=math.sqrt(D2/row);
            g=abs(m-n);
            for j in range(0,row-g):
                M=(res[j][m]-mean[m])*(res[j+g][n]-mean[n])/(D1*D2);
                fea400[m][n]+=M;
            fea400[m][n]/=(row-g);
    fea380=list();
    for m in range(0,20):
        for n in range(0,20):
            if m!=n:
                fea380.append(fea400[m][n]);
    return fea380;
def get20kFea(res):
    fea20k=np.zeros(100);
    row = res.shape[0];
    for k in range(1,6):
        iMax=row-k;
        for j in range(0,20):
            for i in range(0,iMax):
                fea20k[20*(k-1)+j]+=res[i][j]*res[i+k][j];
            fea20k[20 * (k - 1) + j]/=iMax;
    return fea20k;
def featureSelection(fea380):
    fea38=list();
    index=list([0,20,23,39,40,42,43,59,62,63,88,96,97,100,102,116,117,139,142,148,155,159,162,175,217,218,222,223,237,242,246,259,265,300,302,317,319,334]);
    for i in index:
        fea38.append(fea380[i]);
    return fea38;
def get138Fea(FileName):
    res=getPSSM(FileName);
    res=1/(1+math.e**(-res));
    mean=list();
    for i in range(0,20):
        mean.append(np.mean(res[:,i]));
    fea38=featureSelection(get380Fea(res,mean));
    fea100=get20kFea(res);
    return np.hstack((fea100,fea38));
#!/usr/bin/env python
import numpy as np
from cluster import Cluster
from celib import calc_cv,read_cluster_function,read_quantity
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split,ShuffleSplit,KFold
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model
import argparse

parser=argparse.ArgumentParser(description='Cluster model selection',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p',default='energy',dest='property',help="The property to expand")
parser.add_argument('-n',type=int,default=10,dest='times',help="The averaging times")
parser.add_argument('-t',type=str,default='',dest='title',help="The title of the plot")
args=parser.parse_args()

cluster_function=read_cluster_function()
quantity=np.loadtxt('all%s.out' %(args.property))
cluster=Cluster()

radius_cutoff=cluster.get_cluster_radius()

data=KFold(n_splits=10)
#data=ShuffleSplit(n_splits=args.times, train_size=120, test_size=30, random_state=None)
#X_train,X_test,y_train,y_test=train_test_split(cluster_function,quantity,train_size=120,test_size=40)
#models=cluster.construct_candidates(4,100)

model_size=[]
output=file('cv_test3_%s.dat' %(args.property),'w')

for i in radius_cutoff:
    model=cluster.get_cutoff_clusters(i)
    print i
    print model
    cv_list=np.array([])
    MAE_list=np.array([])
    test_list=np.array([])
    model_size.append(len(model))
    
    for train_index,test_index in data.split(cluster_function):
        if len(model)<len(train_index):
            X_train=cluster_function[train_index]
            X_test=cluster_function[test_index]
            y_train=quantity[train_index]
            y_test=quantity[test_index]

            eci=np.linalg.lstsq(X_train[:,model],y_train,rcond=None)[0]
            #ridge=linear_model.RidgeCV().fit(X_train[:,model],y_train)
            #eci=ridge.coef_
            #eci[0]=ridge.intercept_

            #u,s,vh=np.linalg.svd(X_train[:,model])

            cv=calc_cv(X_train[:,model],eci,y_train)
            MAE=mean_absolute_error(np.dot(X_train[:,model],eci),y_train)
            test_error=mean_absolute_error(np.dot(X_test[:,model],eci),y_test)

            #print 'Number of cluster:',len(eci)
            #print 'ECI:',eci
            #print 'Norm:',np.linalg.norm(eci)
            #print 'CV:',cv
            #print 'Singular Values:',s
            #print '\n'

            cv_list=np.append(cv_list,cv)
            MAE_list=np.append(MAE_list,MAE)
            test_list=np.append(test_list,test_error)

    e1=plt.errorbar(len(model),cv_list.mean(),yerr=cv_list.std(),fmt='s',color='C0',elinewidth=1,ecolor='C0',capsize=5,capthick=1)
    e2=plt.errorbar(len(model),MAE_list.mean(),yerr=MAE_list.std(),fmt='s',color='C1',elinewidth=1,ecolor='C1',capsize=5,capthick=1)
    e3=plt.errorbar(len(model),test_list.mean(),yerr=test_list.std(),fmt='s',color='C7',elinewidth=1,ecolor='C7',capsize=5,capthick=1)
    output.write('%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' %(len(model),cv_list.mean(),cv_list.std(),MAE_list.mean(),MAE_list.std(),test_list.mean(),test_list.std()))

output.close()

plt.xlabel('Number of Clusters')
plt.ylabel('CV/eV')
plt.legend((e1,e2,e3),('LOOCV','Training Score','Test Score'),loc='best',fontsize=12)
plt.title('%s-%s' %(args.title,args.property))
plt.show()

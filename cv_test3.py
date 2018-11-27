#!/usr/bin/env python
import numpy as np
from cluster import Cluster
from celib import calc_cv,read_cluster_function,read_quantity
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split,ShuffleSplit,KFold
from sklearn.metrics import mean_squared_error
from sklearn import linear_model
import argparse,math

parser=argparse.ArgumentParser(description='Cluster model selection. Null and point clusters are alway included. Scan two-body clusters first, the two-body cluster model with the smallest LOOCV is selected. Then three-body clustes are taken into consideration. Likewise, the model with the smallest LOOCV is selected. If it is larger than that of the two-body cluster model, then three-body clusters should not be added, else the similar procedure is performed for four-body clusters. K-fold CV is used to evaluate different models.\n Note: Before using this script, generate a large cluster pool for candidates by ATAT command \'corrdump\'.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p',default='energy',dest='property',help="The property to expand")
parser.add_argument('-k',type=int,default=10,dest='kfold',help="The k fold cross validation")
parser.add_argument('-t',type=str,default='',dest='title',help="The title of the plot")
args=parser.parse_args()

cluster_function=read_cluster_function()
quantity=np.loadtxt('all%s.out' %(args.property))
cluster=Cluster()

data=KFold(n_splits=args.kfold)

model_size=[]
#output=file('cv_test3_%s.dat' %(args.property),'w')
#output.write('#Number\tRadius\tLOOCV_min\tLOOCV_std\tTrain_mean\tTrain_std\tTest_mean\tTest_std\n')

initial_model=range(cluster.get_exclude_cluster_number(2)) #index of null and point clusters
print initial_model
for i in cluster.generate_candidates(2):
    model=initial_model+i
    #print 'Cut-off radius of clusters: %.3f' %(i) 
    print 'Index of clusters:\n',model
    print '\n'
    loocv_list=np.array([])
    train_list=np.array([])
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

            loocv=calc_cv(X_train[:,model],eci,y_train)
            train_error=math.sqrt(mean_squared_error(np.dot(X_train[:,model],eci),y_train))
            test_error=math.sqrt(mean_squared_error(np.dot(X_test[:,model],eci),y_test))

            #print 'Number of cluster:',len(eci)
            #print 'ECI:',eci
            #print 'Norm:',np.linalg.norm(eci)
            #print 'CV:',cv
            #print 'Singular Values:',s
            #print '\n'

            loocv_list=np.append(loocv_list,loocv)
            train_list=np.append(train_list,train_error)
            test_list=np.append(test_list,test_error)

    e1=plt.errorbar(len(model),loocv_list.mean(),yerr=loocv_list.std(),fmt='s',color='C0',elinewidth=1,ecolor='C0',capsize=5,capthick=1)
    e2=plt.errorbar(len(model),train_list.mean(),yerr=train_list.std(),fmt='s',color='C1',elinewidth=1,ecolor='C1',capsize=5,capthick=1)
    e3=plt.errorbar(len(model),test_list.mean(),yerr=test_list.std(),fmt='s',color='C7',elinewidth=1,ecolor='C7',capsize=5,capthick=1)
    #output.write('%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n' %(len(model),i,loocv_list.mean(),loocv_list.std(),train_list.mean(),train_list.std(),test_list.mean(),test_list.std()))

#output.close()

plt.xlabel('Number of Clusters')
plt.ylabel('RMSD/eV')
plt.legend((e1,e2,e3),('LOOCV','Training Score','Test Score'),loc='best',fontsize=12)
plt.title('%s-%s' %(args.title,args.property))
plt.show()

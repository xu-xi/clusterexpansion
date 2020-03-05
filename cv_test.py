#!/usr/bin/env python3
import numpy as np
from cluster import Cluster
from celib import calc_cv,read_cluster_function,read_quantity
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split,ShuffleSplit,RepeatedKFold
from sklearn.metrics import mean_squared_error
from sklearn import linear_model
import argparse,math

parser=argparse.ArgumentParser(description='Cluster model selection. Increase cluster radius gradually, and clusters (up to 4-body) within the cut-off radius are added into models. K-fold CV is used to evaluate different models.\n Note: Before using this script, generate a large cluster pool for candidates by ATAT command \'corrdump\'.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p',default='energy',dest='property',help="The property to expand")
parser.add_argument('-k',type=int,default=10,dest='kfold',help="Number of folds. Must be at least 2.")
parser.add_argument('-n',type=int,default=1,dest='repeated',help="Number of times cross-validator needs to be repeated.")
parser.add_argument('-t',type=str,default='',dest='title',help="The title of the plot")
args=parser.parse_args()

cluster_function=read_cluster_function()
quantity=np.loadtxt('all%s.out' %(args.property))
cluster=Cluster()

radius_cutoff=cluster.get_cluster_radius()

data=RepeatedKFold(n_splits=args.kfold,n_repeats=args.repeated)
#data=ShuffleSplit(n_splits=args.times, train_size=120, test_size=30, random_state=None)
#X_train,X_test,y_train,y_test=train_test_split(cluster_function,quantity,train_size=120,test_size=40)
#models=cluster.construct_candidates(4,100)

model_size=[]
output=open('cv_test_%s.dat' %(args.property),'w')
output.write('#Number\tRadius\tTrain_mean\tTrain_std\tTest_mean\tTest_std\n')

for i in radius_cutoff:
    model=cluster.get_cutoff_clusters(i)
    print('Cut-off radius of clusters: %.3f' %(i))
    print('Index of clusters:\n', model, '\n')
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

            train_error=math.sqrt(mean_squared_error(np.dot(X_train[:,model],eci),y_train))
            test_error=math.sqrt(mean_squared_error(np.dot(X_test[:,model],eci),y_test))

            #print 'Number of cluster:',len(eci)
            #print 'ECI:',eci
            #print 'Norm:',np.linalg.norm(eci)
            #print 'CV:',cv
            #print 'Singular Values:',s
            #print '\n'

            train_list=np.append(train_list,train_error)
            test_list=np.append(test_list,test_error)

    e2=plt.errorbar(i,train_list.mean(),yerr=train_list.std(),fmt='s',color='C0',elinewidth=1,ecolor='C0',capsize=5,capthick=1)
    e3=plt.errorbar(i,test_list.mean(),yerr=test_list.std(),fmt='s',color='C3',elinewidth=1,ecolor='C3',capsize=5,capthick=1)
    output.write('%i\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' %(len(model),i,train_list.mean(),train_list.std(),test_list.mean(),test_list.std()))

output.close()

#plt.xlabel('Number of Clusters')
plt.xlabel(r'Diameter of Clusters($\AA$)')
plt.ylabel('RMSD(eV)')
plt.legend((e2,e3),('Training Score','Test Score'),loc='best',fontsize=12)
#plt.title('%s-%s' %(args.title,args.property))
plt.show()

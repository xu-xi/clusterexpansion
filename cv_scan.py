#!/usr/bin/env python
import numpy as np
from cluster import Cluster
from celib import calc_cv,BIC
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split,ShuffleSplit,KFold
from sklearn.metrics import mean_squared_error
from sklearn import linear_model
import argparse,math,subprocess

parser=argparse.ArgumentParser(description='Model selection for cluster expansion. Null and point clusters are alway included. Scan two-body clusters first, the two-body cluster model with the smallest score is selected. Then three-body clustes are taken into consideration. Likewise, the model with the smallest score is selected. If it is larger than that of the two-body cluster model, then three-body clusters should not be added, else the similar procedure is performed for four-body clusters. K-fold score is used to evaluate different models.\n Note: Before using this script, generate a large cluster pool for candidates by ATAT command \'corrdump\'.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p',default='energy',dest='property',help="The property to expand")
parser.add_argument('-s',default='LOOCV',dest='score',help="Use BIC or LOOCV to assess different models")
args=parser.parse_args()


def optimal_model(score_type,cluster_order,initial_model,initial_score):
    'add n-order cluster to initial_model stepwise and return the optimal model with the smallest BIC or LOOCV score'
    model_size = []
    cutoff = []
    score_list = []
    best_score = initial_score
    
    subprocess.check_call('getclus > clusters.tmp',shell=True)
    radius = np.loadtxt('clusters.tmp')[:,1]

    print '\nAdd %s-body clusters:' %(cluster_order)
    for i in cluster.generate_candidates(cluster_order):
        model = initial_model + i
        print model
        model_size.append(len(model))
        cutoff.append(radius[model][-1])

        eci = np.linalg.lstsq(X[:,model],y,rcond=None)[0]
        if score_type == 'BIC':
            score = BIC(X[:,model],eci,y)
        elif score_type == 'LOOCV':
            score = calc_cv(X[:,model],eci,y)
        else:
            raise Exception("Only BIC and LOOCV are available right now.")

        score_list.append(score)

        datafile.write('%.2f\t%.5f\n' %(radius[model][-1],score))

        if score < best_score:
            best_score = score
            best_model = model

    datafile.write('\n')

    if 'best_model' in vars():
        return best_model,best_score,model_size,cutoff,score_list
    else:
        return initial_model,initial_score,model_size,cutoff,score_list


subprocess.check_call('clusterexpand -e %s' %(args.property),shell=True)

global X, y, cluster, datafile
X = np.loadtxt('allcorr.out')
y = np.loadtxt('all%s.out' %(args.property))
cluster = Cluster()

datafile = file('score.dat','w')
datafile.write('#Cluster_diameter\t%s\n' %(args.score))

#use a simple model as the initial model
simple_model = range(cluster.get_exclude_cluster_number(2)) #a simple model having only null and point clusters
eci = np.linalg.lstsq(X[:,simple_model],y,rcond=None)[0]
if args.score == "BIC":
    initial_score = BIC(X[:,simple_model],eci,y)
elif args.score == "LOOCV":
    initial_score = calc_cv(X[:,simple_model],eci,y)

best_model_2,best_score_2,model_size_2,cutoff_2,score_list_2 = optimal_model(args.score,2,simple_model,initial_score)
best_model_3,best_score_3,model_size_3,cutoff_3,score_list_3 = optimal_model(args.score,3,best_model_2,best_score_2)
best_model_4,best_score_4,model_size_4,cutoff_4,score_list_4 = optimal_model(args.score,4,best_model_3,best_score_3)

datafile.close()

plt.plot(cutoff_2,score_list_2,'+--',markersize=15,fillstyle='none',label='2-body')
plt.plot(cutoff_3,score_list_3,'^--',markersize=15,fillstyle='none',label='3-body')
plt.plot(cutoff_4,score_list_4,'s--',markersize=15,fillstyle='none',label='4-body')
plt.xlabel(r'Cluster Diameter/$\AA$')
plt.ylabel('%s Score/eV' %(args.score)) 
plt.legend(loc=0)
plt.show()

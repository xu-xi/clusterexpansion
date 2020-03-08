#!/usr/bin/env python
import math, argparse, subprocess
import numpy as np
import matplotlib.pyplot as plt
from cluster import Cluster
from celib import BIC, read_cluster_function

parser = argparse.ArgumentParser(description='Calculate and plot BIC for different cluster expansion models. Models are constructed hierarchically according to cluster radius.',
                formatter_class = argparse.ArgumentDefaultsHelpFormatter) 
#parser.add_argument('-p',default='energy',dest='property',help="The property to expand")
args = parser.parse_args()

#subprocess.check_call('clusterexpand -e energy',shell=True)

cluster = Cluster()
cluster_function = read_cluster_function()
str_num,clus_num = cluster_function.shape

y = np.loadtxt('allenergy.out')

models = cluster.construct_candidates(4,100) #cluster order & max cluster number

bic_set = []
parameter_num = []

#ECI = np.linalg.lstsq(cluster_function[:,:2],y,rcond=None)[0]
#min_bic = BIC(cluster_function[:,:2],ECI,y)
min_bic = 0

for model in models:
    ECI = np.linalg.lstsq(cluster_function[:,model],y,rcond=None)[0]
    parameter_num.append(len(ECI))
    bic = BIC(cluster_function[:,model],ECI,y)
    print(model, bic)
    bic_set.append(bic)
    if bic < min_bic:
        best_model = model
        min_bic = bic

#print(best_model)
print(parameter_num[bic_set.index(min(bic_set))], min_bic)

plt.plot(parameter_num,bic_set,'o--',ms=12,fillstyle='none')
plt.plot(parameter_num[bic_set.index(min(bic_set))],min_bic,'o',ms=12,color='red',label='The optimal model')
plt.legend(loc=0)
plt.xlabel('Number of Clusters')
plt.ylabel('BIC')
plt.show()

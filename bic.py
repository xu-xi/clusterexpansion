#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt
from cluster import Cluster
from celib import BIC,read_cluster_function

cluster=Cluster()

cluster_function=read_cluster_function()
str_num,clus_num=cluster_function.shape

energy=np.loadtxt('allenergy.out')
mbj=np.loadtxt('allmbj.out')

models=cluster.construct_candidates(4,100)

BIC_energy=[]
BIC_mbj=[]
parameter_num=[]

for model in models:
    ECI_energy=np.linalg.lstsq(cluster_function[:,model],energy,rcond=None)[0]
    ECI_mbj=np.linalg.lstsq(cluster_function[:,model],mbj,rcond=None)[0]
    parameter_num.append(len(ECI_energy))
    print model
    BIC_energy.append(BIC(cluster_function[:,model],ECI_energy,energy))
    BIC_mbj.append(BIC(cluster_function[:,model],ECI_mbj,mbj))
    #pd_energy.append(math.exp(-0.5*BIC(cluster_function[:,model],ECI_energy,energy))) #approximated posterior probability for CE model
    #pd_mbj.append(math.exp(-0.5*BIC(cluster_function[:,model],ECI_mbj,mbj)))

#print pd_mbj

#plt.plot(parameter_num,pd_energy/sum(pd_energy))
plt.plot(parameter_num,BIC_energy,'o--')
plt.xlabel('Number of Clusters')
plt.ylabel('BIC')
plt.show()
plt.plot(parameter_num,BIC_mbj,'o--')
plt.xlabel('Number of Clusters')
plt.ylabel('BIC')
plt.show()

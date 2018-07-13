#!/usr/bin/env python
import argparse
import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
from celib import read_cluster_function

parser=argparse.ArgumentParser(description='Bayesian Cluster Expansion',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
parser.add_argument('-p',default='energy',dest='property',help="The property to expand")
args=parser.parse_args()

X=read_cluster_function()
Y=np.loadtxt('all%s.out' %(args.property))

basic_model=pm.Model()

_,clus_num=X.shape
J=[0]*clus_num #ECI

with basic_model:
    for i in range(len(J)):
        J[i]=pm.Normal('J%i' %(i),mu=0,sd=10)

    #J[0]=pm.Normal('J0',mu=0,sd=10)
    #J[1]=pm.Normal('J1',mu=0,sd=10)

    sigma=pm.HalfNormal('sigma',sd=1)

    mu=0
    for i in range(len(J)):
        mu+=J[i]*X[:,i]
    #mu=np.dot(J,X.T) #expected value
    #mu=J[0]*X[:,0]+J[1]*X[:,1]

    Y_obs=pm.Normal('Y_obs',mu=mu,sd=sigma,observed=Y) #Likelihood

with basic_model:
    trace=pm.sample(500)

print pm.summary(trace).round(6)

pm.traceplot(trace,varnames=['J0'])
plt.show()
    
pm.plot_posterior(trace,varnames=['J0'])
plt.show()




#!/usr/bin/env python
import argparse
import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
from celib import read_cluster_function
plt.style.use('seaborn-darkgrid')

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
    J[0]=pm.Normal('J0',mu=0,sd=1)
    for i in range(len(J)-1):
        J[i+1]=pm.Normal('J%i' %(i+1),mu=0,sd=0.1)

    #J[0]=pm.Normal('J0',mu=0,sd=10)
    #J[1]=pm.Normal('J1',mu=0,sd=10)

    sigma=pm.HalfNormal('sigma',sd=1)

    mu=0
    for i in range(len(J)):
        mu+=J[i]*X[:,i]
    #mu=np.dot(J,X.T) #expected value
    #mu=J[0]*X[:,0]+J[1]*X[:,1]

    Y_obs=pm.Normal('Y_obs',mu=mu,sd=sigma,observed=Y) #Likelihood

#map_estimate=pm.find_MAP(model=basic_model)
#print(map_estimate)

with basic_model:
    trace=pm.sample(5000,chains=4,tune=1000,nuts_kwargs=dict(target_accept=0.9))

print(pm.summary(trace).round(6))

#pm.traceplot(trace)
#plt.show()


for i in range(len(J)):
    pm.traceplot(trace,varnames=['J%i' %(i)],figsize=(16,9))
    plt.savefig('gnzo-%s-bce-mc-J%i' %(args.property,i))
    pm.plot_posterior(trace,varnames=['J%i' %(i)])
    plt.savefig('gnzo-%s-bce-pd-J%i' %(args.property,i))

pm.plot_posterior(trace,varnames=['sigma'])
plt.savefig('gnzo-%s-bce-mc-sigma' %(args.property))

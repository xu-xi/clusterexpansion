#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy,subprocess,argparse

parser=argparse.ArgumentParser(description='plot the dependence of band gaps on temperature basing on the result of MC simulation',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
parser.add_argument('-f',default='mc.out',dest='mc_file',help="the Monte Carlo data file ")
parser.add_argument('-e',default='bandgap.ecimult',dest='eci_file',help="The ECIs for cluster expansion of bandgaps")
parser.add_argument('-x',type=float,required=True,dest='concentration',help="The concentration of component atom")
parser.add_argument('-t',type=str,dest='title',default='',help="The title of the plot")

args=parser.parse_args()


mcdata=numpy.loadtxt(args.mc_file)
eci=numpy.loadtxt(args.eci_file)
subprocess.call('getclus > clusters.tmp',shell=True)
multi=list(numpy.loadtxt('clusters.tmp')[:,0]) #multi of clusters

rd_bg=0
for i in range(len(eci)):
    rd_bg+=eci[i]*(2*args.concentration-1)**multi[i]

noc=len(file('bandgap.ecimult').readlines()) #number of clusters
mc_temps=list(mcdata[:,0]) #temperatures of MC simulation
clus_corr_funcs=mcdata[:,-noc:] #cluster correlation functions
bandgap=numpy.dot(clus_corr_funcs,eci)

plt.title(args.title)
plt.xlabel('Temperature/K')
plt.ylabel('Band Gap/eV')
plt.xlim(min(mc_temps),max(mc_temps))
plt.scatter(mc_temps,bandgap,c='.25',s=50,edgecolors='none',label='avg.')
plt.axhline(y=rd_bg,c='k',linestyle='dashed',label='random')
plt.legend()
plt.tight_layout()
plt.savefig('bg-mc.png')

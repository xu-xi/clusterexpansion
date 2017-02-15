#!/usr/bin/env python
import numpy,argparse
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description='draw Short-Range Order parameters from MC simulation and plot it',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
parser.add_argument('-f',default='mc.out',dest='mcfile',help="the Monte Carlo data file ")
parser.add_argument('-e',type=int,required=True,dest='exclude',help="The number of cluster are excluded,from null cluster,point cluster to many-body clusters")
parser.add_argument('-i',type=int,required=True,dest='include',help="The number of cluster are included")
parser.add_argument('-x',type=float,required=True,dest='concentration',help="The concentration of component atom")
parser.add_argument('-n',type=int,dest='order',default=2,help="The order of clusters")

args=parser.parse_args()

mcdata=numpy.loadtxt(args.mcfile)
noc=len(file('eci.out').readlines())-args.exclude #noc:number of clusters, except emtpy cluster and point cluster
#print 'The number of clusters is %s\n' % (noc) 

mc_temps=list(mcdata[:,0])
clus_corr_funcs=[]
for i in range(args.include):
    clus_corr_funcs.append(list(mcdata[:,-noc+i]))

random_corr=(2*args.concentration-1)**args.order
#print mc_temperatures
#clus_corr_funcs=(numpy.array(clus_corr_funcs)+1)/2
#print clus_corr_funcs
plt.xlabel('Temperature/K')
plt.ylabel('Cluster Correlation Functions')
plt.xlim(min(mc_temps),max(mc_temps))
plt.plot(mc_temps,[random_corr]*len(mc_temps),linestyle='dashed',c='k',label='random')
'''
plt.plot(mc_temps,clus_corr_funcs[0],c='r',label='2,1')
plt.plot(mc_temps,clus_corr_funcs[1],c='b',label='2,2')
plt.plot(mc_temps,clus_corr_funcs[2],c='g',label='2,3')
plt.plot(mc_temps,clus_corr_funcs[3],c='c',label='2,4')
plt.plot(mc_temps,clus_corr_funcs[4],c='m',label='2,5')
plt.plot(mc_temps,clus_corr_funcs[5],c='y',label='2,6')
'''
#plt.legend(bbox_to_anchor=(1.05,1),loc=2,borderaxespad=0.)

for i in range(len(clus_corr_funcs)):
    plt.plot(mc_temps,clus_corr_funcs[i],label='%s,%s' %(args.order,i+1) )
plt.legend()
plt.savefig('ccf-%s-%s.png' % (args.exclude,args.include))

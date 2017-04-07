#!/usr/bin/env python
import numpy,argparse,subprocess,os,sys
import matplotlib.pyplot as plt
from myfunc import read_clusters

def Main(Arglist):
    parser=argparse.ArgumentParser(description='plot Short-Range Order parameters from MC simulation',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-f',default='mc.out',dest='mcfile',help="the Monte Carlo data file ")
    parser.add_argument('-x',type=float,required=True,dest='concentration',help="The concentration of component atom")
    parser.add_argument('-o',type=int,dest='order',default=2,help="The order of clusters to plot")
    parser.add_argument('-t',type=str,dest='title',default='',help="The title of the plot")
    parser.add_argument('-T',type=int,dest='temp0',help="The initial temperature of the plot")
    parser.add_argument('--ft',type=str,dest='filetype',default='png',help="any filetype supported by Matplotlib")
    parser.add_argument('-n',type=int,dest='clus_number',help="The number of cluster to plot")
    args=parser.parse_args()

    #read data file
    mcdata=numpy.loadtxt(args.mcfile)
    mc_temps=list(mcdata[:,0])


    #plot cluster correlation functions from MC simulations
    clus_include,clus_exclude=read_clusters(args.order)
    if args.clus_number!=None and args.clus_number <= clus_include:
        clus_include=args.clus_number

    noc=len(file('eci.out').readlines())-clus_exclude

    clus_corr_funcs=[]
    for i in range(clus_include):
        clus_corr_funcs.append(list(mcdata[:,-noc+i]))
    for i in range(len(clus_corr_funcs)):
        plt.plot(mc_temps,clus_corr_funcs[i],label='%s,%s' %(args.order,i+1))

    #plot the short-range orders of random state
    random_corr=(2*args.concentration-1)**args.order
    plt.plot(mc_temps,[random_corr]*len(mc_temps),linestyle='dashed',linewidth=1,c='k',label='random')

    #details of the plot
    plt.title(args.title)
    plt.xlabel('Temperature/K')
    plt.ylabel('Cluster Correlation Functions')
    if args.temp0==None:
        plt.xlim(min(mc_temps),max(mc_temps))
    else:
        plt.xlim(args.temp0,max(mc_temps))
    plt.legend(loc=1,ncol=2,fontsize=10)
    #plt.tight_layout()
    plt.savefig('ccf-%s.%s' % (args.order,args.filetype))

if __name__ =="__main__":
    Main(sys.argv)

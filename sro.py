#!/usr/bin/env python
import numpy,argparse,subprocess,os,sys
import matplotlib.pyplot as plt

def read_clusters(order):
    subprocess.check_call('getclus > clusters.tmp',shell=True)
    clus_data=list(numpy.loadtxt('clusters.tmp')[:,0])
    os.remove('clusters.tmp')
    exclude_clus_num=0
    for i in range(order):
        exclude_clus_num+=clus_data.count(i)
    return clus_data.count(order),exclude_clus_num

def Main(Arglist):
    parser=argparse.ArgumentParser(description='plot Short-Range Order parameters from MC simulation',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-f',default='mc.out',dest='mcfile',help="the Monte Carlo data file ")
    #parser.add_argument('-i',type=int,dest='include',help="The number of cluster to plot")
    parser.add_argument('-x',type=float,required=True,dest='concentration',help="The concentration of component atom")
    parser.add_argument('-n',type=int,dest='order',default=2,help="The order of clusters to plot")
    parser.add_argument('-t',type=str,dest='title',default='',help="The title of the plot")
    args=parser.parse_args()

    clus_include,clus_exclude=read_clusters(args.order)
    mcdata=numpy.loadtxt(args.mcfile)
    noc=len(file('eci.out').readlines())-clus_exclude

    mc_temps=list(mcdata[:,0])
    clus_corr_funcs=[]
    for i in range(clus_include):
        clus_corr_funcs.append(list(mcdata[:,-noc+i]))

    random_corr=(2*args.concentration-1)**args.order
    plt.title(args.title)
    plt.xlabel('Temperature/K')
    plt.ylabel('Cluster Correlation Functions')
    plt.xlim(min(mc_temps),max(mc_temps))
    plt.plot(mc_temps,[random_corr]*len(mc_temps),linestyle='dashed',linewidth=1,c='k',label='random')

    for i in range(len(clus_corr_funcs)):
        plt.plot(mc_temps,clus_corr_funcs[i],label='%s,%s' %(args.order,i+1))
    plt.legend()
    plt.tight_layout()
    plt.savefig('ccf-%s.png' % (args.order))

if __name__ =="__main__":
    Main(sys.argv)

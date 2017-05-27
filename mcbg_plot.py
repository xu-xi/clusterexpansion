#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy,subprocess,argparse,sys

def Main(Arglist):
    parser=argparse.ArgumentParser(description='plot the dependence of band gaps on temperature basing on the result of MC simulation',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-f',default='mc.out',dest='mc_file',help="the Monte Carlo data file ")
    parser.add_argument('-e',default='bandgap.ecimult',dest='eci_file',help="The ECIs for cluster expansion of bandgaps")
    parser.add_argument('-x',type=float,required=True,dest='concentration',help="The concentration of component atom")
    parser.add_argument('-T',type=int,dest='temp0',help="The initial temperature of the plot")
    parser.add_argument('-t',type=str,dest='title',default='',help="The title of the plot")
    parser.add_argument('--ft',type=str,dest='filetype',default='png',help="The output filetype of the plot by Matplotlib")
    args=parser.parse_args()

    #read data file
    try:
        mcdata=numpy.loadtxt(args.mc_file)
        eci=numpy.loadtxt(args.eci_file)
    except IOError:
        print 'ERROR: can not read MC file or ECI file.';sys.exit(1)

    subprocess.call('getclus > clusters.tmp',shell=True)
    order=list(numpy.loadtxt('clusters.tmp')[:,0]) #order of clusters

    rd_bg=0
    for i in range(len(eci)):
        rd_bg+=eci[i]*(2*args.concentration-1)**order[i]

    noc=len(file(args.eci_file).readlines()) #number of clusters
    mc_temps=list(mcdata[:,0]) #temperatures of MC simulation
    clus_corr_funcs=mcdata[:,-noc:] #cluster correlation functions
    bandgap=numpy.dot(clus_corr_funcs,eci)

    plt.scatter(mc_temps,bandgap,c='.25',s=50,edgecolors='none',label='avg.')
    plt.axhline(y=rd_bg,c='k',linestyle='dashed',label='random')

    #details of the plot
    plt.title(args.title)
    plt.xlabel('Temperature/K')
    plt.ylabel('Band Gap/eV')
    if args.temp0==None:
        plt.xlim(min(mc_temps),max(mc_temps))
    else:
        plt.xlim(args.temp0,max(mc_temps))
    plt.legend()
    plt.tight_layout()
    plt.savefig('bg-mc.%s' % (args.filetype))

if __name__ =="__main__":
    Main(sys.argv)

#!/usr/bin/env python
import numpy,argparse,subprocess,os,sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from celib import read_clusters

#plt.style.use('ggplot')


def Main(Arglist):
    parser=argparse.ArgumentParser(description='plot Short-Range Order parameters from MC simulation',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-f',default='mc.out',dest='mcfile',help="The Monte Carlo data file ")
    parser.add_argument('-x',type=float,required=True,dest='concentration',help="The concentration of component atom")
    parser.add_argument('-o',type=int,dest='order',default=2,help="The order of clusters to plot")
    parser.add_argument('-t',type=str,dest='title',default='',help="The title of the plot")
    parser.add_argument('-T',type=int,dest='temp0',help="The initial temperature of the plot")
    #parser.add_argument('--ft',type=str,dest='filetype',default='png',help="Any filetype supported by Matplotlib")
    parser.add_argument('-n',type=int,dest='clus_number',help="The number of cluster to plot")
    parser.add_argument('-E',action='store_true',dest='plot_averaged_energy',help="Plot averaged energy vs temperature as well")
    args=parser.parse_args()

    #read data file
    k_b=8.617e-5  # reduced Boltzman's constant 
    mcdata=numpy.loadtxt(args.mcfile)
    mc_temps=list(mcdata[:,0])
    if args.plot_averaged_energy:
        if os.path.isfile('mcheader.out'):
            line_number=int(subprocess.check_output('grep Egc_mc mcheader.out',shell=True).split(':')[0]) 
            E_avg=list(mcdata[:,line_number-1])
            F=[]
            for i in range(len(E_avg)):
                Free_energy=E_avg[i]+mc_temps[i]*k_b*(args.concentration*numpy.log(args.concentration)+(1-args.concentration)*numpy.log(1-args.concentration))
                F.append(Free_energy)
            plt.plot(mc_temps,E_avg,'o--',c='C3')
            plt.plot(mc_temps,F,'^--',c='C0')
            plt.show()
            plt.close()
        else:
            E_avg=list(mcdata[:,2])
            E2=list(mcdata[:,5]*10**7) #variance of the energy
            #F=list(mcdata[:,4])

            fig,ax1=plt.subplots()
            ax1.plot(mc_temps,E_avg,'o--',linewidth=0.5,color='C3')
            ax1.set_xlabel('Temperature/K')
            ax1.set_ylabel('Averaged Energy/eV',color='C3')
            ax1.tick_params('y',colors='C3')

            ax2=ax1.twinx()
            ax2.plot(mc_temps,E2,'^--',linewidth=0.5,color='C0')
            ax2.set_ylabel(r'Variance of Energy/($10^{-7}$eV$^2$)',color='C0')
            ax2.tick_params('y',colors='C0')
            if args.temp0==None:
                plt.xlim(min(mc_temps),max(mc_temps))
            else:
                plt.xlim(args.temp0,max(mc_temps))
            plt.show()
            #plt.xlim(min(mc_temps),max(mc_temps))
            #plt.tight_layout()
            #plt.savefig('averaged_energy.eps')
            plt.close()

    #plt.scatter(mc_temps,Egc)
    #plt.savefig('egc.png')
    #plt.close()

    #plot cluster correlation functions from MC simulations
    clus_include,clus_exclude=read_clusters(args.order)
    if clus_include > 10:
        print('WARNING: too many clusters to display!')
        #clus_include=10
    if args.clus_number!=None and args.clus_number < clus_include:
        clus_include=args.clus_number

    #noc=len(file('eci.out').readlines())-clus_exclude
    noc=int(subprocess.check_output('getclus | wc -l',shell=True))-clus_exclude

    clus_corr_funcs=[]
    for i in range(clus_include):
        clus_corr_funcs.append(list(mcdata[:,-noc+i]))
    for i in range(len(clus_corr_funcs)):
        if args.order==1:
            plt.plot(mc_temps,[(j+1)/2 for j in clus_corr_funcs[i]],label='%s' %(i+1))
        else:
            plt.plot(mc_temps,clus_corr_funcs[i],label='%s' %(i+1))

    #plot the short-range orders of random state
    if args.order==1:
        plt.plot(mc_temps,[args.concentration]*len(mc_temps),linestyle='dashed',c='k',label='random')
    else:
        random_corr=(2*args.concentration-1)**args.order
        plt.plot(mc_temps,[random_corr]*len(mc_temps),linestyle='dashed',c='k',label='random')

    #details of the plot
    plt.title(args.title)
    plt.xlabel('Temperature/K')
    if args.order==1:
        plt.ylabel('Concentration x')
    else:
        plt.ylabel('Cluster Correlation Functions')

    if args.temp0==None:
        plt.xlim(min(mc_temps),max(mc_temps))
    else:
        plt.xlim(args.temp0,max(mc_temps))
    plt.legend(loc=0,ncol=2,fontsize=10)
    #plt.tight_layout()
    #plt.savefig('ccf-%s.%s' % (args.order,args.filetype))
    plt.show()

if __name__ =="__main__":
    Main(sys.argv)

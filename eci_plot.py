#!/usr/bin/env python
import numpy,subprocess,argparse,sys
import matplotlib.pyplot as plt
from myfunc import read_cluster_number

def eci_at_certain_temperature(ecifile,temperature):
    'return ECIs at given temperature from teci.out'
    if ecifile=='teci.out':
        teci=list(numpy.loadtxt(ecifile,usecols=(0)))[1:]
    else:
        teci=list(numpy.loadtxt(ecifile,usecols=(0)))
    trange=numpy.loadtxt('Trange.in')
    interval=int(trange[0]/(trange[1]-1))
    cluster_number=int(subprocess.check_output('getclus | wc -l',shell=True))
    index=temperature/interval*cluster_number
    return teci[index:(index+cluster_number)]

def Main(ArgList):
    parser=argparse.ArgumentParser(description='Plot ECIs against the index of clusters, including two-body clusters and higher order clusters. Temperature-dependent ECIs plotting is also supported.',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-f',type=str,dest='ecifile',default='eci.out',help="The ECI file")
    parser.add_argument('-o',dest='order',type=int,default=2,help="plot ECI from which order of clusters.")
    parser.add_argument('-T',type=int,dest='temp',nargs='+',help="temperature(s)")
    parser.add_argument('--ft',type=str,dest='filetype',default='png',help="any filetype supported by Matplotlib")
    #parser.add_argument('-m',type=int,dest='mode',default=0,help="the display mode of ECI")
    args=parser.parse_args()

    cluster_number=read_cluster_number()
    clus_exclude=0
    for i in range(args.order):
        clus_exclude+=cluster_number[i]
    for i in range(args.order):
        cluster_number.pop(0)

    def new_cluster_number(alist):
        blist=[]
        for i in range(len(alist)):
            k=0
            for j in range(i):
                k+=alist[j]
            blist.append(k)
        return blist
            

    #if args.mode==0:
   # 
   #     f,ax=plt.subplots(len(args.temp),sharex=True,sharey=True)
   #     
   #     for i in range(len(args.temp)):
   #         eci=eci_at_certain_temperature(args.ecifile,args.temp[i])[clus_exclude:]
   #         ax[i].bar(numpy.arange(1,len(eci)+1),eci,color='.25')
   #         ax[i].axhline(y=0,linewidth=.5,c='k')
   #         #ax[i].legend(['%sK' %(args.temp[i])],loc=1,fancybox=None)
   #         ax[i].text(0.8,0.6,'%sK' %(args.temp[i]),transform=ax[i].transAxes)
   #     
   #     plt.setp([a.get_xticklabels() for a in f.axes[:-1]],visible=False)
   #     f.subplots_adjust(hspace=0.0)
    if args.temp == None:
        eci=list(numpy.loadtxt('%s' %(args.ecifile)))[clus_exclude:]
        plt.bar(numpy.arange(1,len(eci)+1),eci,color='.25')
        #plt.xticks(numpy.arange(1,len(eci)))
        #plt.xticks([0,1,2,10,13],['null','point','pair','trip','quad'])
    else:
        for i in range(len(args.temp)):
            eci=eci_at_certain_temperature(args.ecifile,args.temp[i])[clus_exclude:]
            plt.plot(numpy.arange(1,len(eci)+1),eci,'o--',linewidth=1,label='T=%sK' %(args.temp[i]))
            plt.legend(loc=0,fontsize=10)

    for i in new_cluster_number(cluster_number):
        plt.axvline(x=i+0.5,linestyle='dashed',linewidth=.5,c='k')
    plt.axhline(y=0,linewidth=.5,c='k')
    plt.xlim(0.5,len(eci)+0.5)
    plt.xticks(range(1,len(eci),len(eci)/5))
    plt.xlabel('Index of clusters')
    plt.ylabel('ECI/eV')
    plt.show()

if __name__=='__main__':
    Main(sys.argv)

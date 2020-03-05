#!/usr/bin/env python3
import os,subprocess,numpy,argparse,sys,subprocess,math
import matplotlib.pyplot as plt
from celib import read_clusters,data_io
from sklearn.metrics import mean_squared_error,mean_absolute_error 


def Main(ArgList):
    parser=argparse.ArgumentParser(description='Plot the results of cluster expansion for the given property and output the data file.',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-p',default='energy',dest='property',help="The property to plot")
    parser.add_argument('-n',dest='name',help="The name of the property to show on the axis")
    parser.add_argument('--pa',action='store_true',dest='average',help="The quantity is per atom already")
    #parser.add_argument('-t',dest='title',type=str,default='',help="add the title of the plot of ECI")
    args=parser.parse_args()

    if args.average:
        cv=float(subprocess.check_output('clusterexpand -e -pa -cv %s | tail -1' %(args.property),shell=True))
    else:
        cv=float(subprocess.check_output('clusterexpand -e -cv %s | tail -1' %(args.property),shell=True))

    #read data file if there is one
    #if not os.path.isfile('%s.dat' %(args.property)):
    data_io(args.property,args.average)
        
    data=numpy.loadtxt('%s.dat' % (args.property),usecols=(2,3))
    vasp_data=list(data[:,0])
    ce_data=list(data[:,1])
    
    RMSD=math.sqrt(mean_squared_error(ce_data,vasp_data))
    MAE=mean_absolute_error(ce_data,vasp_data)
    #MAE=0 #Mean Absolute Error
    #for i in range(len(ce_data)):
    #    MAE+=abs(ce_data[i]-vasp_data[i])
    #MAE/=len(ce_data)

    #plot calculated and fitted values
    plt.scatter(vasp_data,ce_data,c='.25',s=50)

    a=max(max(ce_data),max(vasp_data))
    b=min(min(ce_data),min(vasp_data))
    d=(a-b)*0.2

    x=numpy.linspace(b-d,a+d,100)
    plt.plot(x,x,linestyle='dashed',linewidth=.5,c='k')

    plt.xlim(b-d,a+d)
    plt.ylim(b-d,a+d)
    if args.name == None:
        plt.xlabel('Calculated %s/eV' % (args.property))
        plt.ylabel('Fitted %s/eV' % (args.property))
    else:
        plt.xlabel('Calculated %s/eV' % (args.name))
        plt.ylabel('Fitted %s/eV' % (args.name))

    plt.text(a-2*d,b,'CV = %.4f eV\nRMSD = %.4f eV' %(cv,RMSD))
    #plt.savefig('%s-ce.%s' %(args.property,args.filetype))
    plt.show()
    
if __name__=='__main__':
    Main(sys.argv)

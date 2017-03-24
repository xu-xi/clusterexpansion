#!/usr/bin/env python
import os,subprocess,numpy,argparse,sys,subprocess
import matplotlib.pyplot as plt
from myfunc import read_clusters

def data_io(data,average):
    'collect successfully finished results and output them to a data file'
    lat_atom_number=len(file('lat.in').readlines())-6
    datafile=file('%s.dat' % (data),'w')
    datafile.write('#index\t%s-ce\t%s-vasp\n' % (data,data))
    for item in os.listdir(os.environ['PWD']):
        fullpath=os.path.join(os.environ['PWD'],item)
        if os.path.isdir(fullpath):
            os.chdir(fullpath)
            if os.path.isfile(data) and not os.path.isfile('error'):
                datafile.write(str(os.path.basename(fullpath))+'\t')
                sc_atom_number=len(file('str.out').readlines())-6
                vasp_data=float(subprocess.check_output('cat %s' %(data),shell=True))
                if average==True:
                    ce_data=float(subprocess.check_output('corrdump -c -mi -l=../lat.in -cf=../clusters.out -eci=../%s.ecimult' %(args.property),shell=True))
                else:
                    ce_data=(sc_atom_number/lat_atom_number)*float(subprocess.check_output('corrdump -c -l=../lat.in -cf=../clusters.out -eci=../%s.eci' %(data),shell=True))
                datafile.write('%.5f\t%.5f\n' %(ce_data,vasp_data))
                os.chdir('../')
    datafile.close()

def Main(ArgList):
    parser=argparse.ArgumentParser(description='plot results of cluster expansion for the given property.',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-p',default='energy',dest='property',help="the property to plot")
    parser.add_argument('-a',action='store_true',dest='average',help="decide whether the property to be averaged")
    parser.add_argument('--ft',dest='filetype',default='png',help="any filetype supported by matplotlib")
    parser.add_argument('-t',dest='title',type=str,default='',help="the title of the plot of ECI")
    parser.add_argument('-o',dest='order',type=int,default=1,help="plot ECI from which order of clusters. The default value is 1 since the ECI of null cluster is just a shift")
    args=parser.parse_args()

    #read data file
    data_io(args.property,args.average)
    data=numpy.loadtxt('%s.dat' % (args.property))
    index=list(data[:,0])
    ce_data=list(data[:,1])
    vasp_data=list(data[:,2])

    #plot calculated and fitted values
    plt.scatter(ce_data,vasp_data,c='.25',s=50)

    a=max(max(ce_data),max(vasp_data))
    b=min(min(ce_data),min(vasp_data))
    d=(a-b)*0.2

    x=numpy.linspace(b-d,a+d,100)
    plt.plot(x,x,linestyle='dashed',linewidth=.5,c='k')

    plt.xlim(b-d,a+d)
    plt.ylim(b-d,a+d)
    plt.xlabel('Fitted %s/eV' % (args.property))
    plt.ylabel('Calculated %s/eV' % (args.property))
    plt.savefig('%s-ce.%s' %(args.property,args.filetype))
    plt.close()
    
    #plot ECIs
    clus_include,clus_exclude=read_clusters(args.order)
    eci=list(numpy.loadtxt('%s.eci' %(args.property)))[clus_exclude:]
    plt.bar(numpy.arange(1,len(eci)+1),eci,color='.25')
    plt.axhline(y=0,linewidth=.5,c='k')

    plt.title(args.title)
    plt.xlabel('Index of clusters')
    plt.ylabel('ECI/eV')
    plt.xlim(0.5,len(eci)+0.5)
    #plt.xticks(numpy.arange(1,len(eci)))
    #plt.xticks([0,1,2,10,13],['null','point','pair','trip','quad'])
    plt.savefig('%s-eci.%s' %(args.property,args.filetype))
    plt.close()

if __name__=='__main__':
    Main(sys.argv)

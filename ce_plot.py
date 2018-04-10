#!/usr/bin/env python
import os,subprocess,numpy,argparse,sys,subprocess
import matplotlib.pyplot as plt
from celib import read_clusters

def data_io(data,average=False):
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
                    ce_data=float(subprocess.check_output('corrdump -c -mi -l=../lat.in -cf=../clusters.out -eci=../%s.ecimult' %(data),shell=True))
                else:
                    ce_data=(sc_atom_number/lat_atom_number)*float(subprocess.check_output('corrdump -c -l=../lat.in -cf=../clusters.out -eci=../%s.eci' %(data),shell=True))
                datafile.write('%.5f\t%.5f\n' %(ce_data,vasp_data))
            os.chdir('../')
    datafile.close()

def Main(ArgList):
    parser=argparse.ArgumentParser(description='Plot the results of cluster expansion for the given property and output the data file.',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-p',default='energy',dest='property',help="The property to plot")
    parser.add_argument('-n',dest='name',help="The name of the property to show on the axis")
    parser.add_argument('--pa',action='store_true',dest='average',help="The quantity is per atom already")
    parser.add_argument('--cv',dest='print_cv',action='store_true',help="Print Cross Validation and Mean Standard Error within the plot")
    #parser.add_argument('--ft',dest='filetype',default='png',help="any filetype supported by matplotlib")
    #parser.add_argument('-t',dest='title',type=str,default='',help="add the title of the plot of ECI")
    args=parser.parse_args()

    if args.average:
        cv=float(subprocess.check_output('clusterexpand -e -pa -cv %s | tail -1' %(args.property),shell=True))
    else:
        cv=float(subprocess.check_output('clusterexpand -e -cv %s | tail -1' %(args.property),shell=True))
    #read data file
    data_io(args.property,args.average)
    data=numpy.loadtxt('%s.dat' % (args.property),usecols=(1,2))
    ce_data=list(data[:,0])
    vasp_data=list(data[:,1])
    
    MAE=0 #Mean Absolute Error
    for i in range(len(ce_data)):
        MAE+=abs(ce_data[i]-vasp_data[i])
    MAE/=len(ce_data)

    #plot calculated and fitted values
    plt.scatter(ce_data,vasp_data,c='.25',s=50)

    a=max(max(ce_data),max(vasp_data))
    b=min(min(ce_data),min(vasp_data))
    d=(a-b)*0.2

    x=numpy.linspace(b-d,a+d,100)
    plt.plot(x,x,linestyle='dashed',linewidth=.5,c='k')

    plt.xlim(b-d,a+d)
    plt.ylim(b-d,a+d)
    if args.name == None:
        plt.xlabel('Fitted %s/eV' % (args.property))
        plt.ylabel('Calculated %s/eV' % (args.property))
    else:
        plt.xlabel('Fitted %s/eV' % (args.name))
        plt.ylabel('Calculated %s/eV' % (args.name))

    if args.print_cv:
        plt.text(a-2*d,b,'CV = %.3f eV\nMAE = %.3f eV' %(cv,MAE))
    else:
        print 'CV = %.4f eV' %(cv)
        print 'MAE = %.4f eV' %(MAE)
    #plt.savefig('%s-ce.%s' %(args.property,args.filetype))
    plt.show()
    
if __name__=='__main__':
    Main(sys.argv)

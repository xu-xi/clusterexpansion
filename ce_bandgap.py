#!/usr/bin/env python
import os,subprocess,numpy,argparse
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description='plot results of cluster expansion for given property',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
parser.add_argument('-p',default='bandgap',dest='property',help="the property to plot")
args=parser.parse_args()

lat_atom_number=len(file('lat.in').readlines())-6
datafile=file('%s.dat' % (args.property),'w')
datafile.write('#index %s-ce %s-vasp\n' % (args.property,args.property))
for item in os.listdir(os.environ['PWD']):
        fullpath=os.path.join(os.environ['PWD'],item)
        if os.path.isdir(fullpath):
                os.chdir(fullpath)
                if os.path.isfile(args.property) and not os.path.isfile('error'):
                        #datafile.write(str(os.path.basename(fullpath))+' ')
                        sc_atom_number=len(file('str.out').readlines())-6
                        vasp_data=float(subprocess.check_output('cat %s' %(args.property),shell=True))
                        #ce_energy=(sc_atom_number/lat_atom_number)*float(subprocess.check_output('corrdump -c -l=../lat.in -cf=../clusters.out -eci=../bandgap.eci',shell=True))
                        ce_data=float(subprocess.check_output('corrdump -c -mi -l=../lat.in -cf=../clusters.out -eci=../%s.ecimult' %(args.property),shell=True))
                        datafile.write('%.5f %.5f\n' %(ce_data,vasp_data))
                        os.chdir('../')
datafile.close()

#os.chdir('../')
data=numpy.loadtxt('%s.dat' % (args.property))
#index=list(data[:,0])
ce_data=list(data[:,0])
vasp_data=list(data[:,1])

plt.xlabel('Fitted Bandgap/eV')
plt.ylabel('Calculated Bandgap/eV')
plt.scatter(ce_data,vasp_data,c='.25',s=50)

a=max(max(ce_data),max(vasp_data))
b=min(min(ce_data),min(vasp_data))
d=(a-b)*0.2
plt.xlim(b-d,a+d)
plt.ylim(b-d,a+d)

x=numpy.linspace(b-d,a+d,100)
plt.plot(x,x,linestyle='dashed',linewidth=.5,c='k')

plt.savefig('%s-ce.png' %(args.property))
plt.close()

eci=list(numpy.loadtxt('%s.eci' %(args.property)))
plt.xlabel('Index of clusters')
#plt.ylabel(r'$J_{\alpha}m_{\alpha}$/eV')
plt.ylabel('ECI/eV')
plt.xlim(-0.5,len(eci)-0.5)
#plt.xticks(numpy.arange(len(eci)))
plt.xticks([2,10,13],['pair','trip','quad'])

plt.bar(numpy.arange(len(eci)),eci,color='.25',width=1)
plt.axhline(y=0,linewidth=.5,c='k')
plt.savefig('%s-eci.png' %(args.property))

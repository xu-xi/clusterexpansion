#!/usr/bin/env python
import os,subprocess,numpy
import matplotlib.pyplot as plt

pku_red=(0.545,0,0.07)

lat_atom_number=len(file('lat.in').readlines())-6
datafile=file('energy_data.out','w')
datafile.write('#index\tvasp_energy\tce_energy\n')
for item in os.listdir(os.environ['PWD']):
	fullpath=os.path.join(os.environ['PWD'],item)
	if os.path.isdir(fullpath):
		os.chdir(fullpath)
		if os.path.isfile('energy') and not os.path.isfile('error'):
			#datafile.write(str(os.path.basename(fullpath))+'\t')
			sc_atom_number=len(file('str.out').readlines())-6
			vasp_energy=float(subprocess.check_output('cat energy',shell=True))
			ce_energy=(sc_atom_number/lat_atom_number)*float(subprocess.check_output('corrdump -c -l=../lat.in -cf=../clusters.out -eci=../energy.eci',shell=True))
			datafile.write('%.5f\t%.5f\n' % (ce_energy,vasp_energy))
			os.chdir('../')
datafile.close()

data=numpy.loadtxt('energy_data.out')
#index=list(data[:,0])
ce_data=list(data[:,0])
vasp_data=list(data[:,1])

plt.xlabel('Fitted Energy/eV')
plt.ylabel('Calculated Energy/eV')
plt.scatter(ce_data,vasp_data,s=50,c='.25',edgecolors='none')
#plt.plot(x,x,linestyle='-',c='.25')
#plt.scatter(ce_data,vasp_data,c=pkr,s=50,edgecolors='none')


a=max(max(ce_data),max(vasp_data))
b=min(min(ce_data),min(vasp_data))
d=(a-b)*0.2
plt.xlim(b-d,a+d)
plt.ylim(b-d,a+d)

x=numpy.linspace(b-d,a+d,100)
plt.plot(x,x,linestyle='dashed',linewidth=.5,c='k')

plt.savefig('energy-ce.png')
plt.close()

#Null cluster and point cluster are excluded
eci=list(numpy.loadtxt('energy.eci'))[2:]
plt.xlabel('Index of clusters')
plt.ylabel('ECI/eV')
plt.xlim(0.5,len(eci)+0.5)
#plt.xticks([1,9,12],['pair','trip','quad'])
plt.bar(numpy.arange(1,len(eci)+1),eci,color='.25',width=1)
plt.axhline(y=0,linewidth=.5,c='k')

plt.savefig('energy-eci.png')


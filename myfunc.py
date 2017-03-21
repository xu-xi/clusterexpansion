#!/usr/bin/env python
import sys,re,math,commands,itertools,shutil,os,random,subprocess
from numpy import *
from fractions import Fraction
from compiler.ast import flatten

_author_='Xu Xi'

def listwrite(afile,alist):
        'write a two-dimension list to file' 
	for i in alist:
		for j in i:
			afile.write(str(j)+' ')
		afile.write('\n')

def projection(vector,bases):
	'return a list of the coefficients of a vector on a base'
	return list(dot(linalg.inv(transpose(bases)),vector))

def translation(site,platvect,slatvect):
	'generate a list of sites according to translation symmetry'
	sites=[]
	a_project=[projection(slatvect[x],platvect)[0] for x in range(3)]
	b_project=[projection(slatvect[x],platvect)[1] for x in range(3)]
	c_project=[projection(slatvect[x],platvect)[2] for x in range(3)]
        for i in range(int(min(a_project))-2,int(max(a_project))+2,1):
		for j in range(int(min(b_project))-2,int(max(b_project))+2,1):
			for k in range(int(min(c_project))-2,int(max(c_project))+2,1):
				tsite=array(site)+platvect[0]*i+platvect[1]*j+platvect[2]*k
				project=projection(tsite,slatvect)
				if min(project)>=0 and max(project)<1:
					tsite=list(tsite)
					tsite.append('')
					sites.append(tsite)
	return sites

def occupy(strfile,sublattice,atomlist):
	i=0
	for sites in sublattice:
		if len(sites)!=len(atomlist[i]):
			print 'ERROR: the number of site do not equal that of atoms'
			print sites
			print atomlist[i]
			sys.exit(1)
		j=0
		for site in sites:
			site[3]=atomlist[i][j]
			j+=1
		listwrite(strfile,sites)
		i+=1

def occupysaj(strfile,sublattice,atomlist):
	for i in range(len(sublattice)):
		for j in range(len(sublattice[i])):
			sublattice[i][j].append(atomlist[i][j])
	listwrite(strfile,sublattice)

def initstr(lattice,superlattice):
	'generate a random occupied structure'
	try:
		latfile=file(lattice,'r')
		supercell=file(superlattice,'r')
	except:
		raise IOError

	strfile=file('str.out','w')
	for i in range(3):
		strfile.write(latfile.readline())
	latvect=[]
	strvect=[]
	for i in range(3):
		latvect.append(map(float,latfile.readline().split()))
		strvect.append(map(float,supercell.readline().split()))
	listwrite(strfile,strvect)
	strvect=array(strvect)
	latvect=array(latvect)
	pvolume=abs(dot(latvect[0],cross(latvect[1],latvect[2])))
	svolume=abs(dot(strvect[0],cross(strvect[1],strvect[2])))
	n=int(svolume/pvolume)
	
	isublattices=[]
	iatomlist=[]
	vsublattices=[]
	sublattice_atoms=[]
	sublattice_c=[]
	sitecount=[]
	vatomlist=[]
	for line in latfile:
		g=re.findall(r'(\w+)=(\d+)/(\d+)',line.split()[3])
		site=map(float,line.split()[:3])
		sites=translation(site,latvect,strvect)
		if g==[]:
			isublattices.append(sites)
			atoms=[line.split()[3]]*n
			iatomlist.append(atoms)
		else:
			site_atom=[]
			site_c=[]
			for i in range(len(g)):
				site_atom.append(g[i][0])
				site_c.append(Fraction(int(g[i][1]),int(g[i][2])))
			if sum(site_c)!=1:
				print 'ERROR: the occupy probabilities of one site must sum up to 1\n'
				print line
				sys.exit(1)
			#atoms=[]
			#for i in range(len(a)):
			#		b=[a[i]]*int(n*c[i])
			#		atoms+=b
			#random.shuffle(atoms)
			#if a in vatoms:
			#	vatomlist[vatoms.index(a)]+=atoms
			#	vsublattices[vatoms.index(a)]+=sites
			#else:
			#	vatoms.append(a)
			#	vatomlist.append(atoms)
			#	vsublattices.append(sites)
			if site_atom in sublattice_atoms:
				sitecount[sublattice_atoms.index(site_atom)]+=1
				vsublattices[sublattice_atoms.index(site_atom)]+=sites
			else:
				sublattice_atoms.append(site_atom)
				sublattice_c.append(site_c)
				sitecount.append(1)
				vsublattices.append(sites)
	for i in range(len(sublattice_atoms)):
		atoms=[]
		for j in range(len(sublattice_atoms[i])):
			b=[sublattice_atoms[i][j]]*int(n*sublattice_c[i][j]*sitecount[i])
			atoms+=b
		random.shuffle(atoms)
		vatomlist.append(atoms)

	#print isublattices,iatomlist
	#print vsublattices,vatomlist
	occupy(strfile,isublattices,iatomlist)
	occupy(strfile,vsublattices,vatomlist)
	latfile.close()
	strfile.close()
	return isublattices,iatomlist,vsublattices,vatomlist,sublattice_c,n

def ce_energy(ecifile,strfile):
	'get energy from cluster expansion by the ATAT code corrdump'
	cmd='corrdump -c -eci=%s -s=%s' %(ecifile,strfile)
	failure,output=commands.getstatusoutput(cmd)
	if failure:
		print 'ERROR: running corrdump failed\n%s\n%s' %(cmd,output)
		sys.exit(1)
	return float(output)

def mdenergy(mcstep):
	'calculate Madelung energy' 
	#os.chdir(str(mcstep))
	commands.getstatusoutput('runstruct_vasp -nr&&py_conv -f POSCAR -i vasp -o w2k')
	os.rename('POSCAR.struct',str(mcstep)+'.struct')
	commands.getstatusoutput('py_initmad -f '+str(mcstep)+'.struct')
	cmd='calcmad '+str(mcstep)+'.inmad'+'|grep Energy|head -1|awk \'{print $3}\''
	failure,output=commands.getstatusoutput(cmd)
	os.chdir('../')
	return float(output)

def sort(dictionary,number):
	dictionary=sorted(dictionary.iteritems(),key=lambda d:d[0])
	ground_index=dictionary[0][1]
	ground_energy=dictionary[0][0]
	os.rename(str(ground_index),'0')
	print '#the energy of ground state:',ground_energy,'\n','#the number of inequal structures:',len(dictionary)
	if number!=0:
		dictionary=dictionary[:number]
	dictionary=dict(dictionary)
	dictionary[ground_energy]='0'
	return dictionary

def swap(alist):
	'swap the position of two randomly choosen elements from a list'
	import random
	while 1:
		l=len(alist)
		e1=random.randrange(l)
		e2=random.randrange(l)
		if alist[e1]!=alist[e2]:
			alist[e1],alist[e2]=alist[e2],alist[e1]
			return alist
			break

def enumer(alist):
	'generate all possible combinations of a list of two elements'
	atoms=list(set(alist))
	sites=itertools.combinations(range(len(alist)),alist.count(atoms[0]))
	#allist=[]
	for i in sites:
		alist=[atoms[1]]*len(alist)
		for j in i:
			alist[j]=atoms[0]
		yield alist

def enumer_test(alist):
	'generate all posible combinations of a list'
	elements=list(set(alist))
	element_num=[]
	length=len(alist)
	for i in elements:
		element_num.append(alist.count(elements[i]))
	for j in range(len(elements)-1):
		selected_sites=itertools.combinations(range(length),element_num[i])
	
	
def getbandgap():
	'return a dictionary whose keys are index and values are bandgaps'
	bandgap={}
	rootdir=os.environ['PWD']
        filelist=os.listdir(rootdir)
        for ifile in filelist:
		fullpath=os.path.join(rootdir,ifile)
		if os.path.isdir(fullpath):
			os.chdir(fullpath)
			if os.path.isfile('OUTCAR'):
				failure,output=commands.getstatusoutput('pv_gap.py OUTCAR|tail -1|awk \'{print $3}\'')
				if failure:
					print 'Error in '+str(fullpath)
				else:
					bandgap[str(ifile)]=float(output)
			else:
				print 'No OUTCAR in'+str(fullpath) 
	os.chdir(rootdir)
	return bandgap

def touch(filename,value):
	'write a file with a value for cluster expansion'
	with open(filename) as afile:
		afile.write('%s' %(value))

def objective_correlation_functions(concentration,cluster):
	'generate the objective cluster correlation functions or read it from input file'
	obj_corr=[]
	for i in range(len(cluster)):
		obj_corr.append((2*concentration[0][0]-1)**cluster[i])
	if os.path.isfile('corr_func.in'):
		read_corr=[]
		for line in file('corr_func.in'):
			read_corr.append(float(line.split()[0]))
		for j in range(len(read_corr)):
			obj_corr[j]=read_corr[j]
	return obj_corr

def Boltzmann_factor(energy,temperature):
	if temperature == -1: #-1 stands for the infinity temperature case
		return 1
	elif temperature == 0:
		print 'I cannot handle the zero temperature case.';sys.exit(1)
	else:
		exponent=-(energy*1.161E4)/temperature
		if exponent <=-50:
			return 0
		else:
			return math.pow(math.e,exponent)


def calc_correlation_functions(weight,obj_corr,radius):
	'calculate correlation functions'
	failure,output=commands.getstatusoutput('corrdump -c -l=../lat.in -s=str.out')
	if failure:
                print 'ERROR: running corrdump failed!'
                print output;sys.exit(1)
        else:
                str_corr=map(float,output.split())
	correlations=[]
	if len(str_corr)!=len(obj_corr):
		print 'ERROR: the number of clusters whose correlation functions to be calculated do not match';sys.exit(1)
	else:
		for i in range(len(obj_corr)):
			correlations.append(str_corr[i]-obj_corr[i])
	if correlations==[0.0]*len(correlations):
		objective='perfect' 
		return objective,correlations
	else:
		for i in range(len(correlations)):
			if correlations[:i+1]==(i+1)*[0.0] and correlations[i+1]!=0.0:
				index=i
	if 'index' in vars():
		objective=sum(map(abs,correlations))-weight*radius[index]
	else:
		objective=sum(map(abs,correlations))
	return objective,correlations

def equilibration_time(quantity,tolerance):
    t_max=len(quantity)
    diff=abs(average(quantity[:t_max/2])-average(quantity[t_max/2:]))
    if diff>tolerance:
        if t_max/2 >= 100:
            equilibration_time(quantity[t_max/2:],tolerance)
        else:
            print 'Not converged';sys.exit(1)
    else:
        print t_max/2

def skipchoice(alist,number):
    import random
    while 1:
        selected=random.sample(alist,number)
        j=0
        for i in selected:
            if i+1 not in selected and i-1 not in selected:
                j+=1
        if j==len(selected):
            return selected
            break

def read_clusters(order):
    subprocess.check_call('getclus > clusters.tmp',shell=True)
    clus_data=list(loadtxt('clusters.tmp')[:,0])
    os.remove('clusters.tmp')
    exclude_clus_num=0
    for i in range(order):
        exclude_clus_num+=clus_data.count(i)
    return clus_data.count(order),exclude_clus_num


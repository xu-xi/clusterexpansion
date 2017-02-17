#!/usr/bin/env python
import argparse,sys,os,shutil,subprocess,random,copy,numpy
from myfunc import initstr,occupy,ce_energy

parser=argparse.ArgumentParser(description='cluster expansion construction basing on given clusters and supercell. Configurations will be generated randomly to reach the target cv value. The default cv is 0, thus the construction process will never terminate. First-principles energies from VASP can be expanded in this version',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-c',default='clusters.out',dest='clusters',help="pre-defined clusters file")
parser.add_argument('-l',default='str.in',dest='lattice',help="the lattice file")
parser.add_argument('-v',type=float,default=0.0,dest='cv',help="the objective cross validation value")
parser.add_argument('-s',default='supercell.in',dest='supercell',help="the input supercell file")
#parser.add_argument('-n',type=int,default=5,dest='num',help="the minimum number of structures to bulid cluster expansion")
parser.add_argument('-p',default='energy',dest='property',help="the property to expand")

args=parser.parse_args()

if os.path.isfile('vasp.wrap'):
	pass
else:
	print 'ERROR: There is no vasp.wrap in current dir';sys.exit(1)

def random_structure_generator(index,energy_list):
	'generate structures randomly and call VASP to calculate its energy'
        while 1:
                str_file=file('str.out','r')
                new_str_file=file('new_str.out','w')
                for i in range(6):
                        new_str_file.write(str_file.readline())
                occupy(new_str_file,isublattices,iatomlist)
                vatomlist_trial=copy.deepcopy(vatomlist)
                i=random.randrange(len(vatomlist_trial))
                for i in range(len(vatomlist_trial)):
                        random.shuffle(vatomlist_trial[i])
                occupy(new_str_file,vsublattices,vatomlist_trial)
                new_str_file.close()
                str_file.close()
                os.remove('str.out')
                os.rename('new_str.out','str.out')
                if index==1:
                        break
                else:
                        new_str_energy=ce_energy('energy.eci','str.out')
                        if new_str_energy not in energy_list:
                                break
        os.mkdir(str(index))
	shutil.copy('str.out',str(index))
	os.chdir(str(index))
	subprocess.check_call("runstruct_vasp")
        os.chdir('../')
        '''
        subprocess.check_call("gzip -d OUTCAR.static.gz")
	subprocess.check_call("mpirun -np 8 vasp > vasp.out",shell=True)
	if subprocess.call("grep required OUTCAR > /dev/null",shell=True) == 1:
		os.mknod('error')
		celog.write('Calulation is not converged in %s\n' % (index))
		os.chdir('../')
	else:
		subprocess.call("grep without OUTCAR|tail -1|awk '{print $4}' > energy",shell=True)
                #subprocess.call("pv_gap.py OUTCAR|tail -1|awk '{print $3}' > bandgap",shell=True)
		os.chdir('../')
        '''

celog=file('ce.log','w')
index=0
#scan current directory to collect finished calculation
for item in os.listdir(os.environ['PWD']):
	fullpath=os.path.join(os.environ['PWD'],item)
	if os.path.isdir(fullpath):
		os.chdir(fullpath)
		if os.path.isfile('energy') and not os.path.isfile('error'):
			index+=1
			os.chdir('../')

celog.write(subprocess.check_output("date"))
celog.write(str(index)+' structures have alreadt been calulated in current directory.\n')
celog.flush()

try:
        energy_list=list(numpy.loadtxt('allenergy.out'))
except:
        energy_list=[]
        
#generate initial structure
isublattices,iatomlist,vsublattices,vatomlist,concentration,n=initstr(args.lattice,args.supercell)
celog.write('Initialization done.\n')
celog.flush()

#minimum cluster expansion
ce_min=int(subprocess.check_output("getclus | wc -l",shell=True))-int(index)
if ce_min>0:
        celog.write(str(ce_min)+' structures remain to be calculated in order to build the minimum cluster expansion...\n')
        for i in range(ce_min):
                index+=1
                random_structure_generator(index,energy_list)
                cv=float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" % (args.property),shell=True))
                energy_list=numpy.loadtxt('allenergy.out').tolist()
	        celog.write(subprocess.check_output("date"))
                celog.write('Cycle %s CV: %s\n' % (index,cv))
                celog.flush()

#build cluster expansion until convengence
while 1:
	index+=1
	random_structure_generator(index,energy_list)
	cv=float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" % (args.property),shell=True))
        energy_list=numpy.loadtxt('allenergy.out').tolist()
	celog.write(subprocess.check_output("date"))
	celog.write('Cycle %s CV: %s\n' % (index,cv))
	celog.flush()
        if os.path.isfile('stop'):
                subprocess.call('checkrelax')
                celog.close()
                os.remove('stop')
                break
	if cv<=args.cv:
		celog.write('Cluster expansion finished successfully')
		celog.close()
                subprocess.call('checkrelax')
		break
		

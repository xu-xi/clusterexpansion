#!/usr/bin/env python3
import argparse,sys,os,shutil,subprocess,random,copy,socket,ase.io,time
from celib import initstr,occupy,ce_energy
import numpy as np

def vasp_no_error(index):
    'check if current calculation is converged or not'
    os.chdir(str(index))

    if os.path.isfile('vasp.out.relax'):
        if subprocess.call('grep required vasp.out.relax',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))
            return 0
    elif os.path.isfile('vasp.out.static'):
        if subprocess.call('grep required vasp.out.static',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))
            return 0
    elif os.path.isfile('vasp.out'):
        if subprocess.call('grep required vasp.out',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))
            return 0
    else:
        return 0

def generate_wien2k_machines(processor_number=4):
    'generate .machines file for wien2k'
    hostname=socket.gethostname()
    para_file=open('.machines','w')
    para_file.write('#\nlapw0:')
    for i in range(processor_number):
        para_file.write(hostname+' ')
    for i in range(processor_number):
        para_file.write('\n1:'+hostname)
    para_file.write('\n\ngranularity:1\nextrafine:1')
    para_file.close()

def run_wien2k_scf():
    'call WIEN2k to do scf calculation and one-shot mBJ'
    os.mkdir('scf')
    shutil.copy('CONTCAR.static','./scf/scf')
    os.chdir('./scf')
    structure = ase.io.read('CONTCAR.static',format='vasp')
    ase.io.write('scf.struct',structure,format='struct') 
    generate_wien2k_machines()
    #subprocess.check_call('x sgroup',shell=True)
    #shutil.copy('scf.struct_sgroup','scf.struct')
    subprocess.check_call('init_lapw -b -ecut -8.0 -numk 100',shell=True)
    subprocess.check_call('run_lapw -p -ec 0.00001',shell=True)
    subprocess.check_call('w2k_save -d mbj',shell=True)
    subprocess.check_call('w2k_clean -p',shell=True)
    shutil.copy('./.machines','./mbj')
    os.chdir('./mbj')
    subprocess.check_call('w2k_tbbj -p',shell=True)


def random_structure_generator(lattice,supercell,index,struct_number,mbj=False):
    'generate structures randomly and call VASP to calculate its energy. The inner product of two cluster functions is used to judge their structural similarity'
    t0 = time.clock()
    while 1:
        initstr(lattice,supercell)
        if struct_number >= 1:
            corr_func = np.loadtxt('allcorr.out',ndmin=2)
            new_str = np.array(subprocess.check_output('corrdump -c -s=str.out',shell=True).split(),dtype=float)
            ip = []
            for i in range(len(corr_func)):
                ip.append(np.dot(corr_func[i],new_str)/(np.linalg.norm(corr_func[i])*np.linalg.norm(new_str)))
            if 1.0 - max(ip) > 0.01:
                break
            elif time.clock() - t0 > 100:
                raise RuntimeError('Could not find new structure. Please enlarge the supercell.')
        else:
            break
                
    while 1:
        try:
            os.mkdir(str(index))
        except:
            index += 1
        else:
            break

    shutil.copy('str.out',str(index))
    os.chdir(str(index))
    subprocess.check_call("runstruct_vasp")
    if mbj==True:
        run_wien2k_scf()
    os.chdir('../')



def Main(ArgList):
    parser=argparse.ArgumentParser(description=
'''
Cluster expansion construction for given clusters and supercell. Configurations will be generated randomly to reach the target cv value. Note that the default cv is 0.001 and the default structure number is 100. Change these default values by needs. create an empty file named \'stop\' can cleanly finish running job.

Four input files are required: lat.in, str.in, supercell.in and vasp.wrap, where lat.in, clusters.out and vasp.wrap are just the same as ATAT.

An example of str.in:

4.100000 0.000000 0.000000
0.000000 4.100000 0.000000
0.000000 0.000000 4.100000
1.000000 0.000000 0.000000
0.000000 1.000000 0.000000
0.000000 0.000000 1.000000
0.000000 0.500000 0.500000 O=2/3,N=1/3
0.500000 0.000000 0.500000 O=2/3,N=1/3
0.500000 0.500000 0.000000 O=2/3,N=1/3
0.000000 0.000000 0.000000 Ba
0.500000 0.500000 0.500000 Ta

An example of supercell.in (2*2*2 simple supecell):

2 0 0
0 2 0
0 0 2
''',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-l',default='str.in',dest='lattice',help="the lattice file")
    parser.add_argument('-v',type=float,default=0.001,dest='cv',help="the objective cross validation value")
    parser.add_argument('-n',type=int,dest='number',default=100,help="the maximum mumber of structures to fit")
    parser.add_argument('-s',default='supercell.in',dest='supercell',help="the input supercell file")
    parser.add_argument('-p',default='energy',dest='property',help="the property to expand")
    parser.add_argument('--mbj',action='store_true',dest='mbj',help="call WIEN2k to calculate mBJ band gaps")
    #parser.add_argument('--version',action='version',version='2017.2.23',help="output the version of the program")
    args=parser.parse_args()

    index_number = 0
    struct_number = 0
    wd = os.getcwd() #working directory

    #if not os.path.isfile('clusters.out'):
    subprocess.check_call('corrdump -clus -2=6 -3=5 -4=4',shell=True)

    #scan current directory to collect finished calculation
    for item in os.listdir(wd):
        fullpath = os.path.join(wd,item)
        if os.path.isdir(fullpath):
            index_number += 1
            os.chdir(fullpath)
            if os.path.isfile('energy') and not os.path.isfile('error'):
                struct_number += 1
            os.chdir(wd)
       
    celog = open('ce.log','w')

    if struct_number != 0:
        celog.write(str(struct_number)+' structures have already been calulated in current directory.\n')
        subprocess.check_call('clusterexpand -e -cv energy',shell=True)
        celog.flush()

    #minimum cluster expansion
    '''
    if args.clusters!=None:
        ce_min=int(subprocess.check_output("getclus | wc -l",shell=True))-int(struct_number)
        if ce_min>0:
            celog.write(str(ce_min)+' structures remain to be calculated in order to build the minimum cluster expansion...\n')
            for i in range(ce_min):
                random_structure_generator(args.lattice,args.supercell,index,struct_number,energy_list)
                cv=float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" % (args.property),shell=True))
                energy_list=np.loadtxt('allenergy.out').tolist()
                celog.write(subprocess.check_output("date"))
                celog.write('Cycle %s CV: %s\n' % (struct_number,cv))
                celog.flush()
                if os.path.isfile('stop'):
                    celog.close()
                    os.remove('stop')
                    break
   '''

    #build cluster expansion until convengence
    while 1:
        index_number += 1
        random_structure_generator(args.lattice,args.supercell,index_number,struct_number,args.mbj)
        if vasp_no_error(index_number) != 0:
            struct_number += 1
        os.chdir(wd)
        cv = float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" %(args.property),shell=True))
        energy_list = np.loadtxt('allenergy.out',ndmin=1)
        celog.write(subprocess.check_output("date"))
        celog.write('Cycle %s CV: %s\n' % (struct_number,cv))
        celog.flush()
        if os.path.isfile('stop') or cv < args.cv:
            celog.close()
            try:
               os.remove('stop')
            except:
                pass
            break
        if args.number != None:
            if index_number >= args.number:
                celog.close()
                break
		
if __name__=="__main__":
    Main(sys.argv);

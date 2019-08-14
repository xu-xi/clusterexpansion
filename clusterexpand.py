#!/usr/bin/env python
import argparse,sys,os,shutil,subprocess,random,copy,numpy,socket,ase.io
from celib import initstr,occupy,ce_energy
from clusters_opt import clusters_optimizer

def vasp_no_error(index):
    'check if current calculation is converged or not'
    os.chdir(str(index))
    if os.path.isfile('vasp.out.relax'):
        if subprocess.call('grep required vasp.out.relax',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))
            os.chdir('../')
            return 0
    elif os.path.isfile('vasp.out.static'):
        if subprocess.call('grep required vasp.out.static',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))
            os.chdir('../')
            return 0
    elif os.path.isfile('vasp.out'):
        if subprocess.call('grep required vasp.out',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))
            os.chdir('../')
            return 0
    else:
        return 0

def generate_wien2k_machines(processor_number=4):
    'generate .machines file for wien2k'
    hostname=socket.gethostname()
    para_file=file('.machines','w')
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


def random_structure_generator(lattice,supercell,index,struct_number,energy_list,mbj=False):
    'generate structures randomly and call VASP to calculate its energy'
    while 1:
        initstr(lattice,supercell)
        if index<=2:
            break
        else:
            new_str_energy=float(subprocess.check_output('corrdump -c -eci=energy.eci -s=str.out',shell=True))
            if new_str_energy not in energy_list:
                break
    while 1:
        try:
            os.mkdir(str(index))
        except:
            index+=1
        else:
            break

    shutil.copy('str.out',str(index))
    os.chdir(str(index))
    subprocess.check_call("runstruct_vasp")
    if mbj==True:
        run_wien2k_scf()
    os.chdir('../')


def Main(ArgList):
    parser=argparse.ArgumentParser(description="Cluster expansion construction basing on given clusters and supercell. Configurations will be generated randomly to reach the target cv value. Note that the default cv is 0, thus the construction process will never terminate by default. Set a finite CV by -v or specify -n or create an empty file named \'stop\' to finish current job. \n\nFour input files are required: lat.in, str.in, supercell.in and vasp.wrap.",

                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c',default='clusters.out',dest='clusters',help="pre-defined clusters file")
    parser.add_argument('-l',default='str.in',dest='lattice',help="the lattice file")
    parser.add_argument('-v',type=float,default=0.0,dest='cv',help="the objective cross validation value")
    parser.add_argument('-n',type=int,dest='number',help="the maximum mumber of structures to fit")
    parser.add_argument('-s',default='supercell.in',dest='supercell',help="the input supercell file")
    parser.add_argument('-p',default='energy',dest='property',help="the property to expand")
    parser.add_argument('--mbj',action='store_true',dest='mbj',help="call WIEN2k to calculate mBJ band gaps")
    #parser.add_argument('--version',action='version',version='2017.2.23',help="output the version of the program")
    args=parser.parse_args()

    index_number=0
    struct_number=0
    #scan current directory to collect finished calculation
    for item in os.listdir(os.getcwd()):
        fullpath=os.path.join(os.getcwd(),item)
        if os.path.isdir(fullpath):
            index_number+=1
            os.chdir(fullpath)
            if os.path.isfile('energy') and not os.path.isfile('error'):
                struct_number+=1
            os.chdir('../')

    #just optimize the choice of clusters and do not construct cluster expansion
    #if args.opt_only==True:
    #    cluster_optimizer()
    #    sys.exit(1)
        
    celog=file('ce.log','w')
    celog.write(subprocess.check_output("date"))
    celog.write(str(struct_number)+' structures have already been calulated in current directory.\n')
    celog.flush()

    #if two structures have different ce energies, they are definitly different, therefore the energy_list is used to avoid repeated calculation 
    try:
        energy_list=list(numpy.loadtxt('allenergy.out'))
    except:
        energy_list=[]

    #minimum cluster expansion
    '''
    if args.clusters!=None:
        ce_min=int(subprocess.check_output("getclus | wc -l",shell=True))-int(struct_number)
        if ce_min>0:
            celog.write(str(ce_min)+' structures remain to be calculated in order to build the minimum cluster expansion...\n')
            for i in range(ce_min):
                random_structure_generator(args.lattice,args.supercell,index,struct_number,energy_list)
                cv=float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" % (args.property),shell=True))
                energy_list=numpy.loadtxt('allenergy.out').tolist()
                celog.write(subprocess.check_output("date"))
                celog.write('Cycle %s CV: %s\n' % (struct_number,cv))
                celog.flush()
                if os.path.isfile('stop'):
                    celog.close()
                    os.remove('stop')
                    break
    '''
    if struct_number < 5: 
        for i in range(5-struct_number):
            index_number+=1
            random_structure_generator(args.lattice,args.supercell,index_number,struct_number,energy_list,args.mbj)
            if vasp_no_error(index_number)!=0:
                struct_number+=1
            cv=float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" % (args.property),shell=True))
            energy_list=numpy.loadtxt('allenergy.out').tolist()
            celog.write(subprocess.check_output("date"))
            celog.write('Cycle %s CV: %s\n' % (struct_number,cv))
            celog.flush()
            if os.path.isfile('stop'):
                celog.close()
                os.remove('stop')
                sys.exit(1)
                break

    #build cluster expansion until convengence
    while 1:
        index_number+=1
        random_structure_generator(args.lattice,args.supercell,index_number,struct_number,energy_list,args.mbj)
        if vasp_no_error(index_number)!=0:
            struct_number+=1
        #cv=clusters_optimizer(args.property,struct_number-1,args.step)
        cv=float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" %(args.property),shell=True))
        energy_list=numpy.loadtxt('allenergy.out').tolist()
        celog.write(subprocess.check_output("date"))
        celog.write('Cycle %s CV: %s\n' % (struct_number,cv))
        celog.flush()
        if os.path.isfile('stop') or cv<=args.cv:
            celog.close()
            os.remove('stop')
            break
        if args.number!=None:
            if index_number>=args.number:
                celog.close()
                break
		
if __name__=="__main__":
    Main(sys.argv);

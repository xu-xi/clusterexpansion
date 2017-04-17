#!/usr/bin/env python
import argparse,sys,os,shutil,subprocess,random,copy,numpy
from myfunc import initstr,occupy,ce_energy

def check_vasp_error(index):
    'check if current calculation is converged or not'
    if os.path.isfile('vasp.out.relax'):
        if subprocess.call('grep required vasp.out.relax',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))
    elif os.path.isfile('vasp.out.static'):
        if subprocess.call('grep required vasp.out.static',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))
    elif os.path.isfile('vasp.out'):
        if subprocess.call('grep required vasp.out',shell=True)!=0:
            os.mknod('error')
            sys.stderr.write('ERROR in %s, check the output' % (index))

def random_structure_generator(lattice,supercell,index,energy_list):
    'generate structures randomly and call VASP to calculate its energy'
    while 1:
        initstr(lattice,supercell)
        if index<=2:
            break
        else:
            new_str_energy=ce_energy('energy.eci','str.out')
            if new_str_energy not in energy_list:
                break
    os.mkdir(str(index))
    shutil.copy('str.out',str(index))
    os.chdir(str(index))
    subprocess.check_call("runstruct_vasp")
    check_vasp_error(index)
    os.chdir('../')

def Main(ArgList):
    parser=argparse.ArgumentParser(description='Cluster expansion construction basing on given clusters and supercell. Configurations will be generated randomly to reach the target cv value. Note that the default cv is 0, thus the construction process will never terminate by default. Set a finite CV by -v or specify -n or create an empty file named \'stop\' to finish current job. Five input files are required: lat.in,str.in,supercell.in,clusters.out,vasp.wrap',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c',default='clusters.out',dest='clusters',help="pre-defined clusters file")
    parser.add_argument('-l',default='str.in',dest='lattice',help="the lattice file")
    parser.add_argument('-v',type=float,default=0.0,dest='cv',help="the objective cross validation value")
    parser.add_argument('-n',type=int,dest='number',help="the number of maximum structures to fit")
    parser.add_argument('-s',default='supercell.in',dest='supercell',help="the input supercell file")
    parser.add_argument('-p',default='energy',dest='property',help="the property to expand")
    #parser.add_argument('--version',action='version',version='2017.2.23',help="output the version of the program")

    args=parser.parse_args()

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
    celog.write(str(index)+' structures have already been calulated in current directory.\n')
    celog.flush()

    #if two structures have different ce energies, they are definitly different, therefore energy_list is used to avoid repeated calculation 
    try:
        energy_list=list(numpy.loadtxt('allenergy.out'))
    except:
        energy_list=[]

    #minimum cluster expansion
    ce_min=int(subprocess.check_output("getclus | wc -l",shell=True))-int(index)
    if ce_min>0:
        celog.write(str(ce_min)+' structures remain to be calculated in order to build the minimum cluster expansion...\n')
        for i in range(ce_min):
            index+=1
            random_structure_generator(args.lattice,args.supercell,index,energy_list)
            cv=float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" % (args.property),shell=True))
            energy_list=numpy.loadtxt('allenergy.out').tolist()
            celog.write(subprocess.check_output("date"))
            celog.write('Cycle %s CV: %s\n' % (index,cv))
            celog.flush()
            if os.path.isfile('stop'):
                celog.close()
                os.remove('stop')
                break

    #build cluster expansion until convengence
    while 1:
        index+=1
        random_structure_generator(args.lattice,args.supercell,index,energy_list)
        cv=float(subprocess.check_output("clusterexpand -e -cv %s | tail -1" % (args.property),shell=True))
        energy_list=numpy.loadtxt('allenergy.out').tolist()
        celog.write(subprocess.check_output("date"))
        celog.write('Cycle %s CV: %s\n' % (index,cv))
        celog.flush()
        if os.path.isfile('stop') or cv<=args.cv:
            celog.close()
            os.remove('stop')
            break
        if args.number!=None:
            if index>=args.number:
                celog.close()
                break
		
if __name__=="__main__":
    Main(sys.argv);

#!/usr/bin/env python
import ase.io,os,math,cmath,sys,re,subprocess,argparse
import numpy as np
import matplotlib.pyplot as plt
from celib import read_cluster_function,read_quantity,plot_fit,calc_cv
from sklearn.metrics import mean_squared_error

class Reciprocal_Space_Cluster_Expansion():
    def __init__(self):
        try:
            self.lattice = np.loadtxt('lat.in',usecols=(0,1,2))
        except IOError,error_msg:
            print error_msg; sys.exit(1)

        self.cell = np.dot(self.lattice[:3,:].T,self.lattice[3:6]) #lattice vectors for the primitive cell 
        self.G = 2*math.pi*np.linalg.inv(self.cell.T) #reciprocal lattice vectors for the primitive cell

        self.mixed_atoms = str(subprocess.check_output('grep , lat.in | head -1 | awk \'{print $4}\'',shell=True)).strip().split(',') #get the chemical symbols for mixed occupied atoms.
        print 'Mixed occupied atoms:',self.mixed_atoms
    
    def set_kpoints(self,klist):
        'set k mesh'
        kpoints = []
        for k in klist:
            kpoints.append(np.dot(self.G.T,k))

        return kpoints

    def get_spin_variables(self,strfile):
        'get the spin variables and positions of the mixed occupied atoms'
        self.atoms = ase.io.read(strfile,format = 'vasp') #use the relaxed structure for a more precise discription?
        atom_list = self.atoms.get_chemical_symbols()
        spin_variables = []
        positions = []
        for i in range(len(atom_list)):
            if atom_list[i] == self.mixed_atoms[0]: #the spin variable of 1st one is defined to be -1
                spin_variables.append(-1)
                positions.append(self.atoms.get_positions()[i])
            elif atom_list[i] == self.mixed_atoms[1]: #the spin variable of 2nd one is defined to be +1
                spin_variables.append(1)
                positions.append(self.atoms.get_positions()[i])
            else: #other atoms are ignored
                pass
        #print spin_variables,positions
        return spin_variables,positions

    def get_spin_wave(self,kpoints,strfile):
        'get the Fourier transfrom of the spin variables of a structure'
        #positions = self.get_positions()
        spin_variables, positions = self.get_spin_variables(strfile)
        M = len(spin_variables)
        N = len(kpoints)
        F = np.zeros(shape=(M,N),dtype=complex) #the Fourier transform matrix 
        for i in range(M):
            for j in range(N):
                F[i,j] = cmath.exp(np.dot(kpoints[j],positions[i])*1j)
        if M > N:
            print 'WARNING: The number of k points is less than the number of sites.'
        elif M == N:
            if abs(np.linalg.det(F)) < 1e-1:
                print 'WARNING: singular Fourier matrix.'

        spin_wave = np.dot(F.T,spin_variables)/M
        S = map(lambda x: abs(x)**2,spin_wave)
        #print sum(S)
        return S


def Main(ArgList):
    parser=argparse.ArgumentParser(description='Reciprocal-space Cluster Expansion. Works only for binary one-site unitcell.',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('-f',type=str,default='POSCAR.static',dest='strfile',help="The structure file to read")
    parser.add_argument('-p',type=str,default='energy',dest='property',help="The property to expand")
    parser.add_argument('-k',type=int,default=None,dest='knum',help="The number of k points")
    parser.add_argument('-b',type=str,default=None,dest='B',help="The B atom")
    args=parser.parse_args()

    rce = Reciprocal_Space_Cluster_Expansion()

    #set k mesh
    klist = []
    for i in [0,0.5]:
        for j in [0,0.5]:
            for k in [0,0.5]:
                klist.append([i,j,k])
    print 'k mesh',klist
    kpoints = rce.set_kpoints(klist)
    
    #get reciprocal cluster function
    cluster_function_in_kspace = []
    for item in os.listdir(os.getcwd()):
        fullpath = os.path.join(os.getcwd(),item)
        if os.path.isdir(fullpath):
            os.chdir(fullpath)
            spin_wave = rce.get_spin_wave(kpoints,args.strfile)
            cluster_function_in_kspace.append(spin_wave)
            os.chdir('../')

    cluster_function = read_cluster_function()
    X = np.hstack((cluster_function,cluster_function_in_kspace))
    print 'The shape of mixed cluster function:',X.shape

    if args.property == 'energy':
        y = read_quantity(args.property,average=True)
    else:
        y = read_quantity(args.property,average=False) #For band gaps

    eci = np.linalg.lstsq(X,y,rcond=None)[0]
    print 'ECI:',eci
    rmsd = math.sqrt(mean_squared_error(y,np.dot(X,eci)))
    cv = calc_cv(X,eci,y)
    plot_fit(y,np.dot(X,eci),cv,rmsd)
    plt.xlabel('Calculated %s(eV)' %(args.property))
    plt.ylabel('Fitted %s(eV)' %(args.property))
    plt.show()

if __name__=='__main__':
    Main(sys.argv)

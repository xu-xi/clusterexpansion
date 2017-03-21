#!/usr/bin/env python
import sys,os,itertools,shutil,argparse
from myfunc import *

parser=argparse.ArgumentParser(description='GENerate all STRuctures according to translation symmetry',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
parser.add_argument('-f','--file',default='str.in',dest='latfile',help="the input lattice file")
parser.add_argument('-s','--supercell',default='supercell.in',dest='supercell_file',help="the input supercell file")
args=parser.parse_args()

isublattices,iatomlist,vsublattices,vatomlist,concentrations,n=initstr(args.latfile,args.supercell_file)
strfile=file('str.out','r')
lattinfo=[]
weights={}
index={}
for i in range(6):
	lattinfo.append(map(float,strfile.readline().split()))
step=0
for vatomsites in itertools.product(*[enumer(x) for x in vatomlist]):
	strfile=file('str.out','w')
	listwrite(strfile,lattinfo)
	occupy(strfile,isublattices,iatomlist)
	occupy(strfile,vsublattices,list(vatomsites))
	strfile.close()
	energy0=ce_energy('eci.out','str.out')
	step+=1
	if energy0 in weights:
		weights[energy0]+=1
	else:
		weights[energy0]=1
		index[energy0]=step
		os.mkdir(str(step))
		shutil.copy('str.out',str(step))

datafile=file('data.out','w')
datafile.write('#index\tenergy\tweights\n')
energy_list=sorted(weights)
for i in energy_list:
	datafile.write('%s\t%s\t%s\n' %(index[i],i,weights[i]))
datafile.close()


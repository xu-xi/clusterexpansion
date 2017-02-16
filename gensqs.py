#!/usr/bin/env python
import sys,os,shutil,argparse,scitools.filetable
from myfunc import *

parser=argparse.ArgumentParser(description='GENerate the best SQS according to enumeration method',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
parser.add_argument('-f','--file',default='str.in',dest='latfile',help="the input lattice file ")
#parser.add_argument('-c','--corrfunc',default='corr_func.in',dest='corr_func_file',help="the input file defining correlation functions")
parser.add_argument('-s','--supercell',default='supercell.in',dest='supercell_file',help="the input supercell file ")
parser.add_argument('-w','--weight',type=int,default=1,dest='weight',help="Weight assigned to range of perfect correlation match in objective function ")
args=parser.parse_args()


sqslog=file('gensqs.log','w')
isublattices,iatomlist,vsublattices,vatomlist,concentration,n=initstr(args.latfile,args.supercell_file)
sqslog.write('Initialization done.\n')
sqslog.flush()
strfile=file('str.out','r')
lattinfo=[]
for i in range(6):
	lattinfo.append(map(float,strfile.readline().split()))
min_value=10

failure,output=commands.getstatusoutput('getclus > clusters.tmp')
if failure:
	print 'ERROR: check the cluster.out file or run the corrdump command to generate it';sys.exit(1)
cluster,radius,multi=scitools.filetable.read_columns(file('clusters.tmp'))
os.remove('clusters.tmp')
obj_corr=objective_correlation_functions(concentration,cluster)

for vatomsites in itertools.product(*[enumer(i) for i in vatomlist]):
	strfile=file('str.out','w')
	listwrite(strfile,lattinfo)
	occupy(strfile,isublattices,iatomlist)
	occupy(strfile,vsublattices,list(vatomsites))
	strfile.close()
	objective,correlations=calc_correlation_functions(args.weight,obj_corr,radius)
	#print objective
	if objective=='perfect':
		shutil.copy('str.out','bestsqs.out')
		sqslog.write('Objective_funtion = Perfect!\n')
		sqslog.write('Correlations_mismatched = ')
		for i in correlations:
			sqslog.write('%.5f\t' % (i))
		sqslog.write('\n')
		break
	elif objective<min_value:
		min_value=objective
		shutil.copy('str.out','bestsqs.out')
		sqslog.write('Objective_funtion = %.5f\n' % (objective))
		sqslog.write('Correlations_mismatched = ')
		for i in correlations:
			sqslog.write('%.5f\t' % (i))
		sqslog.write('\n')
		sqslog.flush()
sqslog.close()

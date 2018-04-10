#!/usr/bin/env python
import numpy,argparse
import matplotlib.pyplot as plt

parser=argparse.ArgumentParser(description='plot phase diagram according to the result of phb code',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
parser.add_argument('-f',default='phb.out',dest='phb',nargs='+',help="The output file of code phb, support multiple files")
args=parser.parse_args()

#read data file
for i in args.phb:
    datafile=numpy.loadtxt(i)
    T=list(datafile[:,0])
    x1=list((datafile[:,2]+1)/2)
    x2=list((datafile[:,3]+1)/2)

    T+=list(reversed(T))
    x=x1+list(reversed(x2))

    plt.plot(x,T,linewidth=2.5)
    #plt.scatter(x,T)
    #plt.scatter(x1,T)
    #plt.scatter(x2,T)

plt.xlabel('Concentration x')
plt.ylabel('Temperature/K')
plt.xlim(0,1)
plt.ylim(ymin=0)
plt.show()

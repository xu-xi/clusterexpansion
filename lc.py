#!/usr/bin/env python
import argparse,sys,subprocess
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
from sklearn.model_selection import RepeatedKFold
from celib import plot_learning_curve

parser = argparse.ArgumentParser(description='Plot learning curves of least-squares fitting algorithm for a given cluster expansion model.',
                formatter_class = argparse.ArgumentDefaultsHelpFormatter) 
parser.add_argument('-p',default='energy',dest='property',help="The property to expand")
args = parser.parse_args()

subprocess.check_call('corrdump -c -s=allstr.out > allcorr.out',shell=True)

X = np.loadtxt('allcorr.out')
y = np.loadtxt('all%s.out' %(args.property))

str_num, clus_num = X.shape

if clus_num >= str_num*0.33:
    print "Reduce the number of clusters or increase the number of structures.";sys.exit(1)

CV = RepeatedKFold(5,10)
lstsq = linear_model.LinearRegression()
plot_learning_curve(lstsq,X,y,cv=CV)


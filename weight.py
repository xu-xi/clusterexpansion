#!/usr/bin/env python
import numpy
import matplotlib.pyplot as plt

datafile=numpy.loadtxt('data.out')
index=list(datafile[:,0])
ce_energy=list(datafile[:,1])
weights=list(datafile[:,2])

print int(index[ce_energy.index(min(ce_energy))])
plt.xlabel('CE Energy')
plt.ylabel('Geometrical Probability')
#plt.ylim(0,0.003)
plt.scatter(ce_energy,weights/sum(weights),c='k',s=20)
plt.savefig('bton-weights-energy.png')

#!/usr/bin/env/python

import sys
import os
import commands
import numpy as np
import matplotlib.pyplot as plt
import math

f = open('Errors_1.dat')
file1  = f.readlines()
f.close()
f = open('Errors_2.dat')
file2 = f.readlines()
f.close()
f = open('Errors_3.dat')
file3 = f.readlines()
f.close()


N1 = len(file1)
N2 = len(file2)
N3 = len(file3)

h_1 = np.zeros(N1)
e_p_1 = np.zeros(N1)
e_E_1 = np.zeros(N1)
order_1 = np.zeros(N1)

h_2 = np.zeros(N2)
e_p_2 = np.zeros(N2)
e_E_2 = np.zeros(N2)

h_3 = np.zeros(N3)
e_p_3 = np.zeros(N3)
e_E_3 = np.zeros(N3)


for i in range(0,N1):
	h_1[i] = float(file1[i].split()[0])
	e_p_1[i] = float(file1[i].split()[1])
	e_E_1[i] = float(file1[i].split()[2])


for i in range(0,N1):
	h_2[i] = float(file2[i].split()[0])
	e_p_2[i] = float(file2[i].split()[1])
	e_E_2[i] = float(file2[i].split()[2])

for i in range(0,N3):
	h_3[i] = float(file3[i].split()[0])
	e_p_3[i] = float(file3[i].split()[1])
	e_E_3[i] = float(file3[i].split()[2])



slope, intersept = np.polyfit(np.log(h_1), np.log(e_p_1), 1)
print('k=0, primary, slope = ' + str(slope))
slope, intersept = np.polyfit(np.log(h_1), np.log(e_E_1), 1)
print('k=0, auxillary, slope = ' + str(slope))
slope, intersept = np.polyfit(np.log(h_2), np.log(e_p_2), 1)
print('k=1, primary, slope = ' + str(slope))
slope, intersept = np.polyfit(np.log(h_2), np.log(e_E_2), 1)
print('k=1, auxillary, slope = ' + str(slope))
slope, intersept = np.polyfit(np.log(h_3), np.log(e_p_3), 1)
print('k=2, primary, slope = ' + str(slope))
slope, intersept = np.polyfit(np.log(h_3), np.log(e_E_3), 1)
print('k=2, auxillary, slope = ' + str(slope))

plt.clf()
plt.title('$L^{2}$-errors for MFEM on Poisson')
plt.loglog(h_1, e_p_1, 'b', linewidth=2, label='$k=1$,  primary')
plt.loglog(h_1, e_E_1, 'r', linewidth=2, label='$k=1$,  auxillary')
plt.loglog(h_2, e_p_2, 'g', linewidth=2, label='$k=2$,  primary')
plt.loglog(h_2, e_E_2, 'm', linewidth=2, label='$k=2$,  auxillary')
plt.loglog(h_3, e_p_3, 'k', linewidth=2, label='$k=3$,  primary')
plt.loglog(h_3, e_E_3, 'c', linewidth=2, label='$k=3$,  auxillary')

plt.legend(('$k=1$, Potential', '$k=1$, Electric Field', '$k=2$, Potential', 
						'$k=2$, Electric Field', '$k=3$, primary', '$k=3$, auxillary'), ncol=2, loc=4,
						fontsize=15)

#plt.legend(('$k=0$, Density', '$k=0$, Current', '$k=1$, Density', 
#						'$k=1$, Current', '$k=2$, Density', '$k=2$, Current'), ncol=2, loc=4,
#						fontsize=15)
plt.ylabel('$\log ( \Vert \cdot \Vert_{L^{2}} )$')
plt.xlabel('$\log(h)$')

plt.savefig('errors.eps')
plt.show()


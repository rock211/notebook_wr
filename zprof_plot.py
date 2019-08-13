import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
import copy
from mpl_toolkits import axes_grid1
unit = pa.set_units(muH=1.4271)

from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle
labell = ('No ICM','P1', 'P3','P3h', 'P7','P7h','P14')  # r'No ICM',
#C = ('gray', 'mediumturquoise', 'dodgerblue','mediumblue' ,'goldenrod','salmon', 'firebrick','darkmagenta','goldenrod','royalblue','crimson') # 'plum','orchid','purple'
C = ('gray', 'lightskyblue', 'dodgerblue','mediumblue' ,'goldenrod','salmon', 'firebrick')
Model = [0,8.63*1e3,3.46*1e4,3.46*1e4,6.92*1e4,6.92*1e4,1.38*1e5]
lw = (2.7, 1.3, 2.4, 3.5, 2.4, 3.5, 1.3)
crit = 94
k=1
z = range(895)

for j in (1,2,3,4,5,6):

    if j == 5 or j == 6:
        stop = 474
    else:
        stop = 499

    frac = np.genfromtxt('./proj/fraction_%s.txt' % labell[j])
    plt.plot(range(250, stop) * unit['time'], frac, c=C[j], linewidth=lw[j], label=labell[j])
plt.ylim(0, 1)
plt.xlim(250 * unit['time'].value, 499 * unit['time'].value)
plt.tick_params(labelsize=12, direction='in')
plt.ylabel('Fraction', fontsize=13)
plt.xlabel('Time [Myr]', fontsize=13)
plt.legend(loc=0, fontsize=11.5)
# plt.savefig('D:/yeongu/plots/paperplot/new/surf_frac_all_allphase.png')# % labell[j])
# plt.close()
plt.show()
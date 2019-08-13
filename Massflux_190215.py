import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
from scipy.interpolate import spline
import pyathena as pa
from multiprocessing import Pool
from matplotlib.ticker import MultipleLocator

from mpl_toolkits import axes_grid1

from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl

unit = pa.set_units(muH=1.4271)
print(unit)
# print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
# print unit['density']
kb = 1.3806504 * 1e-16  # boltzmann constant erg/K
volpercell = 7168. * 1024 * 1024 / (128 * 128 * 896)
vpc = volpercell

# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

simid_t = (
'RPS_8pc_noICM_newacc', 'RPS_8pc_ICM0_newacc', 'RPS_8pc_ICM1_newacc', 'RPS_4pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc',
'RPS_4pc_ICM2_newacc', 'RPS_8pc_ICM3_newacc')  # 'MHD_8pc_new' ,
# labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # r'No ICM',
labell = ('No ICM', 'P1', 'P3', 'P3h', 'P7', 'P7h', 'P14')
labelll = ('Cold', 'Unstable', 'Warm', 'Ionized', 'Hot')
C = ('goldenrod', 'royalblue', 'firebrick')
C2 = ('darkblue', 'orange', 'goldenrod', 'red', 'firebrick')
C3 = ('darkred', 'red', 'salmon')
S = (':', '--', '-')
hh = [0.00628, 0.00628, 0.00628, 0.021, 0.1585, 0.1585, 0.25]
ylim = [0.00628, 0.00628, 0.00628, 0.021, 0.1585, 0.1585, 0.25]
ylim_c = [0.14, 0.14, 0.14, 0.35, 11, 11, 11]
inflow = [0, 0.005002 / 2.828, 0.005002, 0.005002, 0.005002 * 2, 0.005002 * 2, 0.005002 * 2.828]
k = 0
# overplot Starformation rate of three different simulations
# plt.figure(figsize=(14,8))
# fig =plt.figure(figsize=(8.5,12))
jj = (0,3,5)
if len(jj)==3:
    ysize=10
    bot = 0.07
elif len(jj)==4:
    ysize=15
    bot = 0.055
elif len(jj)==5:
    ysize=19
    bot = 0.04
elif len(jj) == 6:
    ysize = 23
    bot = 0.025
fig, axs = plt.subplots(len(jj), 1, figsize=(13, ysize), sharex=True)
icm_c = '#03d803ff'
stop = 499
k = 0
for j in jj:  # range(1,7) :

    cuw = np.genfromtxt('./flux/Mom_cuw_%s.txt' % labell[j])
    ion = np.genfromtxt('./flux/Mom_i_%s.txt' % labell[j])
    hot = np.genfromtxt('./flux/Mom_h_%s.txt' % labell[j])

    if j!=0:
        icm = np.genfromtxt('./flux/Mom_icm_%s.txt' % labell[j])
    if j == 5 or j == 6:
        time = np.arange(250 * unit['time'].value, 473 * unit['time'].value, unit['time'].value)
        time = time[0:-1]
    else:
        time = np.arange(250 * unit['time'].value, stop * unit['time'].value, unit['time'].value)

    #### cumsum overlap ####
    ax1 = axs[k].twinx()
    ax1.plot(time, np.cumsum(cuw), c=C2[0], lw=1.5, ls='--')
    ax1.plot(time, np.cumsum(ion), c='salmon', lw=1.5, ls='--') #'#ef6fffff'
    ax1.plot(time, np.cumsum(hot), c=C2[3], lw=1.5, ls='--')
    ax1.plot([],[],ls='-',label='Mass flux',c='k')
    ax1.plot([],[],ls='--',label='Cumulative',c='k')
    ax1.set_ylim(-ylim_c[j] * 0.02, ylim_c[j])
    ax1.tick_params(which='major', direction='in', labelsize=20)
    ax1.locator_params(axis='y', nbins=5)
    if len(jj)==3 or len(jj)==4:
        if k==1:
            ax1.legend(loc='upper left',fontsize=15)
    elif len(jj)==5:
        if k==2:
            ax1.legend(loc='upper left',fontsize=15)
    # if k == 1:
    #    ax1.set_ylabel('Cumulative mass [M$_{\odot}$ pc$^{-2}$]',fontsize=24)
    ########################
    print cuw.max(), ion.max(), hot.max()
    #### Mass flux #########
    axs[k].plot(time, cuw, label='Cold/Warm', color=C2[0], lw=2.5)
    axs[k].plot(time, ion, label='Ionized', color='salmon', lw=2.5) #
    axs[k].plot(time, hot, label='Hot', color=C2[3], lw=2.5)
    if j != 0:
        axs[k].axhline(inflow[j], label='ICM inflow', ls='-.', color='k')
        axs[k].plot(time, icm, label='ICM outflow', color=icm_c, ls='-.')

    axs[k].text(482, hh[j] * 0.9, '%s' % labell[j], fontsize=20, ha='right', va='bottom')
    axs[k].set_ylim(-ylim[j] * 0.02, ylim[j])
    axs[k].tick_params(which='major', direction='in', labelsize=20)

    #########################
    #plt.tight_layout()
    k = k + 1

if len(jj)==3 or len(jj)==4:
    axs[1].legend(loc='center left',fontsize=15)
elif len(jj)==5:
    axs[2].legend(loc='center left',fontsize=15)

# axs[1].set_ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]',fontsize=24)
axs[k-1].set_xlabel('Time [Myr]', fontsize=20)
plt.xlim(250 * unit['time'].value, 498 * unit['time'].value)

######### y label ###########
ax = fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
ax.tick_params(which='both', labelcolor='none', top=False, bottom=False, left=False, right=False)
ax.set_ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]', fontsize=24, labelpad=45)

ax2 = fig.add_subplot(111, sharex=ax, frameon=False)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.tick_params(which='both', labelcolor='none', top=False, bottom=False, left=False, right=False)
ax2.set_ylabel('Cumulative mass [M$_{\odot}$ pc$^{-2}$]', fontsize=24, labelpad=35)
#############################
plt.subplots_adjust(bottom=bot, top=0.99, hspace=.01, left=0.17, right=0.86)
plt.savefig('D:/yeongu/plots/paperplot/new/massflux_test_%s_p.png' % len(jj),dpi=300) # ,transparent=True
plt.show()



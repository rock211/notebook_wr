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

from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl

# In[8]:

unit = pa.set_units(muH=1.4271)
print(unit)
# print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
# print unit['density']
kb = 1.3806504 * 1e-16  # boltzmann constant erg/K / erg = g cm2/s2K
vpc = 7168. * 1024 * 1024 / (128 * 128 * 896)  # volume per cell

z = np.arange(-3588, 3588, 8)
g = 4.5181 * 1e-30  # gravitational constant : pc3/solmass*s2
sig_star = 42  # solmass/pc2
z_s = 245  # pc
r_0 = 8000  # pc
rho_dm = 0.0064  # solmass/pc3
km = 3.24078 * 1e-14  # 1km in parsec
cm = 3.24078 * 1e-19  # 1cm in parsec
gram = 5.02785 * 1e-34  # 1g in solar mass

gext = (2. * np.pi * g * sig_star * (z ) / ((z ) ** 2 + z_s ** 2) ** (0.5) + 4. * np.pi * g * rho_dm * (
            (z ) / (1 + (z ) ** 2 / r_0 ** 2)))  # pc/s2
gext = gext/km/km # km/s2

# plt.plot((z-3584.)/1000,g_z)
# plt.show()


stop = 500

simid_t = (
'RPS_8pc_noICM_newacc', 'RPS_8pc_ICM0_newacc', 'RPS_8pc_ICM1_newacc', 'RPS_4pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc',
'RPS_4pc_ICM2_newacc', 'RPS_8pc_ICM3_newacc')
labell = ('No ICM', 'P1', 'P3', 'P3h', 'P7', 'P7h', 'P14')  # r'No ICM',
Model = [8.63 * 1e3, 3.46 * 1e4, 6.92 * 1e4, 1.38 * 1e5]
Modell = [r'8.63*1e3', r'$3.46x10^4$', r'$6.92x10^4$', r'1.38*1e5']
#C = ('dimgray', 'lightskyblue','mediumblue' ,'','goldenrod','', 'firebrick','darkmagenta','goldenrod','royalblue','crimson') # 'plum','orchid','purple'
C = ('gray', 'lightskyblue', 'dodgerblue','mediumblue','salmon', 'firebrick' ,'goldenrod')

lw = (2.7, 1.3, 2.4, 3.5, 2.4, 3.5, 1.3)

S = ('-.', '--', '-')

# overplot Starformation rate of three different simulations
k = 1

crit_0 = 0.8146 # solar mass per pc2
crit_1 = 3.2576 # K/cm3
crit_2 = 6.5152
crit_3 = 13.034
crit = (0,0.8146,3.2576,3.2576,6.5152,6.5152,13.034)# solar mass per pc2
#plt.figure(figsize=(6, 15))
for j in (2,3,4,5):
    basedir = 'G:/yeongu/'
    if j == 5 or j == 6:
        stop = 474
    else:
        stop = 499
    simid = simid_t[j]
    #frac = []
    #plt.figure(figsize=(6, 10))
    print j
    frac = np.genfromtxt('./proj/fraction_%s.txt' %labell[j])
    plt.plot(range(250,stop)*unit['time'],frac,c=C[j],linewidth=lw[j],label=labell[j])
    #plt.plot([],[],c=C[j],label=labell[j],lw=3)
plt.ylim(0,1)
plt.xlim(250*unit['time'].value,499*unit['time'].value)
plt.tick_params(labelsize=12,direction='in')
plt.ylabel('Fraction',fontsize=13)
plt.xlabel('Time [Myr]',fontsize=13)
plt.legend(loc=0,fontsize=11.5)
plt.savefig('D:/yeongu/plots/paperplot/new/surf_frac_all_allphase.png',dpi=300)# % labell[j])
plt.savefig('D:/yeongu/plots/paperplot/new/surf_frac_all_allphase.eps',format='eps',dpi=300)# % labell[j])
#plt.close()
plt.show()

'''
    for tidx in range(250, stop):  # time step 251, 331, 411, 501
        # plt.figure(figsize=(6, 5))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        #print(ds.field_list)

        # print(ds.derived_field_list)

        # full domain information
        # print ds.domain
        # rs = ds.read_all_data('reynold_stress')
        # print rs
        # information of grid #0
        # print ds.grids[0]


        T1 = ds.read_all_data('T1');
        coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)
        # print temp
        if j != 0:
            scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1
            s_cut = 0.5

        d1 = ds.read_all_data('density') * unit['density'].value  # ISM density # solar mass per pc3
        #d2 = copy.copy(d1)  # ds.read_all_data('density')*unit['density'].value  # ICM density
        d_mks_ism = d1#*8*8*8 # Mass , unit : Solar mass
        #d_mks_icm = d2#*8*8*8 # Mass , unit : Solar mass
        if j != 0:
            d_mks_ism[scalar > s_cut] = 0  # select only ISM
            #d_mks_icm[scalar < s_cut] = 0  # select only ICM
        
        d_c = copy.copy(d_mks_ism)
        d_c[temp > 20000] = 0  # cold ISM
        d_mks_ism[temp < 20000] = 0;
        d_h = d_mks_ism  # hot ISM
        
        #surf_c = (np.sum(d_c*8,axis=0))
        if j==3 or j==5:
            surf_c = (np.sum(d_mks_ism * 4, axis=0))
        else:
            surf_c = (np.sum(d_mks_ism * 8, axis=0))
        surf_c = surf_c.ravel()
        #print len(surf_c[surf_c < crit[j]])
        #print float(len(surf_c))
        lower = len(surf_c[surf_c < crit[j]])/float(len(surf_c))
        frac.append(lower)
        print tidx,lower

        #print range(250,stop)*unit['time']
    fraction = np.array(frac)
    np.savetxt('fraction_%s.txt' % labell[j], fraction)
    plt.plot(range(250,stop)*unit['time'],frac,c=C[j],linewidth=lw[j],label=labell[j])
plt.ylim(0,1)
plt.xlim(250*unit['time'].value,499*unit['time'].value)
plt.tick_params(labelsize=12,direction='in')
plt.ylabel('Fraction',fontsize=13)
plt.xlabel('Time [Myr]',fontsize=13)
plt.legend(loc=0,fontsize=11.5)
#plt.savefig('D:/yeongu/plots/paperplot/new/surf_frac_all_allphase.png')# % labell[j])
#plt.close()
plt.show()
    #plt.show()
        #plt.hist((surf_c),bins=np.arange(-3,3,0.05))
        #plt.xlim(-3,3)
        #plt.show()
        #surf_c = d_c
        #plt.imshow(surf_c,norm=LogNorm(),origin='lower')
        #plt.colorbar()
        #plt.clim(1,1e2)

        #plt.show()
'''
'''
    ##################### Plot #############################
    zz = range(895)
    plt.subplot(10, 1, 2 * k - 1)
    plt.plot(zz,Medi_e, 'k-', label='Early_%s' % labell[j])
    plt.plot(zz,np.nanmedian(Ext_e, axis=0),'r--',label='Ext')
    plt.plot(zz, np.nanmedian(Sg_e, axis=0), 'g--', label='Sg')
    plt.fill_between(zz, Min_e, Max_e, facecolor=C[j], alpha=0.5)

    plt.xticks([73, 198, 323, 448, 573, 698, 823], [])
    #plt.ylabel(r'Pressure $[K cm^{-3}]$')
    plt.ylabel(r'Weight')
    plt.yscale('log', nonposy='clip')
    plt.ylim(1e2, 1e5)
    plt.xlim(0, 896)
    plt.tick_params(which='major', direction='in')
    plt.tick_params(which='minor', direction='in')
    plt.legend(loc='upper right')

    plt.subplot(10, 1, 2 * k)
    plt.plot(zz,Medi_l, 'k-', label='Late_%s' % labell[j])
    plt.plot(zz,np.nanmedian(Ext_l, axis=0),'r--',label='Ext')
    plt.plot(zz, np.nanmedian(Sg_l, axis=0), 'g--', label='Sg')
    plt.fill_between(zz,Min_l, Max_l, facecolor=C[j], alpha=0.5)

    plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    #plt.ylabel(r'Pressure $[K cm^{-3}]$')
    plt.ylabel(r'Weight')
    if k == 5:
        plt.xlabel(r'z [kpc]')
    plt.yscale('log', nonposy='clip')
    plt.ylim(1e2, 1e5)
    plt.xlim(0, 896)
    plt.tick_params(which='major', direction='in')
    plt.tick_params(which='minor', direction='in')
    plt.legend(loc='upper right')

    k = k + 1

plt.tight_layout()
plt.savefig('D:/yeongu/plots/paperplot/new/Weight_5all_3phase_.png')
plt.show()
    '''
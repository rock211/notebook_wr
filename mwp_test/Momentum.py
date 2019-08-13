import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
sys.path.insert(0, 'C:/Users/woorak choi/Desktop/rps-tigress/notebook')
import pyathena as pa
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
import copy
from mpl_toolkits import axes_grid1
import time
from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl

# In[8]:


unit = pa.set_units(muH=1.4271)
print(unit)
# print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
# print unit['density']
kb = 1.38064852 * 1e-16  # boltzmann constant erg/K / erg = g cm2/s2
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
# z = np.arange(-3588, 3588, 8)

# meter = 3.24078*1e-17 # pc
# kg = 5.02785*1e-31 # solar mass

simid_t = (
'RPS_8pc_noICM_newacc', 'RPS_8pc_ICM0_newacc', 'RPS_8pc_ICM1_newacc', 'RPS_4pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc',
'RPS_4pc_ICM2_newacc', 'RPS_8pc_ICM3_newacc')
# labell = ('No ICM','Very Weak' ,'Weak', 'Strong', 'Very Strong','ICM1', 'ICM2', 'ICM3', 'ICM4') #'NonICM',
labell = ('No ICM', 'P1', 'P3', 'P3h', 'P7', 'P7h', 'P14', 'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
Model = [0, 8.63 * 1e3, 3.46 * 1e4, 3.46 * 1e4, 6.92 * 1e4, 6.92 * 1e4, 1.38 * 1e5]
Modell = [r'8.63*1e3', r'$3.46x10^4$', r'$6.92x10^4$', r'1.38*1e5']
S = ('-.', '--', '-')
C = (
'k', 'salmon', 'mediumblue', 'deepskyblue', 'darkgreen', 'lime', 'magenta', 'darkmagenta', 'goldenrod', 'royalblue',
'crimson')  # 'plum','orchid','purple'
# overplot Starformation rate of three different simulations

Myr = unit['time'].to('Myr').value
Msun = unit['mass'].to('Msun').value
agebin = 10
crit = 3
coolftn = pa.coolftn()
s_cut = 0.5

for j in (4,6):
    basedir = 'G:/yeongu/'
    if j == 5 or j == 6:
        stop = 473
    else:
        stop = 499

    simid = simid_t[j]

    ism_mom_c=[]
    ism_mom_h = []
    icm_mom=[]

    for tidx in range(250, stop):  # time step 251, 331, 411, 501
        # plt.figure(figsize=(6, 5))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # print(ds.field_list)

        # print(ds.derived_field_list)

        ############################## general properties ####################################
        T1 = ds.read_all_data('T1')
        temp = coolftn.get_temp(T1)

        if j != 0:
            scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1

        if j == 3 or j == 5:
            Area = 1024 * 1024 / 4 / 4.
            dz = 4
        else:
            Area = 1024 * 1024 / 8 / 8.  # np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
            dz = 8
            # Area = np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
        #######################################################################################

        d1 = ds.read_all_data('density') * unit['density'].value  # ISM density # solarmass / pc3

        d2 = copy.copy(d1)  # ds.read_all_data('density')*unit['density'].value  # ICM density
        d_mks_ism = d1  # *8*8*8 # Mass , unit : Solar mass
        d_mks_icm = d2  # *8*8*8 # Mass , unit : Solar mass
        if j != 0:
            d_mks_ism[scalar > s_cut] = 0  # select only ISM
            d_mks_icm[scalar < s_cut] = 0  # select only ICM
        d_c = copy.copy(d_mks_ism)
        d_c[temp > 20000] = 0  # cold ISM
        d_mks_ism[temp < 20000]=0; d_h = d_mks_ism

        vel = ds.read_all_data('velocity') ; vel2 = copy.copy(vel) #ds.read_all_data('velocity') # km/s
        vel_z = vel[:,:,:,2];  vel_z2 = vel2[:,:,:,2]
        if j != 0:
            vel_z[scalar > s_cut] = 0 # select only ISM
            vel_z2[scalar < s_cut] = 0 # select only ICM

        mom_ism_c = (np.sum(np.sum(d_c*vel_z*km, axis=1), axis=1) / Area)# Solar mass/pc3 * pc/s and horizontal average
        mom_ism_h = (np.sum(np.sum(d_h * vel_z * km, axis=1),axis=1) / Area)  # Solar mass/pc3 * pc/s and horizontal average

        if j !=0:
            mom_icm = (np.sum(np.sum(d_mks_icm*vel_z2*km, axis=1), axis=1) / Area)

        #############################################################################################
        ism_mom_c.append(mom_ism_c)
        ism_mom_h.append(mom_ism_h)

        if j!=0:
            icm_mom.append(mom_icm)

        print j, tidx

    Mom_ism_c = np.array(ism_mom_c)
    Mom_ism_h = np.array(ism_mom_h)
    Mom_icm = np.array(icm_mom)
    np.savetxt('Mom_ism_c_%s.txt' % labell[j], Mom_ism_c)
    np.savetxt('Mom_ism_h_%s.txt' % labell[j], Mom_ism_h)
    np.savetxt('Mom_icm_%s.txt' % labell[j], Mom_icm)


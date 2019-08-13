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
import time
from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl

# In[8]:
#time.sleep(7200)
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
#z = np.arange(-3588, 3588, 8)
#gext = (2. * np.pi * g * sig_star * (z) / ((z) ** 2 + z_s ** 2) ** (0.5) + 4. * np.pi * g * rho_dm * ((z) / (1 + (z) ** 2 / r_0 ** 2)))  # pc/s2
#gext = gext / km / km  # km/s2

# plt.plot((z-3584.)/1000,g_z)
# plt.show()
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

for j in (3,8):
    basedir = 'G:/yeongu/'
    if j == 5 or j == 6:
        stop = 473
    else:
        stop = 499

    simid = simid_t[j]
    Mom_up = []
    P_c = []
    P_h = []
    Weight_c = []
    Weight_h = []
    SFR = []
    Sg = []
    Ext = []
    D_c = []
    D_h = []
    D_icm = []
    Tot_icm = []
    Icm_turb = []
    Icm_ther = []
    ism_c_turb = []
    ism_c_ther = []
    ism_c_mag = []

    ism_h_turb = []
    ism_h_ther = []
    ism_h_mag = []

    if j == 3 or j == 5:
        z = np.arange(-3588, 3588, 4)
        gext = (2. * np.pi * g * sig_star * (z) / ((z) ** 2 + z_s ** 2) ** (0.5) + 4. * np.pi * g * rho_dm * (
                    (z) / (1 + (z) ** 2 / r_0 ** 2)))  # pc/s2
        gext = gext / km / km  # km/s2
    else:
        z = np.arange(-3588, 3588, 8)
        gext = (2. * np.pi * g * sig_star * (z) / ((z) ** 2 + z_s ** 2) ** (0.5) + 4. * np.pi * g * rho_dm * (
                (z) / (1 + (z) ** 2 / r_0 ** 2)))  # pc/s2
        gext = gext / km / km  # km/s2

    # plt.figure(figsize=(6, 5))
    for tidx in range(250, stop):  # time step 251, 331, 411, 501
        # plt.figure(figsize=(6, 5))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # print(ds.field_list)

        # print(ds.derived_field_list)
        '''
        ########################## star formation rate #######################################
        starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)

        star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
        star_clu2 = star_clu[star_clu['age'] * Myr < agebin]

        M_star = sum(star_clu2['mass']) * Msun
        sfr = M_star * unit['time'].value / (1e+6 * agebin * 1.024 * 1.024)
        ######################################################################################
        '''
        ############################## general properties ####################################
        T1 = ds.read_all_data('T1')
        temp = coolftn.get_temp(T1)

        if j != 0:
            scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1

        if j == 3 or j == 5:
            Area = 1024 * 1024 / 4 / 4.
            dz = 4
        else:
            Area = 1024 * 1024 / 8 / 8.
            dz = 8
            # Area = np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
        #######################################################################################


        #################################### pressure calculation ###################################

        B = ds.read_all_data('magnetic_field')
        Pmag1 = B[:,:,:,0]**2/2; Pmag2 = B[:,:,:,1]**2/2; Pmag3 = B[:,:,:,2]**2/2
        mag_p = Pmag1+Pmag2-Pmag3
        mag_p *= (unit['pressure']/kb).cgs.value
        mag_p2 = copy.copy(mag_p)
        if j != 0:
            mag_p[scalar > s_cut] = 0 # ISM
        mag_p_c = copy.copy(mag_p)
        mag_p[temp < 20000] = 0 ; mag_p_h = mag_p
        mag_p_c[temp > 20000] = 0
        #if j != 0:
        #    mag_p2[scalar < s_cut] = 0 # ICM

        d1 = ds.read_all_data('density')*unit['density'].value # ISM density
        d2 = copy.copy(d1)#ds.read_all_data('density')*unit['density'].value  # ICM density
        d_mks_ism = d1*6.76598*1e-11
        d_mks_icm = d2*6.76598*1e-11 # solarmass/pc3 to kg/km3
        if j != 0:
            d_mks_ism[scalar > s_cut] = 0 # select only ISM
            d_mks_icm[scalar < s_cut] = 0 # select only ICM
        #if j != 0:
        #    Area_icm = 1024*1024/8/8.
            #for i in range(896):
            #    Area_icm.append(np.count_nonzero(d_mks_icm[i, :, :]))
            #print Area_icm

        pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # ISM thermal pressure
        pre2 = copy.copy(pre)

        if j != 0:
            pre[scalar > s_cut] = 0 #select only ISM
            pre2[scalar < s_cut] = 0  # select only ICM

        pre_c = copy.copy(pre) # ISM for cold
        pre[temp < 20000] = 0 ; pre_h = pre # hot ISM
        pre_c[temp > 20000] = 0 # cold ISM

        vel = ds.read_all_data('velocity') ; vel2 = copy.copy(vel) #ds.read_all_data('velocity') # km/s
        vel_z = vel[:,:,:,2];  vel_z2 = vel2[:,:,:,2]
        if j != 0:
            vel_z[scalar > s_cut] = 0 # select only ISM
            vel_z2[scalar < s_cut] = 0 # select only ICM

        turb = d_mks_ism*(vel_z**2)*1e-2/kb # ISM turbulence pressure / 1e-2 is conversion factor from kg/km to g/cm
        turb_c = copy.copy(turb)

        turb_c[temp > 20000] = 0 # cold ISM
        turb[temp < 20000] = 0 ; turb_h=turb # hot ISM

        magp_c = (np.sum(np.sum(mag_p_c, axis=1), axis=1)/Area)
        magp_h = (np.sum(np.sum(mag_p_h, axis=1), axis=1)/Area)
        turb_c = (np.sum(np.sum(turb_c,axis=1),axis=1)/Area) # ISM turbulence pressure / cold
        turb_h = (np.sum(np.sum(turb_h,axis=1),axis=1)/Area) # ISM hot
        pre_c = (np.sum(np.sum(pre_c,axis=1),axis=1)/Area) # ISM thermal pressure / cold
        pre_h = (np.sum(np.sum(pre_h,axis=1),axis=1)/Area) # ISM hot

        tot_p_c = turb_c + pre_c + magp_c # Unit : K cm-3
        tot_p_h = turb_h + pre_h + magp_h # Unit : K cm-3
        #plt.semilogy(tot_p)
        #print tot_p
        #print np.max(tot_p)
        #plt.ylim(1e2,1e5)
        #plt.show()
        '''
        ############ ICM pressure ###########
        if j!=0 :
            ICM = d_mks_icm*(vel_z2**2)*1e-2/kb # Ram Pressure : density x velocity^2, 1e-2 : convert to cgs
            ICM = np.sum(np.sum(ICM,axis=1),axis=1)/Area # y projection
            pre2 = np.sum(np.sum(pre2,axis=1),axis=1)/Area # ICM thermal pressure
            tot_icm = pre2 + ICM
        #####################################
        '''
        #############################################################################################


        # SFR.append(sfr)
        P_c.append(tot_p_c[0:-1])
        P_h.append(tot_p_h[0:-1])
        # Icm_ther.append(pre2)
        # Icm_turb.append(ICM)
        # Tot_icm.append(tot_icm)
        ism_c_ther.append(pre_c[0:-1])
        ism_c_turb.append(turb_c[0:-1])
        ism_c_mag.append(magp_c[0:-1])
        ism_h_ther.append(pre_h[0:-1])
        ism_h_turb.append(turb_h[0:-1])
        ism_h_mag.append(magp_h[0:-1])
        # D_icm.append(dicm)
        print j, tidx

    #SFR = np.array(SFR)

    # Sg = np.array(Sg)
    # Ext = np.array(Ext)

    # Icm_ther = np.array(Icm_ther)
    # Icm_turb = np.array(Icm_turb)
    # Tot_icm = np.array(Tot_icm)

    P_c = np.array(P_c)
    Pc_ther = np.array(ism_c_ther)
    Pc_turb = np.array(ism_c_turb)
    Pc_mag = np.array(ism_c_mag)
    P_h = np.array(P_h)
    Ph_ther = np.array(ism_h_ther)
    Ph_turb = np.array(ism_h_turb)
    Ph_mag = np.array(ism_h_mag)


    # np.savetxt('SFR_%s.txt' % labell[j],SFR)
    np.savetxt('Pc_%s.txt' % labell[j],P_c)
    np.savetxt('Pc_ther_%s.txt' % labell[j], Pc_ther)
    np.savetxt('Pc_turb_%s.txt' % labell[j], Pc_turb)
    np.savetxt('Pc_mag_%s.txt' % labell[j], Pc_mag)

    np.savetxt('Ph_%s.txt' % labell[j],P_h)
    np.savetxt('Ph_ther_%s.txt' % labell[j], Ph_ther)
    np.savetxt('Ph_turb_%s.txt' % labell[j], Ph_turb)
    np.savetxt('Ph_mag_%s.txt' % labell[j], Ph_mag)

    # np.savetxt('icm_ther_%s.txt' % labell[j],Icm_ther)
    # np.savetxt('icm_turb_%s.txt' % labell[j],Icm_turb)
    # np.savetxt('icm_tot_%s.txt' % labell[j],Tot_icm)

    # print Weight.shape
    # print SFR.shape
    # print P.shape
    # SFR_e=np.array(SFR[:crit])
    # SFR_l=np.array(SFR[crit:])
    # P_e=np.array(P[:crit])
    # P_l=np.array(P[crit:])
    # print P_l

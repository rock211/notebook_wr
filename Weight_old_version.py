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

#time.sleep(3600)

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

for j in (2,8):
    basedir = 'G:/yeongu/'
    if j == 5 or j == 6:
        stop = 473
    else:
        stop = 499

    simid = simid_t[j]
    Mom_up = []
    P = []
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
    ism_turb = []
    ism_ther = []
    ism_mag = []

    if j == 3 or j == 5:
        z = np.arange(-3588, 3588, 4)
        gext = (2. * np.pi * g * sig_star * (z) / ((z) ** 2 + z_s ** 2) ** (0.5) + 4. * np.pi * g * rho_dm * ((z) / (1 + (z) ** 2 / r_0 ** 2)))  # pc/s2

    else:
        z = np.arange(-3588, 3588, 8)
        gext = (2. * np.pi * g * sig_star * (z) / ((z) ** 2 + z_s ** 2) ** (0.5) + 4. * np.pi * g * rho_dm * ((z) / (1 + (z) ** 2 / r_0 ** 2)))  # pc/s2


    for tidx in range(250, stop):  # time step 251, 331, 411, 501

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        print(ds.field_list)

        print(ds.derived_field_list)
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
            Area = 1024 * 1024 / 8 / 8.  # np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
            dz = 8
            # Area = np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
        #######################################################################################

        ############################## Weight calculation ########################################

        pot = ds.read_all_data('gravitational_potential')
        pot2 = copy.copy(pot)
        if j != 0:
            pot[scalar > s_cut] = 0  # ISM
            pot2[scalar < s_cut] = 0 # ICM

        pot_c = copy.copy(pot)
        pot_c[temp > 20000] = 0  # cold ISM
        pot[temp < 20000] = 0;
        pot_h = pot # hot ISM
        #print np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
        #print np.sum(np.count_nonzero(pot_c,axis=2),axis=1)
        Area_c = np.sum(np.count_nonzero(pot_c, axis=1), axis=1)
        Area_h = np.sum(np.count_nonzero(pot_h, axis=1), axis=1)

        #print Area_c,Area_h
        #pot_icm = np.sum(np.sum(pot2, axis=1), axis=1)
        pot_c = np.sum(np.sum(pot_c,axis=1),axis=1)
        pot_h = np.sum(np.sum(pot_h, axis=1), axis=1)

        #gsg_c = np.diff(pot_c) / dz/Area# km/s2
        #gsg_h = np.diff(pot_h) / dz/Area  # km/s2
        gsg_c = np.gradient(pot_c,dz)/Area
        gsg_h = np.gradient(pot_h,dz)/Area
        #gsg_icm = np.gradient(pot_icm,dz)/Area

        #plt.plot(gsg2)
        #plt.plot(gsg,'k--')
        #plt.plot(gext[0:-1]*Area2/(128*128.))
        #plt.show()

        if j==3 or j==5:
            gext_c = gext[1:-1] * Area_c / Area
            gext_h = gext[1:-1] * Area_h / Area
        else:
            gext_c=gext[0:-1]*Area_c / Area
            gext_h = gext[0:-1] * Area_h / Area

        d1 = ds.read_all_data('density') * unit['density'].value  # ISM density # solarmass / pc3
        d2 = copy.copy(d1)  # ds.read_all_data('density')*unit['density'].value  # ICM density
        d_mks_ism = d1  # *8*8*8 # Mass , unit : Solar mass
        d_mks_icm = d2  # *8*8*8 # Mass , unit : Solar mass
        if j != 0:
            d_mks_ism[scalar > s_cut] = 0  # select only ISM
            d_mks_icm[scalar < s_cut] = 0  # select only ICM
        d_c = copy.copy(d_mks_ism)
        d_c[temp > 20000] = 0  # cold ISM
        d_mks_ism[temp < 20000] = 0;
        d_h = d_mks_ism  # hot ISM

        # d_h1 =  np.nan_to_num(np.sum(np.sum(d_h, axis=1), axis=1)/Area)
        # d_h2 = np.nan_to_num(np.mean(np.mean(d_h, axis=1), axis=1))
        # print d_h1.max(), d_h2.max()

        d_c = np.nan_to_num(np.sum(np.sum(d_c, axis=1), axis=1) / Area)
        d_h = np.nan_to_num(np.sum(np.sum(d_h, axis=1), axis=1) / Area)
        # dicm = np.nan_to_num(np.sum(np.sum(d_mks_icm, axis=1), axis=1)/Area)

        #plt.plot(gext_c,'r-')
        #plt.plot(gext_h,'b-')
        plt.plot(gext_c+gext_h,'k--')
        #plt.plot(gsg_c/(3.086e13)**2,'m--')
        #plt.plot(gsg_h / (3.086e13) ** 2, 'c--')
        plt.plot(gsg_c/(3.086e13)**2+gsg_h / (3.086e13) ** 2, 'g--')
        plt.show()

        # Weight calculation
        # integral (density x gravity) from zmax to z
        # external gravity could be derived from exact equation
        # self gravity could be derived from potential differentiation

        #### Warm/cold ####
        ext_c = (gext_c*d_c/kb)[::-1].cumsum()[::-1]*dz*cm/gram
        sg_c = (np.nan_to_num(gsg_c)/(3.086e13)**2*d_c/kb)[::-1].cumsum()[::-1]*dz*cm/gram
        tot_w_c = sg_c + ext_c # Unit K cm-3

        ### Hot ###
        ext_h = (gext_h*d_h/kb)[::-1].cumsum()[::-1]*dz*cm/gram
        sg_h = (np.nan_to_num(gsg_h)/(3.086e13)**2*d_h/kb)[::-1].cumsum()[::-1]*dz*cm/gram
        tot_w_h = sg_h+ext_h # Unit K cm-3

        ### ICM weight ###
        # ext_icm = (gext/(3.086e13)**2*dicm/kb)[::-1].cumsum()[::-1]*dz*cm/gram
        # sg_icm = (np.nan_to_num(gsg_icm)/(3.086e13)**2*dicm/kb)[::-1].cumsum()[::-1]*dz*cm/gram
        # weight_icm = ext_icm+sg_icm

        #plt.semilogy(dicm,'k--')
        plt.plot(ext_c,'r-')
        plt.plot(sg_c,'g--')
        plt.plot(ext_h,'m-')
        plt.plot(sg_h,'c--')
        #plt.plot(weight_icm,'b-')
        plt.plot(tot_w_c,'k')
        plt.plot(tot_w_h, 'k--')
        plt.show()
        #############################################################################################

        Weight_c.append(tot_w_c)
        Weight_h.append(tot_w_h)
        # Sg.append(sg)
        # Ext.append(ext)

        # SFR.append(sfr)

        #D_c.append(d_c)
        #D_h.append(d_h)
        #D_icm.append(dicm)
        print j, tidx

    Weight_c = np.array(Weight_c)
    Weight_h = np.array(Weight_h)
    #Sg = np.array(Sg)
    #Ext = np.array(Ext)

    #SFR = np.array(SFR)
    #Dc = np.array(D_c)
    #Dh = np.array(D_h)
    #D_icm = np.array(D_icm)


    np.savetxt('Weight_c_%s.txt' % labell[j],Weight_c)
    np.savetxt('Weight_h_%s.txt' % labell[j], Weight_h)
    #np.savetxt('Sg_%s.txt' % labell[j],Sg)
    #np.savetxt('Ext_%s.txt' % labell[j],Ext)

    #np.savetxt('SFR_%s.txt' % labell[j],SFR)

    #np.savetxt('Dc_%s.txt' % labell[j], Dc)
    #np.savetxt('Dh_%s.txt' % labell[j],Dh)
    #np.savetxt('Dicm_%s.txt' % labell[j],D_icm)

    # print Weight.shape
    # print SFR.shape
    # print P.shape
    # SFR_e=np.array(SFR[:crit])
    # SFR_l=np.array(SFR[crit:])
    # P_e=np.array(P[:crit])
    # P_l=np.array(P[crit:])
    # print P_l

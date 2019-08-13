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

gext = (2. * np.pi * g * sig_star * (z) / ((z) ** 2 + z_s ** 2) ** (0.5) + 4. * np.pi * g * rho_dm * (
        (z) / (1 + (z) ** 2 / r_0 ** 2)))  # pc/s2
gext = gext / km / km  # km/s2

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
for j in (0, 1, 2, 4, 6):
    basedir = 'G:/yeongu/'
    if j == 5 or j == 6:
        stop = 473
    else:
        stop = 499

    simid = simid_t[j]
    Mom_up = []
    P = []
    Weight = []
    SFR = []
    Sg = []
    Ext = []
    D_c = []
    D_h = []
    D_icm = []
    Tot_icm = []
    Icm_turb = []
    Icm_ther = []
    M=[]
    # plt.figure(figsize=(6, 5))
    for tidx in range(250, stop):  # time step 251, 331, 411, 501
        # plt.figure(figsize=(6, 5))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # print(ds.field_list)

        # print(ds.derived_field_list)
        ########################## star formation rate #######################################
        starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)

        star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
        M_star = sum(star_clu['mass']) * Msun/(1024*1024)
        #print M_star/(1024*1024)
        M.append(M_star)


        '''
        star_clu2 = star_clu[star_clu['age'] * Myr < agebin]

        M_star = sum(star_clu2['mass']) * Msun
        sfr = M_star * unit['time'].value / (1e+6 * agebin * 1.024 * 1.024)
        ######################################################################################

        ############################## general properties ####################################
        T1 = ds.read_all_data('T1');
        coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)

        if j != 0:
            scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1
            s_cut = 0.5
        if j == 3 or j == 5:
            Area = 1024 * 1024 / 4 / 4.
        else:
            Area = 1024 * 1024 / 8 / 8.  # np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
            # Area = np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
        #######################################################################################

        ############################## Weight calculation ########################################

        pot = ds.read_all_data('gravitational_potential')
        pot2 = copy.copy(pot)
        if j != 0:
            pot[scalar > s_cut] = 0  # ISM
            pot2[scalar < s_cut] = 0  # ICM

        pot_c = copy.copy(pot)
        pot_c[temp > 20000] = 0  # cold ISM
        pot[temp < 20000] = 0;
        pot_h = pot  # hot ISM
        # print np.sum(np.count_nonzero(pot_c,axis=1),axis=1)
        # print np.sum(np.count_nonzero(pot_c,axis=2),axis=1)
        Area2 = np.sum(np.count_nonzero(pot_c, axis=1), axis=1)

        # print Area2
        # pot_icm = np.sum(np.sum(pot2, axis=1), axis=1) / Area
        pot_c = np.sum(np.sum(pot_c, axis=1), axis=1) / Area2

        if j == 3 or j == 5:
            dz = 4
        else:
            dz = 8
        # gsg2 = np.gradient(pot_c,dz)
        gsg = np.diff(pot_c) / dz * Area2[1::] / Area  # km/s2
        # plt.plot(gsg2)
        # plt.plot(gsg,'k--')
        # plt.plot(gext[0:-1]*Area2/(128*128.))

        # plt.show()
        # gsg_icm = np.diff(pot_icm) / dz
        gext2 = gext[0:-1] * Area2 / (128 * 128.)
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
        # dicm = np.nan_to_num(np.sum(np.sum(d_mks_icm, axis=1), axis=1)/Area)
        # plt.plot(d_c)
        # plt.show()
        # Weight calculation
        # integral (density x gravity) from zmax to z
        # external gravity could be derived from exact equation
        # self gravity could be derived from potential differentiation

        ext = (gext2[0:-1] / (3.086e13) ** 2 * d_c[0:-1] / kb)[::-1].cumsum()[::-1] * cm * 8 / gram
        sg = (np.nan_to_num(gsg[0::]) / (3.086e13) ** 2 * d_c[0:-1] / kb)[::-1].cumsum()[::-1] * cm * 8 / gram

        # plt.plot(ext)
        # plt.plot(sg)
        # plt.show()

        # plt.plot(ext)
        # plt.plot(sg)
        # plt.show()
        tot_w = sg + ext

        # ext_icm = (gext[0:-2]/(3.086e13)**2*dicm[0:-1]/kb)[::-1].cumsum()[::-1]*cm*8/gram
        # sg_icm = (np.nan_to_num(gsg_icm[0::])/(3.086e13)**2*dicm[0:-1]/kb)[::-1].cumsum()[::-1]*cm*8/gram
        # weight_icm = ext_icm+sg_icm

        # plt.semilogy(dicm,'k--')
        # plt.plot(ext_icm,'r-')
        # plt.plot(sg_icm,'g--')
        # plt.plot(weight_icm,'b-')
        # plt.plot(tot_w,'k')
        # plt.show()
        #############################################################################################
        '''
        '''
        if j==3 or j==5:
            Area = 1024*1024/4/4.
        else:
            Area = 1024 * 1024 / 8 / 8.
        #################################### pressure calculation ###################################

        B = ds.read_all_data('magnetic_field')
        Pmag1 = B[:,:,:,0]**2/2
        Pmag2 = B[:,:,:,1]**2/2
        Pmag3 = B[:,:,:,2]**2/2
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

        #print pre2
        if j != 0:
            pre[scalar > s_cut] = 0 #select only ISM
            pre2[scalar < s_cut] = 0  # select only ICM
        #print pre2
        pre_c = copy.copy(pre) # ISM for cold
        pre[temp < 20000] = 0 ; pre_h = pre # hot ISM
        pre_c[temp > 20000] = 0 # cold ISM


        vel = ds.read_all_data('velocity') ; vel2 = copy.copy(vel) #ds.read_all_data('velocity') # km/s
        vel_z = vel[:,:,:,2];  vel_z2 = vel2[:,:,:,2]
        if j != 0:
            vel_z[scalar > s_cut] = 0 # select only ISM
            vel_z2[scalar < s_cut] = 0 # select only ICM


        turb = d_mks_ism*(vel_z**2)*1e-2/kb # ISM turbulence pressure 1e-2 is conversion factor from mks to g cm
        turb_c = copy.copy(turb)

        turb_c[temp > 20000] = 0 # cold ISM
        turb[temp < 20000] = 0 ; turb_h=turb # hot ISM

        magp_c = (np.sum(np.sum(mag_p_c, axis=1), axis=1)/Area)
        magp_h = (np.sum(np.sum(mag_p_h, axis=1), axis=1)/Area)
        turb_c = (np.sum(np.sum(turb_c,axis=1),axis=1)/Area) # ISM turbulence pressure / cold
        turb_h = (np.sum(np.sum(turb_h,axis=1),axis=1)/Area) # ISM hot
        pre_c = (np.sum(np.sum(pre_c,axis=1),axis=1)/Area) # ISM thermal pressure / cold
        pre_h = (np.sum(np.sum(pre_h,axis=1),axis=1)/Area) # ISM hot
        #m_pre_c = np.repeat(np.mean(np.mean(m_pre_c, axis=1), axis=1), 8) # magnetic pressure
        #m_pre_h = np.repeat(np.mean(np.mean(m_pre_h, axis=1), axis=1), 8) # magnetic pressure
        tot_p = turb_c + pre_c + magp_c# + turb_h + pre_h + magp_h
        plt.semilogy(tot_p)
        print tot_p
        print np.max(tot_p)
        plt.ylim(1e2,1e5)
        plt.show()
        if j!=0 :
            ICM = d_mks_icm*(vel_z2**2)*1e-2/kb # Ram Pressure : density x velocity^2, 1e-2 : convert to cgs
            ICM = np.sum(np.sum(ICM,axis=1),axis=1)/Area # y projection
            pre2 = np.sum(np.sum(pre2,axis=1),axis=1)/Area # ICM thermal pressure
            tot_icm = pre2 + ICM

        #plt.semilogy(tot_icm,'r')
        #plt.semilogy(pre2,'g')
        #plt.semilogy(ICM,'b')
        #plt.axhline(Model[j],ls='--',c='k')
        #plt.semilogy(tot_p,'k')
        #plt.ylim(1e2,1e5)
        #plt.show()
        '''

        #############################################################################################
        #Sg.append(sg)
        #Ext.append(ext)
        # D_c.append(d_c)
        # D_h.append(d_h)
        #Weight.append(tot_w)
        # SFR.append(sfr)
        # P.append(tot_p[0:-1])
        # Icm_ther.append(pre2)
        # Icm_turb.append(ICM)
        # Tot_icm.append(tot_icm)
        # D_icm.append(dicm)
        print j, tidx
    plt.plot(M,c=C[j],label=labell[j])
    plt.legend(loc=0)

    #Weight = np.array(Weight)
    # SFR = np.array(SFR)
    # P = np.array(P)
    #Sg = np.array(Sg)
    #Ext = np.array(Ext)
    # D_c = np.array(D_c)
    # D_h = np.array(D_h)
    # Icm_ther = np.array(Icm_ther)
    # Icm_turb = np.array(Icm_turb)
    # Tot_icm = np.array(Tot_icm)
    # D_icm = np.array(D_icm)

    #np.savetxt('Weight_%s.txt' % labell[j], Weight)
    # np.savetxt('SFR_%s.txt' % labell[j],SFR)
    # np.savetxt('P_%s.txt' % labell[j],P)
    #np.savetxt('Sg_%s.txt' % labell[j], Sg)
    #np.savetxt('Ext_%s.txt' % labell[j], Ext)
    # np.savetxt('Dc_%s.txt' % labell[j],D_c)
    # np.savetxt('Dh_%s.txt' % labell[j],D_h)
    # np.savetxt('icm_ther_%s.txt' % labell[j],Icm_ther)
    # np.savetxt('icm_turb_%s.txt' % labell[j],Icm_turb)
    # np.savetxt('icm_tot_%s.txt' % labell[j],Tot_icm)
    # np.savetxt('Dicm_%s.txt' % labell[j],D_icm)

    # print Weight.shape
    # print SFR.shape
    # print P.shape
    # SFR_e=np.array(SFR[:crit])
    # SFR_l=np.array(SFR[crit:])
    # P_e=np.array(P[:crit])
    # P_l=np.array(P[crit:])
    # print P_l
    #plt.show()
    '''
    if min(SFR_e)==0 or min(SFR_l)==0 or min(P_e)==0 or min(P_l)==0:
        sfr_e_min=0
        sfr_l_min=0
        p_e_min=0
        p_l_min=0
    else:
        sfr_e_min=np.log10(np.min(SFR_e))
        sfr_l_min=np.log10(np.min(SFR_l))
        p_e_min=np.log10(np.min(P_e))
        p_l_min=np.log10(np.min(P_l))

    ############################### Upper and Lower #########################################
    SFR_e_errors = [np.mean(SFR_e)-np.min(SFR_e),np.max(SFR_e)-np.mean(SFR_e)] # [np.nonzero(SFR_e)] nonzero selection
    print SFR_e_errors
    SFR_l_errors = [np.mean(SFR_l)-np.min(SFR_l),np.max(SFR_l)-np.mean(SFR_l)]
    print SFR_l_errors
    P_e_errors = [np.mean(P_e)-np.min(P_e),np.max(P_e)-np.mean(P_e)]
    print np.array([P_e_errors]).T
    P_l_errors = [np.mean(P_l)-np.min(P_l),np.max(P_l)-np.mean(P_l)]
    print np.array([P_l_errors]).T
    #########################################################################################

    plt.errorbar(np.mean(P_e),np.mean(SFR_e),yerr=np.array([SFR_e_errors]).T,xerr=np.array([P_e_errors]).T,ecolor=C[j],marker='o',mfc=C[j],mec=C[j],capsize=3.5,label=labell[j])
    plt.errorbar(np.mean(P_l),np.mean(SFR_l),yerr=np.array([SFR_l_errors]).T,xerr=np.array([P_l_errors]).T,ecolor=C[j],marker='s',mfc='none',mec=C[j],capsize=3.5)

##### legend #####
s1 = plt.scatter([], [], marker='o',edgecolors='k',facecolors='k', alpha=0.7,label='Early')
s2 = plt.scatter([], [], marker='s',edgecolors='k',facecolors='none', alpha=0.7,label='Late')
##################

##### theoretical line #####
pressure=np.arange(1e2,7e5)
sfr_l = 2.1*1e-3*(pressure/1e4)**1.18
plt.plot((pressure),(sfr_l),ls='--',c='k')
############################

##### legend modify #####
h, l =plt.gca().get_legend_handles_labels()
legend1=plt.legend(h[2:],l[2:],loc='upper left')#,fontsize=14,framealpha=0.3)
plt.legend(h[:2],l[:2],loc='center left')#,fontsize=12)
plt.gca().add_artist(legend1)
#########################

plt.xscale('log', nonposx='clip')
plt.yscale('log', nonposy='clip')
#plt.xlim(3e2,3e5)
plt.ylim(10**(-4.5),10**(-1.5))
plt.ylabel('log(SFR)')
plt.xlabel('log(Weight)')
#plt.legend(loc=0)
plt.tight_layout()
#plt.savefig('D:/yeongu/plots/paperplot/new/SFR_P_wz.png',dpi=300)
#plt.savefig('D:/yeongu/plots/paperplot/new/SFR_P_wz.eps',format='eps',dpi=300)
plt.show()
'''
'''
    cc=range(250,stop)

    plt.scatter(np.log10(P),np.log10(SFR),s=5,c=cc)
    plt.colorbar()
'''
'''
    plt.ylabel('log(SFR)')
    plt.xlabel('log(P/kb)')
    plt.ylim(-4.5,-1.5)
    plt.xlim(2,5)
    plt.savefig('D:/yeongu/plots/%s_SFR_P.png' % simid_t[j])
    plt.show()
    plt.close()
'''
plt.show()
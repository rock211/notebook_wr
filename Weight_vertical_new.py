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
labell = ('No ICM', 'P1', 'P3', 'P3h', 'P7', 'P7h', 'P14', 'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
Model = [8.63 * 1e3, 3.46 * 1e4, 6.92 * 1e4, 1.38 * 1e5]
Modell = [r'8.63*1e3', r'$3.46x10^4$', r'$6.92x10^4$', r'1.38*1e5']
#C = ('k', 'salmon', 'mediumblue', 'deepskyblue', 'darkgreen', 'lime', 'firebrick', 'darkmagenta', 'goldenrod', 'royalblue','crimson')  # 'plum','orchid','purple'
C = ('gray', 'lightskyblue', 'dodgerblue','mediumblue' ,'goldenrod','salmon', 'firebrick')

S = ('-.', '--', '-')

# overplot Starformation rate of three different simulations
k = 1
plt.figure(figsize=(6, 13))
for j in (0, 1, 2, 4,6):
    basedir = 'G:/yeongu/'
    if j == 5 or j == 6:
        stop = 474
    else:
        stop = 499
    simid = simid_t[j]
    Mom_up = []
    Turb_e = []
    Turb_l = []
    Mag_e = []
    Mag_l = []
    Ther_e = []
    Ther_l = []
    Tot_e = []
    Tot_l = []
    Ext_e = []
    Ext_l = []
    Sg_e = []
    Sg_l = []
    #plt.figure(figsize=(6, 10))

    crit = 94
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

        #if j != 0:
        #    scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1
        #    s_cut = 0.5

        pot = ds.read_all_data('gravitational_potential')

        #if j != 0:
        #    pot[scalar > s_cut] = 0  # ISM

        pot_c = copy.copy(pot)
        pot_c[temp > 20000] = 0  # cold ISM
        pot[temp < 20000] = 0;
        pot_h = pot # hot ISM
        Area = 1024*1024/8/8.#np.sum(np.count_nonzero(pot_c,axis=1),axis=1)

        pot_c = np.sum(np.sum(pot_c,axis=1),axis=1)/Area
        dz = 8
        gsg = np.diff(pot_c) / dz # km/s2
        #print gsg

        #print gext.shape
        #plt.plot(gsg)
        #plt.plot(gext/km/km)
        #plt.show()

        d1 = ds.read_all_data('density') * unit['density'].value  # ISM density
        d2 = copy.copy(d1)  # ds.read_all_data('density')*unit['density'].value  # ICM density
        d_mks_ism = d1#*8*8*8 # Mass , unit : Solar mass
        d_mks_icm = d2#*8*8*8 # Mass , unit : Solar mass

        #if j != 0:
        #    d_mks_ism[scalar > s_cut] = 0  # select only ISM
        #    d_mks_icm[scalar < s_cut] = 0  # select only ICM

        d_c = copy.copy(d_mks_ism)
        d_c[temp > 20000] = 0  # cold ISM
        d_mks_ism[temp < 20000] = 0;
        d_h = d_mks_ism  # hot ISM
        d_c = np.nan_to_num(np.sum(np.sum(d_c, axis=1), axis=1)/Area)

        ext = (gext[0:-2]/(3.086e13)**2*d_c[0:-1]/kb)[::-1].cumsum()[::-1]*cm*8/gram
        sg = (np.nan_to_num(gsg[0::])/(3.086e13)**2*d_c[0:-1]/kb)[::-1].cumsum()[::-1]*cm*8/gram
        tot = ext+sg
        #print ext.shape
        #print sg.shape
        print np.nanmax(np.log10(tot))
        #plt.semilogy(ext,c='r')
        #plt.semilogy(sg,c='g')
        #plt.semilogy(tot,c='b')
        #print ext
        #print sg
        #plt.ylim(1e2,1e6)
        #plt.show()
        #d_h = np.repeat(np.sum(np.sum(d_h, axis=1), axis=1), 8)

        if tidx < crit + 250:
            Tot_e.append(tot)
            Sg_e.append(sg)
            Ext_e.append(ext)
            #Turb_e.append(turb_c)
            #Ther_e.append(pre_c)
            #Mag_e.append(magp_c)
        else:
            Tot_l.append(tot)
            Sg_l.append(sg)
            Ext_l.append(ext)
            #Turb_l.append(turb_c)
            #Ther_l.append(pre_c)
            #Mag_l.append(magp_c)
        print tidx

    Tot_e = np.array(Tot_e); Tot_l = np.array(Tot_l)
    Ext_e = np.array(Ext_e); Ext_l = np.array(Ext_l)
    Sg_e = np.array(Sg_e); Sg_l = np.array(Sg_l)
    #Turb_e = np.array(Turb_e);
    #Turb_l = np.array(Turb_l)
    #Ther_e = np.array(Ther_e);
    #Ther_l = np.array(Ther_l)
    #Mag_e = np.array(Mag_e);
    #Mag_l = np.array(Mag_l)

    Medi_e = np.nanmedian(Tot_e, axis=0);    Medi_l = np.nanmedian(Tot_l, axis=0)
    Min_e = np.nanpercentile(Tot_e,25, axis=0);    Min_l = np.nanpercentile(Tot_l,25, axis=0)
    Max_e = np.nanpercentile(Tot_e,75, axis=0);    Max_l = np.nanpercentile(Tot_l,75, axis=0)

    ##################### Plot #############################
    zz = range(895)
    plt.subplot(5, 1, k)
    plt.plot(zz, Medi_e, 'c-', label='Early_%s' % labell[j])
    #plt.plot(zz,np.nanmedian(Ext_e, axis=0),'r--',label='Ext')
    #plt.plot(zz, np.nanmedian(Sg_e, axis=0), 'g--', label='Sg')
    plt.fill_between(zz, Min_e, Max_e, facecolor=C[j], alpha=0.5)

    #plt.xticks([73, 198, 323, 448, 573, 698, 823], [])
    #plt.ylabel(r'Pressure $[K cm^{-3}]$')
    #plt.ylabel(r'Weight')
    #plt.yscale('log', nonposy='clip')
    #plt.ylim(1e2, 1e5)
    #plt.xlim(0, 896)
    #plt.tick_params(which='major', direction='in')
    #plt.tick_params(which='minor', direction='in')
    #plt.legend(loc='upper right')

    #plt.subplot(6, 1, 2 * k)
    plt.plot(zz,Medi_l, 'm--', label='Late_%s' % labell[j])
    #plt.plot(zz,np.nanmedian(Ext_l, axis=0),'r--',label='Ext')
    #plt.plot(zz, np.nanmedian(Sg_l, axis=0), 'g--', label='Sg')
    plt.fill_between(zz,Min_l, Max_l, facecolor=C[j], alpha=0.5)

    plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    #plt.ylabel(r'Pressure $[K cm^{-3}]$')
    #plt.ylabel(r'Weight')
    plt.ylabel('Density')
    if k == 5:
        plt.xlabel(r'z [kpc]')
    plt.yscale('log', nonposy='clip')
    plt.ylim(1e2, 1e5)
    #plt.ylim(1e-4, 1e-1)
    plt.xlim(0, 896)
    plt.tick_params(which='major', direction='in')
    plt.tick_params(which='minor', direction='in')
    plt.legend(loc='upper right')

    k = k + 1

plt.tight_layout()
#plt.savefig('D:/yeongu/plots/paperplot/new/densityw.png')
plt.savefig('D:/yeongu/plots/paperplot/new/Weight_wICM_3phase_new.png')
plt.show()

'''
        ######### restoring pressure #############
        d = np.repeat(d_mks_ism,8,axis=0)
        restP = d*g_z[:,None,None]
        restP = np.mean(np.mean(restP,axis=1),axis=1)*meter*10/kg/kb
        P = []
        for zz in range(7168) :
            P.append(np.sum(restP[zz:7168]))
        ##########################################

        #plt.semilogy(z, P, c='crimson' , label='Restoring P')
        plt.semilogy(z,magp_c,c='g',label='Cold/Warm magnetic pressure')
        plt.semilogy(z,pre_c,c='crimson',label='Cold/Warm thermal pressure')
        #plt.semilogy(z, pre_c, c='crimson', label='Cold/Warm thermal pressure')
        plt.semilogy(z,turb_c,c='salmon',label='Cold/Warm turbulence pressure')
        #plt.semilogy(z,tot_c,c='magenta',label='ISM total pressure')
        plt.semilogy(z, tot_c, c='magenta', label='ISM total pressure')
        if j != 0:
            plt.semilogy(z,tot_icm , c='navy', label='ICM total pressure') # darkblue
            plt.semilogy(z,pre2,c='blue',label='ICM thermal pressure',alpha=0.8)
            plt.semilogy(z,ICM,c='steelblue',label='ICM ram pressure',alpha=0.8)
        #plt.semilogy(xnew, ftot(xnew), c='darkblue', label='Total Pressure')
        #plt.semilogy(xnew,fpre(xnew),c='royalblue',label='Thermal Pressure')
        #plt.semilogy(xnew,fICM(xnew),c='deepskyblue',label='Ram Pressure')
        #plt.semilogy(xx, rp, c='plum', label='Restoring Pressure')
        ml = MultipleLocator(5)
        plt.minorticks_on()
        plt.axes().tick_params(which='minor', direction='in')
        plt.axes().tick_params(which='major', direction='in')
        if j != 0:
            plt.axhline(Model[j-1],linestyle='--',label='ICM Input',linewidth=1,c='k')
            plt.text(775 * 8, Model[j - 1] * 0.68, '%s' % Modell[j - 1])
            plt.text(150, 0.43 * 1e6, r'Scalar cut = %s' % s_cut)
        plt.axvline(448*8,linestyle='--',linewidth=0.5,c='k',alpha=0.7)
        plt.xticks([73*8, 198*8, 323*8, 448*8, 573*8, 698*8, 823*8],['-3', '-2', '-1', '0', '1', '2', '3'])
        plt.text(150,0.65*1e6,r'T= %s Myr' % round(tidx*unit['time'].value,1))
        plt.xlabel(r'z [kpc]')
        plt.ylabel(r'Pressure $[K cm^{-3}]$')
        plt.ylim(1e2,1e6)
        plt.xlim(z[0],z[-1])
        plt.legend(loc=0,labelspacing=0.1)
        plt.tight_layout()
        #plt.title('Pressure Variance_%s_%s' % (labell[j],tidx))
        if j!=0:
            plt.savefig('D:/yeongu/plots/new_data_pre_modi/new/Press_Vary_%s_%s_%s_w_newmag.png' % (s_cut,labell[j], tidx),dpi=200)
        else:
            plt.savefig('D:/yeongu/plots/new_data_pre_modi/new/Press_Vary_%s_%s_w_newmag.png' % (labell[j], tidx), dpi=200)

        #plt.show()
        plt.close()
'''
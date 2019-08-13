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
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
#print unit['density']
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K / erg = g cm2/s2K
vpc = 7168.*1024*1024/(128*128*896) # volume per cell

z = np.arange(-3588, 3588, 8)
g = 4.5181 * 1e-30 # gravitational constant : pc3/solmass*s2
sig_star = 42 # solmass/pc2
z_s = 245 # pc
r_0 = 8000 # pc
rho_dm = 0.0064 # solmass/pc3
g_z = (2. * np.pi * g * sig_star * (z ) / ((z ) ** 2 + z_s ** 2) ** (0.5) + 4. * np.pi * g * rho_dm * ((z ) / (1 + (z ) ** 2 / r_0 ** 2)))  # pc/s2

#plt.plot((z-3584.)/1000,g_z)
#plt.show()
meter = 3.24078*1e-17 # pc
kg = 5.02785*1e-31 # solar mass

stop = 500

simid_t = ('RPS_8pc_noICM_newacc','RPS_8pc_ICM0_newacc','RPS_8pc_ICM1_newacc','RPS_4pc_ICM1_newacc','RPS_8pc_ICM2_newacc','RPS_4pc_ICM2_newacc','RPS_8pc_ICM3_newacc')
labell = ('No ICM','P1', 'P3','P3h', 'P7','P7h','P14' ,'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
Model = [8.63*1e3,3.46*1e4,6.92*1e4,1.38*1e5]
Modell = [r'8.63*1e3',r'$3.46x10^4$',r'$6.92x10^4$',r'1.38*1e5']
C = ('k', 'salmon', 'mediumblue','deepskyblue' ,'darkgreen','lime', 'magenta','darkmagenta','goldenrod','royalblue','crimson') # 'plum','orchid','purple'
S = ('-.','--','-')

# overplot Starformation rate of three different simulations
k=1
plt.figure(figsize=(6, 11.5))
for j in (0,1,2,4,6) :
    basedir = 'G:/yeongu/'
    if j==5 or j==6:
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
    Ther_l= []
    Tot_e = []
    Tot_l = []
    #plt.figure(figsize=(6, 5))
    #stop = 256
    crit = 94
    for tidx in range(250, stop):  # time step 251, 331, 411, 501
        #plt.figure(figsize=(6, 5))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        #print(ds.field_list)

        #print(ds.derived_field_list)

        # full domain information
        #print ds.domain
        #rs = ds.read_all_data('reynold_stress')
        #print rs
        # information of grid #0
        #print ds.grids[0]
        T1 = ds.read_all_data('T1'); coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)
        #print temp
        if j != 0:
            scalar = ds.read_all_data('specific_scalar3') # ism = 0 / icm = 1
            s_cut = 0.5
        B = ds.read_all_data('magnetic_field')
        Pmag1 = B[:,:,:,0]**2/2
        Pmag2 = B[:,:,:,1]**2/2
        Pmag3 = B[:,:,:,2]**2/2
        mag_p = Pmag1+Pmag2-Pmag3
        mag_p *= (unit['pressure']/kb).cgs.value


        mag_p2 = copy.copy(mag_p)
        #print mag_p2
        if j != 0:
            mag_p[scalar > s_cut] = 0 # ISM
        mag_p_c = copy.copy(mag_p)
        #print mag_p[temp < 20000]
        mag_p[temp < 20000] = 0 ; mag_p_h = mag_p
        mag_p_c[temp > 20000] = 0
        if j != 0:
            mag_p2[scalar < s_cut] = 0 # ICM


        d1 = ds.read_all_data('density')*unit['density'].value # ISM density
        d2 = copy.copy(d1)#ds.read_all_data('density')*unit['density'].value  # ICM density
        d_mks_ism = d1*6.76745*1e-11
        d_mks_icm = d2*6.76745*1e-11 # solarmass/pc3 to kg/km3
        if j != 0:
            d_mks_ism[scalar > s_cut] = 0 # select only ISM
            d_mks_icm[scalar < s_cut] = 0 # select only ICM
        if j != 0:
            Area_icm = 1024*1024/8/8.
            #for i in range(896):
            #    Area_icm.append(np.count_nonzero(d_mks_icm[i, :, :]))
            #print Area_icm

        pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # ISM thermal pressure
        pre2 = copy.copy(pre)

        #pre2 = ds.read_all_data('pressure')*unit['pressure'].value/kb # ICM thermal pressure
        if j != 0:
            pre[scalar > s_cut] = 0 #select only ISM
        pre_c = copy.copy(pre) # ISM for cold

        pre[temp < 20000] = 0 ; pre_h = pre # hot ISM
        pre_c[temp > 20000] = 0 # cold ISM
        if j != 0:
            pre2[scalar < s_cut] = 0 # select only ICM

        Area=1024*1024/8/8.
        #for i in range(896):
        #    Area.append(np.count_nonzero(pre_c[i,:,:]))
        #print Area
        #print len(Area)

        vel = ds.read_all_data('velocity') ; vel2 = copy.copy(vel) #ds.read_all_data('velocity') # km/s
        vel_z = vel[:,:,:,2];  vel_z2 = vel2[:,:,:,2]
        if j != 0:
            vel_z[scalar > s_cut] = 0 # select only ISM
            vel_z2[scalar < s_cut] = 0 # select only ICM

        ############### ISM pressure ####################

        turb = d_mks_ism*(vel_z**2)*1e-2/kb # ISM turbulence pressure 1e-2 is conversion factor from mks to g cm
        turb_c = copy.copy(turb)

        turb_c[temp > 20000] = 0 # cold ISM
        turb[temp < 20000] = 0 ; turb_h=turb # hot ISM

        magp_c = np.sum(np.sum(mag_p_c, axis=1), axis=1)/Area
        turb_c = np.sum(np.sum(turb_c,axis=1),axis=1)/Area # ISM turbulence pressure / cold
        turb_h = np.sum(np.sum(turb_h,axis=1),axis=1)/Area # ISM hot
        pre_c = np.sum(np.sum(pre_c,axis=1),axis=1)/Area # ISM thermal pressure / cold
        pre_h = np.sum(np.sum(pre_h,axis=1),axis=1)/Area # ISM hot
        #m_pre_c = np.repeat(np.mean(np.mean(m_pre_c, axis=1), axis=1), 8) # magnetic pressure
        #m_pre_h = np.repeat(np.mean(np.mean(m_pre_h, axis=1), axis=1), 8) # magnetic pressure
        tot_c = turb_c + pre_c + magp_c
        ############### ICM pressure #####################
        if j!=0 :
            ICM = d_mks_icm*(vel_z2**2)*1e-2/kb # Ram Pressure : density x velocity^2, 1e-2 : convert to cgs
            ICM = np.sum(np.sum(ICM,axis=1),axis=1)/Area_icm # y projection
            pre2 = np.sum(np.sum(pre2,axis=1),axis=1)/Area_icm # ICM thermal pressure
            tot_icm = pre2 + ICM

        if tidx < crit+250 :
            Tot_e.append(tot_c)
            Turb_e.append(turb_c)
            Ther_e.append(pre_c)
            Mag_e.append(magp_c)
        else:
            Tot_l.append(tot_c)
            Turb_l.append(turb_c)
            Ther_l.append(pre_c)
            Mag_l.append(magp_c)
        print tidx
    Tot_e = np.array(Tot_e); Tot_l = np.array(Tot_l)
    Turb_e = np.array(Turb_e); Turb_l = np.array(Turb_l)
    Ther_e = np.array(Ther_e); Ther_l = np.array(Ther_l)
    Mag_e = np.array(Mag_e); Mag_l = np.array(Mag_l)

    Medi_e = np.nanmedian(Tot_e,axis=0); Medi_l = np.nanmedian(Tot_l,axis=0)
    Min_e = np.nanpercentile(Tot_e,25,axis=0); Min_l = np.nanpercentile(Tot_l,25,axis=0)
    Max_e = np.nanpercentile(Tot_e,75,axis=0); Max_l = np.nanpercentile(Tot_l,75,axis=0)
    #Min_e = np.min(Tot_e,axis=0); Min_l = np.min(Tot_l,axis=0)
    #Max_e = np.max(Tot_e,axis=0); Max_l = np.max(Tot_l,axis=0)

##################### Plot #############################
    zz = range(896)
    plt.subplot(5,1,k)

    plt.plot(zz,Medi_e,'c-',label='Early_%s' % labell[j])
    plt.fill_between(zz,Min_e,Max_e,facecolor=C[j],alpha=0.5)

    #plt.xticks([73 * 8, 198 * 8, 323 * 8, 448 * 8, 573 * 8, 698 * 8, 823 * 8], [])
    #plt.ylabel(r'Pressure $[K cm^{-3}]$')
    #plt.yscale('log', nonposy='clip')
    #plt.ylim(1e2, 1e5)
    #plt.xlim(0, 7168)
    #plt.tick_params(which='major', direction='in')
    #plt.tick_params(which='minor', direction='in')
    #plt.legend(loc='upper right')

    #plt.subplot(10,1,2*k)

    plt.plot(zz,Medi_l,'m--',label='Late_%s' % labell[j])
    plt.fill_between(zz,Min_l,Max_l,facecolor=C[j],alpha=0.5)

    plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    plt.ylabel(r'Pressure $[K cm^{-3}]$')
    if k==5:
        plt.xlabel(r'z [kpc]')
    plt.yscale('log', nonposy='clip')
    plt.ylim(1e2, 1e5)
    plt.xlim(0,896)
    plt.tick_params(which='major', direction='in')
    plt.tick_params(which='minor', direction='in')
    plt.legend(loc='upper right')

    k=k+1

plt.tight_layout()
plt.savefig('D:/yeongu/plots/paperplot/new/Press_new_w_mag_quartile_5.png')
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
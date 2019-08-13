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
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K / erg = g cm2/s2
vpc = 7168.*1024*1024/(128*128*896) # volume per cell

z = np.arange(0,7168,1)
g = 4.5181 * 1e-30 # gravitational constant : pc3/solmass*s2
sig_star = 42 # solmass/pc2
z_s = 245 # pc
r_0 = 8000 # pc
rho_dm = 0.0064 # solmass/pc3
g_z = 2.*np.pi*g*sig_star*(z-3584)/((z-3584)**2+z_s**2)**(0.5) + 4.*np.pi*g*rho_dm*((z-3584)/(1+(z-3584)**2/r_0**2)) #

meter = 3.24078*1e-17 # pc
kg = 5.02785*1e-31 # solar mass

stop = 500

simid_t = ('RPS_8pc_ICM1','RPS_8pc_ICM2') #'MHD_8pc_new',
labell = ('ICM1','ICM2') #'NonICM',
Model = [3.46*1e4,6.92*1e4]
Modell = [r'$3.46x10^4$',r'$6.92x10^4$']
S = ('-.','--','-')

# overplot Starformation rate of three different simulations

for j in range(len(simid_t)) :
    basedir = 'D:/yeongu/'
    simid = simid_t[j]
    Mom_up = []
    plt.figure(figsize=(5, 5))
    for tidx in range(250, 500):  # time step 251, 331, 411, 501

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # name of original data fields we stored from the simulation
        #print(ds.field_list)

        # It also has predefined data fields can be calculated from the original data.

        #print(ds.derived_field_list)

        # full domain information
        #print ds.domain

        # information of grid #0
        #print ds.grids[0]

        # yet, we didn't read data.
        # let's read each data field in a full domain

        scalar = ds.read_all_data('specific_scalar4') # ism = 0 / icm = 1
        d = ds.read_all_data('density')*unit['density'].value # density
        #d2 = ds.read_all_data('density') * unit['density'].value  # density
        d_ism = d#*6.76745*1e-11
        d_mks_icm = d*6.76745*1e-11 # solarmass/pc3 to kg/km3

        s_cut = 0.5

        d_ism = d_ism[scalar < s_cut]
        d_mks_icm = d_mks_icm[scalar >= s_cut]
        #mass = d # density times volume per cell = mass per cell
        pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
        pre = pre[scalar > s_cut]

        vel = ds.read_all_data('velocity') # km/s
        vel_z = vel[:,:,:,2]
        vel_z = vel_z[scalar > s_cut]
        #v_z_p = vel_z[vel_z > 0]
        #m_p = mass[vel_z> 0]
        #v_z_n = vel_z[vel_z < 0]
        #m_n = mass[vel_z < 0]
        ICM = d_mks_icm*(vel_z**2)*1e-2/kb # Ram Pressure : density x velocity^2, 1e-2 : convert to cgs
        ICM = np.repeat(np.mean(np.mean(ICM,axis=1),axis=1),8) # y projection
        pre = np.repeat(np.mean(np.mean(pre,axis=1),axis=1),8) # thermal pressure
        tot = pre + ICM

        ######### restoring pressure #############
        d = np.repeat(d_ism,8,axis=0)
        restP = d*g_z[:,None,None]
        restP = np.mean(np.mean(restP,axis=1),axis=1)*meter*10/kg/kb
        P = []
        for zz in range(7168) :
            P.append(np.sum(restP[zz:7168]))
        ##########################################

        plt.semilogy(z, P, c='crimson' , label='Restoring P')
        plt.semilogy(z,tot , c='navy', label='Total Pressure') # darkblue
        plt.semilogy(z,pre,c='blue',label='Thermal Pressure',alpha=0.8)
        plt.semilogy(z,ICM,c='steelblue',label='Ram Pressure',alpha=0.8)
        #plt.semilogy(xnew, ftot(xnew), c='darkblue', label='Total Pressure')
        #plt.semilogy(xnew,fpre(xnew),c='royalblue',label='Thermal Pressure')
        #plt.semilogy(xnew,fICM(xnew),c='deepskyblue',label='Ram Pressure')
        #plt.semilogy(xx, rp, c='plum', label='Restoring Pressure')
        ml = MultipleLocator(5)
        plt.minorticks_on()
        plt.axes().tick_params(which='minor', direction='in')
        plt.axes().tick_params(which='major', direction='in')

        plt.axhline(Model[j],linestyle='--',label='Initial Input',linewidth=1,c='k')
        plt.axvline(448*8,linestyle='--',linewidth=0.5,c='k')
        plt.text(775*8,Model[j]*0.68,'%s' % Modell[j])
        plt.xticks([73*8, 198*8, 323*8, 448*8, 573*8, 698*8, 823*8],['-3', '-2', '-1', '0', '1', '2', '3'])
        plt.text(150,0.65*1e6,r'T= %s Myr' % round(tidx*unit['time'].value,1))
        plt.text(150,0.43*1e6,r'Scalar cut = %s' % s_cut)
        plt.xlabel(r'z [kpc]')
        plt.ylabel(r'Pressure $[K cm^{-3}]$')
        plt.ylim(1e2,1e6)
        plt.xlim(z[0],z[-1])
        plt.legend(loc=0,labelspacing=0.1)
        #plt.title('Pressure Variance_%s_%s' % (labell[j],tidx))
        plt.show()
        #plt.savefig('D:/yeongu/plots/new_data_pre_modi/Press_Vary_%s_%s_%s.png' % (s_cut,labell[j], tidx),dpi=400)

        plt.close()
        print tidx
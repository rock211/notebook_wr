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
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K / erg = g cm2/s2
vpc = 7168.*1024*1024/(128*128*896) # volume per cell

z = np.arange(0,7168,1)
g = 4.5181 * 1e-30 # gravitational constant : pc3/solmass*s2
sig_star = 42 # solmass/pc2
z_s = 245 # pc
r_0 = 8000 # pc
rho_dm = 0.0064 # solmass/pc3
g_z = 2.*np.pi*g*sig_star*(z-3584)/((z-3584)**2+z_s**2)**(0.5) + 4.*np.pi*g*rho_dm*((z-3584)/(1+(z-3584)**2/r_0**2)) #

#plt.plot((z-3584.)/1000,g_z)
#plt.show()
meter = 3.24078*1e-17 # pc
kg = 5.02785*1e-31 # solar mass

simid_t = ('RPS_8pc_noICM_newacc','RPS_8pc_ICM0_newacc','RPS_8pc_ICM1_newacc','RPS_4pc_ICM1_newacc','RPS_8pc_ICM2_newacc','RPS_4pc_ICM2_newacc','RPS_8pc_ICM3_newacc')
#labell = ('No ICM','Very Weak' ,'Weak', 'Strong', 'Very Strong','ICM1', 'ICM2', 'ICM3', 'ICM4') #'NonICM',
labell = ('No ICM','P1', 'P3','P3h', 'P7','P7h','P14' ,'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
Model = [8.63*1e3,3.46*1e4,6.92*1e4,1.38*1e5]
Modell = [r'8.63*1e3',r'$3.46x10^4$',r'$6.92x10^4$',r'1.38*1e5']
S = ('-.','--','-')
C = ('k', 'salmon', 'mediumblue','deepskyblue' ,'darkgreen','lime', 'magenta','darkmagenta','goldenrod','royalblue','crimson') # 'plum','orchid','purple'
# overplot Starformation rate of three different simulations

Myr = unit['time'].to('Myr').value
Msun = unit['mass'].to('Msun').value
agebin = 10
crit = 94
for j in (0,1,2,4,6) :
    basedir = 'G:/yeongu/'
    if j==4 or j==6 :
        stop = 473
    elif j==5:
        stop = 389
    else:
        stop = 499

    simid = simid_t[j]
    Mom_up = []
    P=[]
    SFR=[]
    #plt.figure(figsize=(6, 5))
    for tidx in range(250, stop):  # time step 251, 331, 411, 501
        #plt.figure(figsize=(6, 5))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        #print(ds.field_list)

        #print(ds.derived_field_list)

        starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)

        star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
        star_clu2 = star_clu[star_clu['age'] * Myr < agebin]

        M_star = sum(star_clu2['mass']) * Msun
        sfr = M_star * unit['time'].value / (1e+6 * agebin * 1.024 * 1.024)

        # full domain information
        #print ds.domain
        #rs = ds.read_all_data('reynold_stress')
        #print rs
        # information of grid #0
        #print ds.grids[0]
        T1 = ds.read_all_data('T1'); coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)

        if j != 0:
            scalar = ds.read_all_data('specific_scalar3') # ism = 0 / icm = 1
            s_cut = 0.5

        d1 = ds.read_all_data('density')*unit['density'].value # ISM density
        d_mks_ism = d1*6.76745*1e-11
        if j != 0:
            d_mks_ism[scalar > s_cut] = 0 # select only ISM

        vel = ds.read_all_data('velocity')
        vel_z = vel[:,:,:,2];
        if j != 0:
            vel_z[scalar > s_cut] = 0 # select only ISM

        pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # ISM thermal pressure
        if j != 0:
            pre[scalar > s_cut] = 0 #select only ISM

        B = ds.read_all_data('magnetic_field')
        Pmag1 = B[:,:,:,0]**2/2
        Pmag2 = B[:,:,:,1]**2/2
        Pmag3 = B[:,:,:,2]**2/2
        mag_p = Pmag1+Pmag2-Pmag3
        mag_p *= (unit['pressure']/kb).cgs.value
        if j != 0:
            mag_p[scalar > s_cut] = 0 #select only ISM
        ############### ISM pressure ####################

        turb = d_mks_ism*(vel_z**2)*1e-2/kb # ISM turbulence pressure

        P_ism = np.mean(np.mean(turb+pre+mag_p,axis=1),axis=1)

        P.append(P_ism.max())
        SFR.append(sfr)
        print tidx


    SFR_e=np.array(SFR[:crit])
    SFR_l=np.array(SFR[crit:])
    P_e=np.array(P[:crit])
    #print P_e
    P_l=np.array(P[crit:])
    #print P_l
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
    '''
    ############################### Upper and Lower #########################################
    SFR_e_errors = [np.mean(SFR_e)-np.percentile(SFR_e,25),np.percentile(SFR_e,75)-np.mean(SFR_e)] # [np.nonzero(SFR_e)] nonzero selection
    print SFR_e_errors
    SFR_l_errors = [np.mean(SFR_l)-np.percentile(SFR_l,25),np.percentile(SFR_l,75)-np.mean(SFR_l)]
    print SFR_l_errors
    P_e_errors = [np.mean(P_e)-np.percentile(P_e,25),np.percentile(P_e,75)-np.mean(P_e)]
    print np.array([P_e_errors]).T
    P_l_errors = [np.mean(P_l)-np.percentile(P_l,25),np.percentile(P_l,75)-np.mean(P_l)]
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
plt.xlim(3e2,3e5)
plt.ylim(10**(-4.5),10**(-1.5))
plt.ylabel('log(SFR)')
plt.xlabel('log(P/kb)')
#plt.legend(loc=0)
plt.tight_layout()
plt.savefig('D:/yeongu/plots/paperplot/new/SFR_P_wz_quartile.png',dpi=300)
plt.savefig('D:/yeongu/plots/paperplot/new/SFR_P_wz_quartile.eps',format='eps',dpi=300)
plt.show()
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
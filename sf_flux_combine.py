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
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
#print unit['density']
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K
volpercell = 7168.*1024*1024/(128*128*896)
vpc = volpercell

# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

simid_t = ('R8_8pc_metal','RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1','RPS_8pc_ICM2','RPS_8pc_ICM3','RPS_8pc_ICM4') # 'MHD_8pc_new' ,
labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # r'No ICM',
labelll = ('Cold','Unstable','Warm','Ionized','Hot')
#C = ('goldenrod','royalblue','firebrick')
C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')
C3 = ('darkred','red','salmon')
C = ('k','r','g','b','cyan','magenta','y')
S = (':','--','-')

Msun = unit['mass'].to('Msun').value
print Msun
Myr=unit['time'].to('Myr').value
# overplot Starformation rate of three different simulations

for j in range(4,6) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j]
    Mom_up = []

    if j == 5:
        stop = 472
    elif j == 6:
        stop = 401
    else:
        stop = 500

    Mom_1 = []
    Mom_2 = []
    Mom_3 = []
    Mom_c = []
    Mom_u = []
    Mom_w = []
    Mom_h = []
    Mom_i = []
    Mom_cuw = []
    Mom_hi = []
    Mom_t = []
    Mom_icm = []
    Mj = []
    plt.figure(figsize=(8, 5)) # shrink 12,5
    for tidx in range(250, stop):  # time step 251, 331, 411, 501

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

        # this can be original data fields

        #gp = ds.read_all_data('gravitational_potential')
        #print unit['gravitational_potential']
        #print gp
        d = ds.read_all_data('density')*unit['density'].value # density
        #den = ds.read_all_data('density')
        #P = ds.read_all_data('pressure')
        #mass = d # density times volume per cell = mass per cell
        #pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
        vel = ds.read_all_data('velocity')
        T1 = ds.read_all_data('T1')
        vel_z = vel[:,:,:,2]
        coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func
        #temp = ds.read_all_data('temperature')
        #v_z_p = vel_z[vel_z > 0]
        #m_p = mass[vel_z> 0]
        #v_z_n = vel_z[vel_z < 0]
        #m_n = mass[vel_z < 0]

        ###### Mass #######
        height = 895 #698 823 #895 #data set

        d_up = d[height]; d_down = d[0]
        vel_z_up = vel_z[height]; vel_z_down = vel_z[0]
        temp_up = temp[height]; temp_down = temp[0]

        #print np.median(temp_down), np.mean(temp[50])
        #d_c = d_up[temp_up < 184] ; v_c = vel_z_up[temp_up < 184]
        #d_u = d_up[(temp_up > 184) & (temp_up < 5050)] ; v_u = vel_z_up[(temp_up > 184) & (temp_up < 5050)]
        #d_w = d_up[(temp_up > 5050) & (temp_up < 20000)] ; v_w = vel_z_up[(temp_up > 5050) & (temp_up < 20000)]
        #d_h = d_up[(temp_up > 20000) & (temp_up < 500000)] ; v_h = vel_z_up[(temp_up > 20000) & (temp_up < 500000)]
        #d_i = d_up[temp_up > 500000] ; v_i = vel_z_up[temp_up > 500000]
        d_cuw = d_up[temp_up < 20000] ; v_cuw = vel_z_up[temp_up < 20000]
        d_hi = d_up[temp_up > 20000];   v_hi = vel_z_up[temp_up > 20000]

        d_d_cuw = d_down[temp_up < 20000] ; v_d_cuw = vel_z_down[temp_up < 20000]
        d_d_hi = d_down[temp_up > 20000]; v_d_hi = vel_z_down[temp_up > 20000]

        m_cuw = np.sum(np.sum(d_cuw * v_cuw / unit['time'].value)) * 512 / (8 * 1024 * 1024)
        m_hi = np.sum(np.sum(d_hi * v_hi / unit['time'].value)) * 512 / (8 * 1024 * 1024)


        #tot = m_cuw+m_hi
        #Mom_cuw.append(m_cuw); Mom_hi.append(m_hi)
        if j==0:
            m_d_cuw = np.sum(np.sum(d_d_cuw * v_d_cuw / unit['time'].value)) * 512 / (8 * 1024 * 1024)
            m_d_hi = np.sum(np.sum(d_d_hi * v_d_hi / unit['time'].value)) * 512 / (8 * 1024 * 1024)
            #print m_d_cuw+m_d_hi
            Mom_cuw.append(m_cuw - m_d_cuw); Mom_hi.append(m_hi - m_d_hi)
        else:
            scalar = ds.read_all_data('specific_scalar4')  # ism = 0 / icm = 1
            scalar_up = scalar[height]
            m_icm = np.sum(np.sum(scalar_up * d_up * vel_z_up / unit['time'].value)) * 512 / (8 * 1024 * 1024)
            Mom_cuw.append(m_cuw); Mom_hi.append(m_hi); Mom_icm.append(m_icm)

        agebin = 10
        starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)
        star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
        star_clu2 = star_clu[star_clu['age'] * Myr < agebin]
        #M_star += sum(star_clu2['mass'])*Msun # mass sum in time scale / cumulative
        M_star = sum(star_clu2['mass']) * Msun # SFR time variation
        Mj.append(M_star* unit['time'].value / (1e+6*agebin*1.024*1.024)) #* unit['time'].value/ mass divided by time scale & area
        #sfe = M_star* unit['time'].value / (1e+6*agebin*1.024*1.024)/surf_mass

        print tidx
        #print len(vel_z)

    time = np.arange(250*unit['time'].value,stop*unit['time'].value,unit['time'].value)
    #xnew = np.linspace(time.min(),time.max(),5000)
    #smooth = spline(time,Mom_up,xnew)
    plt.subplot(211)
    plt.plot(time, Mom_cuw, label='Cold/Warm', color=C2[0])
    plt.plot(time, Mom_hi, label='Hot', color=C2[3])
    icm_c = 'salmon'
    if j == 1 : # ICM 00
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c)
        icm_inflow = 0.0051002 / 16 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls='--',color='k')
    elif j== 2: # ICM 0
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c)
        icm_inflow = 0.0051002 / 4 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls='--',color='k')
    elif j== 3: # ICM 1
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c)
        icm_inflow = 0.0051002  # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls='--',color='k')
    elif j== 4: # ICM 2
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c)
        icm_inflow = 0.0051002 * 2 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls='--',color='k')
    elif j== 5: # ICM 3
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c)
        icm_inflow = 0.0051002 * 4 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls='--',color='k')
    elif j== 6: # ICM 4
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c)
        icm_inflow = 0.0051002 * 10 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls='--',color='k')

    #plt.plot(time, Mom_t,label='Total', color='dimgrey',ls=':')
    #plt.plot(time, Mom_hi, label='Hot_%s' % labell[j], color=C3[j],ls=S[j]) # Hot only
    #plt.plot(time, Mom_c, label='%s' % labelll[0], color=C2[0])
    #plt.plot(time, Mom_u, label='%s' % labelll[1], color=C2[1])
    #plt.plot(time, Mom_w, label='%s' % labelll[2], color=C2[2])
    #plt.plot(time, Mom_h, label='%s' % labelll[3], color=C2[3])
    #plt.plot(time, Mom_i, label='%s' % labelll[4], color=C2[4])
    #plt.plot(xnew,smooth,label='%s' % labell[j],color=C[j]) ; #plt.plot(time,Mom_3,label='3kpc') ;  plt.plot(time,Mom_2,label='2kpc'), plt.plot(time,Mom_1,label='1kpc')

    #ml = MultipleLocator(5)
    #plt.minorticks_on()
    #plt.axes().tick_params(which='minor', direction='in')
    #plt.axes().tick_params(which='major', direction='in',labelsize=22)

    #plt.xlabel(r'Time [Myr]',fontsize=18)
    plt.xticks([])
    plt.ylabel(r'Mass Flux [$M_{\odot}$ pc$^{-2}$ Myr$^{-1}$]',fontsize=10)
    plt.legend(loc=0)
    plt.xlim(time.min(),time.max())
    #plt.title('%s' % labell[j])

    plt.subplot(212)
    #plt.semilogy(xx,Mj,label='%s' % labell[j],color=C[j],linestyle=S[j],linewidth=L[j])
    plt.plot(time, Mj, label='%s' % labell[j], color=C[j])
    #plt.plot(xx, SFE, label='%s' % labell[j], color=C[j], linestyle=S[j], linewidth=L[j], alpha=alpha[j])
    #ml = MultipleLocator(5)
    #plt.minorticks_on()
    #plt.axes().tick_params(which='minor', direction='in')
    #plt.axes().tick_params(which='major', direction='in', labelsize=16)
    plt.xlabel(r'Time [Myr]', fontsize=14)
    # plt.ylabel(r'SFE',fontsize=17) # For SFE
    plt.ylabel(r'$\Sigma_{SFR}$ [M$_{\odot}$ kpc$^{-2}$ yr$^{-1}$]',fontsize=10) # For SFR
    #plt.ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ pc$^{-2}$]', fontsize=19)  # For cumulative
    # plt.title(r'SFR variation $(Agebin=%s)$' % agebin)
    # plt.yticks([0,0.005,0.01,0.015,0.02]) # For cumulative
    # plt.yticks([0,0.005,0.01,0.015,0.02]) # For SFR
    plt.xlim(time.min(), time.max())
    plt.ylim(0,0.019)
    #plt.ylim(0, 1.1)  # For cumulative
    plt.legend(loc=0, fontsize=10)
    #if j == 1 :
    #    plt.legend(loc=0,fontsize=18)
    #if j != 2 :
    #    plt.ylim(-0.0002,0.0065)

    plt.tight_layout(pad=0.3)
    plt.savefig('D:/yeongu/plots/sf_massflux_%s.png' % labell[j],dpi=500) # ,transparent=True
    #plt.show()
    plt.close()



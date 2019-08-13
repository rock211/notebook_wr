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
Myr=unit['time'].to('Myr').value
# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.





simid_t = ('RPS_8pc_noICM_newacc','RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1_newacc','RPS_8pc_ICM2_newacc','RPS_8pc_ICM3','RPS_8pc_ICM4') # 'MHD_8pc_new' ,
labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # r'No ICM',
labelll = ('Cold','Unstable','Warm','Ionized','Hot')
C = ('goldenrod','royalblue','firebrick')
C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')
C3 = ('darkred','red','salmon')
S = (':','--','-')

k = 1
# overplot Starformation rate of three different simulations
#plt.figure(figsize=(14,8))
#fig =plt.figure(figsize=(8.5,12))
for j in (0,3,4): #range(1,7) :
    basedir = 'G:/yeongu/'
    simid = simid_t[j]
    Mom_up = []

    if j == 5:
        stop = 472
    elif j == 6:
        stop = 401
    else:
        stop = 499
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
    #plt.figure(figsize=(8, 5)) # shrink 12,5

    for tidx in range(250, stop):  # time step 251, 331, 411, 501
        plt.figure(figsize=(8,8))
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
        ###### momentum #######
        #moment_z = d*vel_z/unit['time'].value #SolMass/pc2/Myr
        #m_c = np.sum(np.sum(d_c * v_c / unit['time'].value))*512/(8*1024*1024)
        #m_u = np.sum(np.sum(d_u * v_u / unit['time'].value))*512/(8*1024*1024)
        #m_w = np.sum(np.sum(d_w * v_w / unit['time'].value))*512/(8*1024*1024)
        #m_h = np.sum(np.sum(d_h * v_h / unit['time'].value))*512/(8*1024*1024)
        #m_i = np.sum(np.sum(d_i * v_i / unit['time'].value))*512/(8*1024*1024)
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
            scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1
            scalar_up = scalar[height]
            m_icm = np.sum(np.sum(scalar_up * d_up * vel_z_up / unit['time'].value)) * 512 / (8 * 1024 * 1024)
            Mom_cuw.append(m_cuw); Mom_hi.append(m_hi); Mom_icm.append(m_icm)
        print tidx
        #print len(vel_z)

        plt.plot(Mom_cuw, label='Cold/Warm', color=C2[0],lw=3)
        plt.plot(Mom_hi, label='Hot', color=C2[3],lw=3)
        icm_c = 'green'
        if j == 0:
            #plt.ylim(-0.00031,0.0062)
            plt.ylim(-0.00445, 0.0062)
        if j == 1 : # ICM 00
            plt.plot(Mom_icm, label='ICM outflow', color=icm_c,ls='--')
            icm_inflow = 0.0051002 / 16 # m_d_cuw + m_d_hi
            plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
            plt.ylim(-0.000275, 0.0019)
        elif j== 2: # ICM 0
            plt.plot(Mom_icm, label='ICM outflow', color=icm_c,ls='--')
            icm_inflow = 0.0051002 / 4 # m_d_cuw + m_d_hi
            plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
            plt.ylim(-0.000275, 0.0019)
        elif j== 3: # ICM 1
            plt.plot(Mom_icm, label='ICM outflow', color=icm_c,ls='--')
            icm_inflow = 0.0051002  # m_d_cuw + m_d_hi
            plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
            #plt.ylim(-0.000275, 0.0055)
            #plt.ylim(-0.00031, 0.0062) # for 3 kpc
            plt.ylim(-0.0325, 0.027) # for 2 kpc
        elif j== 4: # ICM 2
            plt.plot(Mom_icm, label='ICM outflow', color=icm_c,ls='--')
            icm_inflow = 0.0051002 * 2 # m_d_cuw + m_d_hi
            plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
            #plt.ylim(-0.0175, 0.39)
            #plt.ylim(-0.00655, 0.131) # for 3 kpc
            plt.ylim(-0.033, 0.131) # for 2 kpc
        elif j== 5: # ICM 3
            plt.plot(Mom_icm, label='ICM outflow', color=icm_c,ls='--')
            icm_inflow = 0.0051002 * 4 # m_d_cuw + m_d_hi
            plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
            plt.ylim(-0.0175, 0.39)
        elif j== 6: # ICM 4
            plt.plot(Mom_icm, label='ICM outflow', color=icm_c,ls='--')
            icm_inflow = 0.0051002 * 10 # m_d_cuw + m_d_hi
            plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
            plt.ylim(-0.0175, 0.39)



        plt.xlabel(r'time [Myr]',fontsize=24)
        plt.ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]', fontsize=28)


        #if j==3:
        #    plt.ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]', fontsize=19)
        plt.yticks(fontsize=19)
        if j!=4:
            plt.ylim(-0.0002,0.0062)
        else:
            plt.ylim(-0.004,0.13)
        plt.xticks([(250 / Myr - 250), (300 / Myr - 250), (350 / Myr - 250), (400 / Myr - 250), (450 / Myr - 250)],
                   [250, 300, 350, 400, 450],fontsize=20)
        plt.xlim(0,250)
        plt.tick_params(direction='in')
        plt.legend(loc=2, fontsize=22)
        plt.tight_layout()
        #plt.show()
        plt.savefig('D:/yeongu/plots/flux_time/massflux_%s_%s.png' % (labell[j],tidx), dpi=200)
        plt.close()

#plt.close()



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
from scipy.interpolate import spline
from matplotlib.ticker import MultipleLocator

sys.path.insert(0, '../')
from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:

import pyathena as pa


unit = pa.set_units(muH=1.4271)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
#print (unit['mass'], unit['time'], unit['magnetic_field'])

Msun = unit['mass'].to('Msun').value
Myr=unit['time'].to('Myr').value

agebin = 10 # unit : Myr, 10 : H-alpha like, 40 : cluster lifetime

M_T_0 = []
M_T_1 = []
M_T_2 = []

plt.subplot()
simid_t = ('R8_8pc_metal','RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1','RPS_8pc_ICM2','RPS_8pc_ICM3','RPS_8pc_ICM4')
labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4')
C = ('k','r','g','b','cyan','magenta','y')
S = ('-','-','-')
L = (1.5,1.5,1.5)

# overplot Starformation rate of three different simulations

plt.figure(figsize=(7, 5))
for j in range(len(simid_t)) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j]

    if j == 5:
        stop = 472
    elif j == 6:
        stop = 401
    else:
        stop = 500
    print simid

    Mj = []
    M_star = 0.0
    o_z = []
    y_z = []
    y_z_t = []
    y_z_b = []

    for tidx in range(250, stop):
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir,simid, simid, tidx)

        # read in domain information
        #ds = pa.AthenaDataSet(vtkfname)

        # full domain information
        #print ds.domain

        # information of grid #0
        #ds.grids[0]

        # yet, we didn't read data.
        # let's read each data field in a full domain

        # this can be original data fields
        #d = ds.read_all_data('density')*unit['density'].value
        #nd = ds.read_all_data('number_density')
        #tem = ds.read_all_data('temperature')
        #coolftn = pa.coolftn()
        #temp = coolftn.get_temp(P / d)  # Temperature derived from P/d & cooling func

        starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)

        #print(sp.shape)
        #print (sp)


        print(tidx)
        star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
        star_clu_y = star_clu[star_clu['age'] * Myr < agebin]
        star_clu_o = star_clu[(star_clu['age'] * Myr > agebin) & (star_clu['age'] * Myr <= 40)]
        #print 'young' , star_clu_y, 'old' , star_clu_o

        young_z = star_clu_y['x3']
        #print young_z.min(), young_z.max()
        young_m = star_clu_y['mass']
        old_z = star_clu_o['x3']
        old_m = star_clu_o['mass']
        if np.sum(young_m)==0 :
            massw_z_y = np.nan
        else :
            massw_z_y = np.sum(young_z*young_m)/np.sum(young_m)

        if np.sum(old_m)==0 :
            massw_z_o = np.nan
        else :
            massw_z_o = np.sum(old_z * old_m) / np.sum(old_m)
        #print young_z
        #print M_star[0]
        #M_star = sum(star_clu2['mass'])*Msun # mass sum in time scale / cumulative
        #M_star = sum(star_clu2['mass']) * Msun # SFR time variation
        #M_t = sum(star_clu['mass']) * Msun

        #Mj.append(M_star / (1e+6*agebin*1.024*1.024))
        #o_z.append(np.mean(old_z))
        o_z.append(massw_z_o)
        #y_z.append(np.mean(young_z))
        #y_z_t.append(massw_z_y/1000.+young_z.max()/1000.)
        #y_z_b.append(massw_z_y/1000.-young_z.min()/1000.)
        y_z.append(massw_z_y/1000.)

    time = np.arange(250 * unit['time'].value, stop * unit['time'].value, unit['time'].value)
    #plt.plot(time,o_z,c=C[j],ls='--',label='Old')
    plt.plot(time,y_z,c=C[j],ls='-',label=labell[j])
    #plt.fill_between(time,y_z_b,y_z_t,color=C[j],alpha=0.5,edgecolor=None)
    plt.title('Agebin = %s Myr' % agebin)
    ml = MultipleLocator(5)
    plt.minorticks_on()
    plt.axes().tick_params(which='minor', direction='in')
    plt.axes().tick_params(which='major', direction='in')
    #plt.yticks([-600,-400,-200,0,200,400,600],[0.6,-0.4,-0.2,0,0.2,0.4,0.6])
    plt.xlabel(r'Time [Myr]')
    plt.ylabel(r'z [kpc]')
    plt.xlim(time.min(), 500)
    plt.ylim(-1,3)
    plt.legend(loc=1)
#plt.savefig('D:/yeongu/plots/new_data/mean_height_weighted_agemodi.png',dpi=500)
    #plt.close()
plt.show()
plt.close()

    #plt.show()
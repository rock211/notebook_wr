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


simid_t = ('RPS_8pc_noICM_newacc','RPS_8pc_ICM0_newacc','RPS_8pc_ICM1_newacc','RPS_8pc_ICM2_newacc','RPS_8pc_ICM3_newacc')
#labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4')
labell = ('No ICM',r'Very Weak','Weak','Strong',r'Very Strong')
C = ('k','salmon','deepskyblue','crimson','royalblue','darkmagenta','goldenrod','royalblue','crimson')
C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')# 'plum','orchid','purple'
M = ((5,1),'o','o','o','o','o','o')
S = ('-','-','-')
L = (1.5,1.5,1.5)
a = (0.8,0.8,0.8,0.8,0.8,0.6,1)

# overplot Starformation rate of three different simulations


for j in (3,4) :
    basedir = 'G:/yeongu/'
    simid = simid_t[j]

    if j == 4:
        stop = 474
    else:
        stop = 499
    print simid

    M_c = []
    M_i = []
    M_h = []
    plt.figure(figsize=(5.25, 5))
    for tidx in range(250, stop):
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir,simid, simid, tidx)

        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # full domain information
        #print ds.domain

        # information of grid #0
        #ds.grids[0]

        # yet, we didn't read data.
        # let's read each data field in a full domain

        # this can be original data fields
        d = ds.read_all_data('density')*unit['density'].value
        vel = ds.read_all_data('velocity')
        vel_z = vel[:,:,:,2]
        T1 = ds.read_all_data('T1'); coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)
        if j != 0:
            scalar = ds.read_all_data('specific_scalar3') # ism = 0 / icm = 1
            s_cut = 0.5
            d[scalar > s_cut] = 0  # select only ISM
            vel_z[scalar > s_cut] = 0
        d_cw = d[temp < 20000]; v_cw = vel_z[temp < 20000]
        d_i = d[(temp > 20000) & (temp < 500000)]; v_i = vel_z[(temp > 20000) & (temp < 500000)]
        d_h = d[temp > 500000]; v_h = vel_z[temp > 500000]
        #print d_cw.max(),d_ih.max()
        #print v_cw.max(),v_ih.max()
        mom_cw = np.mean(d_cw*v_cw)
        mom_i = np.mean(d_i*v_i)
        mom_h = np.mean(d_h*v_h)

        M_c.append(mom_cw); M_i.append(mom_i) ; M_h.append(mom_h)

        print tidx


    time = np.arange(250*unit['time'].value,stop*unit['time'].value,unit['time'].value)
    M_t = np.add(np.add(M_c,M_i),M_h)

    print M_t

    plt.subplot(4, 1, 1)
    plt.plot(time, M_t, label='Total', color=C[0])
    plt.xlim(time.min(), time.max())
    plt.xticks([250, 300, 350, 400, 450], [])
    plt.yticks(fontsize=19)
    plt.legend(loc=0)
    plt.tight_layout()

    plt.subplot(4, 1, 2)
    plt.plot(time, M_c, label='Cold/Warm', color=C2[0])
    #plt.ylabel(r'Mean Momentum', fontsize=24)
    plt.xlim(time.min(), time.max())
    plt.xticks([250, 300, 350, 400, 450], [])
    plt.yticks(fontsize=19)
    plt.legend(loc=0)
    plt.tight_layout()

    plt.subplot(4, 1, 3)
    plt.plot(time, M_i, label='Ion', color=C2[2])
    plt.xlim(time.min(), time.max())
    plt.xticks([250, 300, 350, 400, 450], [])
    plt.yticks(fontsize=19)
    plt.legend(loc=0)
    plt.tight_layout()

    plt.subplot(4, 1, 4)
    plt.plot(time, M_h, label='Hot', color=C2[3])
    plt.xlabel(r'time [Myr]', fontsize=20)
    plt.xticks([250, 300, 350, 400, 450], fontsize=20)
    plt.yticks(fontsize=19)
    plt.xlim(time.min(), time.max())
    plt.legend(loc=0)

    plt.tight_layout()
    plt.savefig('D:/yeongu/plots/%s_net_mom.png' % labell[j] , dpi=300)
    #plt.savefig('D:/yeongu/plots/%s_net_mom.eps' % labell[j],format='eps',dpi=300)
    #plt.show()
    plt.close()
    #


#plt.savefig('D:/yeongu/plots/paperplot/SFR_height_%s_modi_new.png' % j , dpi=300)
#plt.savefig('D:/yeongu/plots/paperplot/SFR_height_new.eps', format='eps',dpi=300)
#plt.show()


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys

sys.path.insert(0, '../')
from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:

import pyathena as pa

# In[8]:

basedir = 'D:/yeongu/'
simid = 'RPS_8pc_ICM1' #('RPS_8pc_ICM1','RPS_8pc_ICM2')

# ## Unit system
#
# The unit system we choose for this simulation is
# * [length] = pc
# * [velocity] = km/s
# * [density] = 1.4271*m_h/cm^3

# In[9]:

# You can retrive the unit system with pa.set_units function.
# To make unit conversion easier, I use astropy's unit and constant.
# You may need to install astropy to use it
# please visit http://www.astropy.org/

unit = pa.set_units(muH=1.4271)
print(unit)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
#print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

# In[10]:

# import pyathena as pa

C = []
U = []
W = []
H = []
I = []
CU = []
HI = []

for tidx in range(250, 500): #time step

    vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
    # read in domain information
    ds = pa.AthenaDataSet(vtkfname)

    # name of original data fields we stored from the simulation
    #print(ds.field_list)

    # It also has predefined data fields can be calculated from the original data.
    #print(ds.derived_field_list)

    # full domain information
    ds.domain

    # information of grid #0
    ds.grids[0]

    # yet, we didn't read data.
    # let's read each data field in a full domain

    # this can be original data fields
    d = ds.read_all_data('density')
    #print(d.shape)
    #print(d)
    # note that it adopts C-like indexing, k (z-index) comes first and i (x-index) comes last
    # vector field has 3 component

    #nd = ds.read_all_data('number_density')
    #print (nd.shape)
    #T = ds.read_all_data('T1')
    #print(T)
    #tem = ds.read_all_data('temperature') # Temperature from data
    #print (tem.shape)

    # calculate sound speed
    #P = ds.read_all_data('pressure')
    #cs = np.sqrt(P / d)
    T1 = ds.read_all_data('T1')
    # calculation of temperature needs additional information about mean molecular weight I adopt
    coolftn = pa.coolftn()
    tem = coolftn.get_temp(T1) # Temperature from P/d
    # print(temp)
    #print(temp[np.argmin(temp)])
    #print(temp[np.argmax(temp)])
    #print((max(temp)[0]),(max(temp)[1]))
    #hh = np.a
    # min(temp,axis=0)
    #
    #print(temp[15])
    #print(tem[15])
    #print(temp[15]/tem[15])

    ##### Phase dividing * Calculate mass ##### Need exact value corresponds to tem. bound


    #print(d_z.shape)
    #print(np.amax(temp[z]),np.amin(temp[z]))
    #print(temp[z])
    #temp[z] = temp[z] - 58000.
    #print(temp[z])
    mc = np.sum(d[tem < 184])
    #print mc
    mu = np.sum(d[(tem > 184) & (tem < 5050)])
    #mcu = np.sum(d_z[tem[z] < 5050])
    mw = np.sum(d[(tem > 5050) & (tem < 20000)])
    mh = np.sum(d[(tem > 20000) & (tem < 500000)])
    #mhi = np.sum(d_z[tem[z] > 20000])
    mi = np.sum(d[tem > 500000])

    Mtot = np.sum(mc + mu + mw + mh + mi) #Total Mass at each z
    #Mtot = np.sum(mcu + mw + mhi)
    C.append(mc/Mtot)
    U.append(mu/Mtot)
    #CU.append(mcu/Mtot)
    W.append(mw/Mtot)
    #HI.append(mhi/Mtot)
    H.append(mh/Mtot)
    I.append(mi/Mtot)
    print tidx
        #print(z)
    #print(C)
    #z = range(448)
stop = 500
xx = np.arange(250*unit['time'].value,stop*unit['time'].value,unit['time'].value)
plt.plot(xx, C,label='Cold')
plt.plot(xx, U,label='Unstable')
#plt.plot(z,CU,label='CU')
plt.plot(xx, W,label='Warm')
#plt.plot(z,HI,label='HI')
plt.plot(xx, H,label='Ionized')
plt.plot(xx, I,label='Hot')
plt.legend(loc=6)
print(tidx)
plt.title('Mass fraction_cool_%s' % simid)
plt.xlabel(r'Time [Myr]')
# plt.xlabel('z_direction')
plt.ylabel('Mass_fraction')
#plt.show()
plt.savefig('D:/yeongu/plots/Mass_frac_n1.png',dpi=500)
#plt.clf()
    #plt.show()
    #print(len(C),len(U),len(W),len(H),len(I))
    #print(Mtot_tidx)


    #print(temp_z.shape)
    #print(np.mean(temp[0]))
    #temp_z1 = temp_z[temp[0] > 58270]
    #temp_z2 = temp_z[temp[0] < np.mean(temp[0])]
    #print(temp_z1)
    #print(temp_z1.shape)
    #print(temp_z2.shape)
    #print(hh[0])
    #print('Temperature from P/d',np.amin(temp),np.amax(temp)) # Temperature from P/d
    #print(np.amax(temp)-np.amin(temp))
    #print('Temperature from P/d', np.amax(temp), np.amin(temp))
    #print (np.amax(d),np.amin(d))
    #print(tem)
    #print(tem.shape)
    #hh = np.amax(tem,axis=0)
    #print(hh)
    #print('Temperature from Data',np.amax(tem,axis=1),np.amin(tem)) # Temperature from P/d

    #cold = tem[tem < 10000]
    #print(cold)

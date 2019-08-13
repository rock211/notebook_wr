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
simid = 'n1e-4_v1414'

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
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
#print (unit['mass'], unit['time'], unit['magnetic_field'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

# In[10]:

# import pyathena as pa

Msun = unit['mass'].to('Msun').value
Myr=unit['time'].to('Myr').value

step = 50
age_cut = 10.
Age = []
M_u = []
M_l = []
M_t = []

M_U1 = []
M_U2 = []
M_U3 = []
M_L1 = []
M_L2 = []
M_L3 = []

M_T = []

H_U =[]
H_L =[]
#Age = pd.DataFrame()
print(Myr)
plt.figure()

for tidx in range(251, 501):
    vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir,simid, simid, tidx)

    # read in domain information
    ds = pa.AthenaDataSet(vtkfname)

    # full domain information
    ds.domain

    # information of grid #0
    ds.grids[0]

    # yet, we didn't read data.
    # let's read each data field in a full domain

    # this can be original data fields
    d = ds.read_all_data('density')
    nd = ds.read_all_data('number_density')
    tem = ds.read_all_data('temperature')

    starfname = vtkfname.replace('id0', 'idstarpar').replace('vtk', 'starpar.vtk')
    sp = pa.read_starvtk(starfname)
    #print(sp.shape)
    # print (sp)

    # SFR comparison between z > 0 & z < 0
    print(tidx)
    star_clu = sp[sp['mass'] != 0]  # SFR comparison between z > 0 & z < 0
    #print(star_clu['x1'])
    M_t = sum(star_clu['mass']) * Msun
    star_upper = star_clu[star_clu['x1'] > 0.]
    star_upper100 = star_upper[star_upper['age'] * Myr < 100]
    star_upper40 = star_upper[star_upper['age'] * Myr < 40]
    star_upper10 = star_upper[star_upper['age'] * Myr < 10]
    M_u1 = sum(star_upper100['mass']) * Msun
    M_u2 = sum(star_upper40['mass']) * Msun
    M_u3 = sum(star_upper10['mass']) * Msun

    #print('UPPER',star_upper['age'])
    star_lower = star_clu[star_clu['x1'] < 0.]
    star_lower100 = star_lower[star_lower['age'] * Myr < 100]
    star_lower40 = star_lower[star_lower['age'] * Myr < 40]
    star_lower10 = star_lower[star_lower['age'] * Myr < 10]
    M_l1 = sum(star_lower100['mass']) * Msun
    M_l2 = sum(star_lower40['mass']) * Msun
    M_l3 = sum(star_lower10['mass']) * Msun

    H_u = np.mean(star_upper100['x1'])
    H_l = np.mean(star_lower100['x1'])

    M_U1.append(M_u1 / (100000*100))
    M_U2.append(M_u2 / (100000*40))
    M_U3.append(M_u3 / (100000*10))
    M_L1.append(M_l1/(100000*100))
    M_L2.append(M_l2 / (100000*40))
    M_L3.append(M_l3 / (100000*10))

    H_U.append(H_u)
    H_L.append(abs(H_l))
    #M_T.append(M_t/100000)
    #print(star_clu['age'])
    #print(star_upper['age'],star_lower['age'])

    #Age.append(star_clu['age'])
    #Age.join(star_clu['age'],how='outer')
    #print('AGE',Age)


#Age.plot.hist(stacked=True, bins=20)
M1 = np.add(M_U1,M_L1)
M2 = np.add(M_U2,M_L2)
M3 = np.add(M_U3,M_L3)
#bins = np.arange(0, 500, 20) # fixed bin size
#plt.hist(Age, bins=bins, alpha=0.5)
#print('0',Age[0],'1',Age[1])
#np.savetxt('age.txt',Age)

xx = range(251,501)

plt.plot(xx,M1,label='100Myr',alpha=0.5)
plt.plot(xx,M2,label='40Myr',alpha=0.5)
plt.plot(xx,M3,label='10Myr',alpha=0.5)
plt.xlabel('time(Myr)')
plt.ylabel('SFR(Msun/yr)')
plt.title('SFR variation')
plt.legend(loc=0)
plt.show()

H_L = np.array(H_L)
H_L[np.isnan(H_L)] = 0
#print(H_L)
LM = np.mean(H_L)
UM = np.mean(H_U)
#print(LM,UM)
plt.plot(xx,H_U,label='Upper',color='navy',alpha = 0.5)
plt.axhline(y=UM, label = 'Upper mean',linestyle='--',color='navy')
plt.plot(xx,H_L,label='Lower',color='r', alpha = 0.5)
plt.axhline(y=LM, label = 'Lower mean',linestyle='--',color='r')


#plt.plot(xx,M_L1,label='100')
#plt.plot(xx,M_L2,label='40')
#plt.plot(xx,M_L3,label='10')
#plt.plot(xx,M_U1,label='100')
#plt.plot(xx,M_U2,label='40')
#plt.plot(xx,M_U3,label='10')


#plt.plot(xx,M_L,label='Lower')
#plt.plot(xx,M_T,label='Total',linestyle='--',alpha=0.5)
plt.xlabel('time')
#plt.ylabel('Msun/yr')
plt.ylabel('Mean Height')
plt.title('Mean Height from z=0 plane')
#plt.title('Upper Part')
#plt.title('age_cut = 10Myr')
#plt.xticks(tick_loc,tick)
plt.legend(loc=0)
plt.show()

'''
sub = sub[1:step]
#print(len(sub))

xxx= range(1,step)
#print(len(xxx))
plt.plot(xxx,sub,marker='o',linestyle='--',alpha=0.5)
plt.title('Mass diff')
plt.xlabel('time')
plt.ylabel('Mass difference')
plt.show()
'''
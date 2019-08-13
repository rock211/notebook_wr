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
simid = 'MHD_8pc_new' # MHD_8pc_new / RPS_8pc_n1e-4_v1414 / RPS_8pc_n2e-4_v1414

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
#print Myr
agebin = 10 # unit : Myr, 10 : H-alpha like, 40 : cluster lifetime

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
#print(Myr)
plt.figure()

for tidx in range(251, 501):
    vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir,simid, simid, tidx)

    # read in domain information
    #ds = pa.AthenaDataSet(vtkfname)

    # full domain information
    #ds.domain

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
    # print (sp)


    print(tidx)
    star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
    #print(star_clu['x1'])
    star_clu2 = star_clu[star_clu['age'] * Myr < agebin] # time scale cut
    M_star = sum(star_clu2['mass'])*Msun # mass sum in time scale

    #M_t = sum(star_clu['mass']) * Msun

    M_U1.append(M_star / (1000000*agebin)) # mass divided by time scale

#Age.plot.hist(stacked=True, bins=20)

#bins = np.arange(0, 500, 20) # fixed bin size
#plt.hist(Age, bins=bins, alpha=0.5)
#print('0',Age[0],'1',Age[1])
#np.savetxt('age.txt',Age)

xx = range(251,501)

plt.plot(xx,M_U1,label='%sMyr' % agebin,alpha=0.5)

plt.xlabel('time(Myr)')
plt.ylabel('SFR(Msun/yr)')
plt.title('SFR variation(nonICM)')
#plt.ylim(0,0.12)
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

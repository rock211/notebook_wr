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
M_Ut = []
M_Lt = []
step = 50
age_cut = 100.
Age = []
#Age = pd.DataFrame()
print(Myr)
plt.figure()
for q in range(0,step):
    M_U = []
    M_L = []
    for tidx in range(251+q*250/step, 251+(q+1)*250/step):

        vtkfname = '%s%s/id0/id0.%04d.vtk' % (basedir, simid, tidx)

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
        #print (d.shape)
        # print (d)
        # note that it adopts C-like indexing, k (z-index) comes first and i (x-index) comes last
        # vector field has 3 component
        v = ds.read_all_data('velocity')
        #print (v.shape)
        vx = v[:, :, :, 0]
        vy = v[:, :, :, 1]
        vz = v[:, :, :, 2]

        nd = ds.read_all_data('number_density')
        #print (nd.shape)

        tem = ds.read_all_data('temperature')
        #print (tem.shape)


        # # Star Particles

        # In[42]:

        # we have star particles, representing star clusters and runaway OB stars
        # the star particle information is stored at /data-directory/id0/
        # read_starvtk return particle information in "pandas.DataFrame", which provides a good data handler
        # please visit http://pandas.pydata.org/
        #print (vtkfname)
        starfname = vtkfname.replace('id0', 'idstarpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)
        # print (sp)


        # SFR comparison between z > 0 & z < 0

        star_clu = sp[sp['mass'] != 0]
        print(star_clu)
        star_upper = star_clu[star_clu['x3'] > 0.]
        star_upper = star_upper[star_upper['age'] * Myr < age_cut]
        M_u = sum(star_upper['mass']) * Msun
        #print('UPPER',star_upper['age'])
        star_lower = star_clu[star_clu['x3'] < 0.]
        star_lower = star_lower[star_lower['age'] * Myr < age_cut]
        M_l = sum(star_lower['mass']) * Msun
        #print(star_clu['age'])
        print(star_upper['age'],star_lower['age'])

        #Age.append(star_clu['age'])
        #Age.join(star_clu['age'],how='outer')
        #print('AGE',Age)
        M_U.append(M_u)
        M_L.append(M_l)

        #print(M_u,M_l)
    #print(len(Age))
    #frames = []
    print((sum(M_U)/(100000*250./(step)),sum(M_L)/(100000*250./(step))),q)
    M_Ut.append((sum(M_U)/(100000*250./(step))))
    M_Lt.append((sum(M_L)/(100000*250./(step))))

print(Age)
#Age.plot.hist(stacked=True, bins=20)

#bins = np.arange(0, 500, 20) # fixed bin size
#plt.hist(Age, bins=bins, alpha=0.5)
#plt.show()
#print(Age)
#print(Age.shape)
#Age0 = Age[0]
#print(Age0(0))
#print(Age0.shape)
#print('0',Age[0],'1',Age[1])
#np.savetxt('age.txt',Age)

xx = range(1,step+1)
Mean = np.add(M_Ut,M_Lt)/2
#print(Mean)
#print(Mean[49])
Mean2 = np.insert(Mean,0,0)
#print(Mean2)
Mean3 = np.append(Mean,0)
#print(Mean3)
sub = np.abs(Mean3 - Mean2)
#print(len(sub))

#tick_loc = range(0,step,step/10)
#tick = range(251,501,25)

#print(M_Ut[0],M_Ut[1])
plt.plot(xx,M_Ut,marker='o',label='Upper',linestyle='--',alpha=0.5)
plt.plot(xx,M_Lt,marker='o',label='Lower',linestyle='--',alpha=0.5)
plt.plot(xx,Mean,marker='o',label='Mean',linestyle='--',alpha=0.5)
plt.xlabel('time')
plt.ylabel('Msun/yr')
#plt.xticks(tick_loc,tick)
plt.legend(loc=0)
plt.show()

sub = sub[1:step]
#print(len(sub))

xxx= range(1,step)
#print(len(xxx))
plt.plot(xxx,sub,marker='o',linestyle='--',alpha=0.5)
plt.title('Mass diff')
plt.xlabel('time')
plt.ylabel('Mass difference')
plt.show()
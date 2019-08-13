import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import matplotlib.animation as animation

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
#print(unit)
print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])
MFu = 0.5476852239548456

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

# In[10]:

# import pyathena as pa
ims=[]
#fig, ax = plt.subplots()
for tidx in range(251, 501):

    vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)

    # read in domain information
    ds = pa.AthenaDataSet(vtkfname)

    # name of original data fields we stored from the simulation
    print(ds.field_list)

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

    # note that it adopts C-like indexing, k (z-index) comes first and i (x-index) comes last
    # vector field has 3 component

    nd = ds.read_all_data('number_density')
    #print (nd.shape)
    T = ds.read_all_data('T1')
    #print(T)
    tem = ds.read_all_data('temperature') # Temperature from data
    #print (tem.shape)

    # calculate sound speed
    P = ds.read_all_data('pressure')
    cs = np.sqrt(P / d)

    # calculation of temperature needs additional information about mean molecular weight I adopt
    coolftn = pa.coolftn()
    temp = coolftn.get_temp(P / d) # Temperature from P/d
    mag = ds.read_all_data('magnetic_field')
    mag1 = ds.read_all_data('magnetic_field1') #z axis comp.
    mag2 = ds.read_all_data('magnetic_field2') #y axis comp.
    mag3 = ds.read_all_data('magnetic_field3') #x axis comp.
    #print(mag.shape)
    mag_str = np.sqrt(mag1**2+mag2**2+mag3**2) # Sum of mag comp.
    #print(mag_str)
    #print(mag_str.shape)
    print(tidx)
    #print(np.amax(mag_str)*MFu,np.amin(mag_str)*MFu)

    #fig = plt.figure()

    #mag_up = mag_str[16:32]
    #mag_up = np.sum(mag_up,axis=0)
    #mag_lw = mag_str[0:16]
    #mag_lw = np.sum(mag_lw,axis=0)
    zlen = 448
    mag_str = np.sum(mag_str,axis = 0)
    d = np.sum(d,axis = 0)
    #print(mag_up.shape)
    #print(mag_lw.shape)
    plt.imshow(mag_str,origin='lower',animated=True,cmap = 'copper',interpolation='bilinear')
    plt.title('Projection Magnetic map_T = %s Myr' % tidx)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'B[$\mu$G]')
    plt.clim(0,1000)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('./magmap/MagMap_%s.png' % tidx)
    plt.clf()
    #plt.show()
    #ims.append([im])
    #plt.show()

    plt.imshow(d,origin='lower',interpolation='bilinear')
    cbar=plt.colorbar()
    cbar.ax.set_ylabel(r'$n_H[cm^{-3}]$')
    plt.clim(0,250)
    plt.title('Projection Density map_T = %s Myr' % tidx)
    #plt.tight_layout()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('./denmap/DenMap_%s.png' % tidx)
    plt.clf()
    #plt.show()

#ani = animation.FuncAnimation(fig, ims, interval = 200,blit=True)
#plt.draw()
#plt.show()
#ani.save('densitimap_z_%s.gif' %z, writer='imagemagick')

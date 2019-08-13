import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa

from mpl_toolkits import axes_grid1


from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl


# In[8]:

basedir = 'D:/yeongu/'
simid = 'MHD_8pc_new' # 'MHD_8pc_new' , 'RPS_8pc_n1e-4_v1414' , 'RPS_8pc_n2e-4_v1414'

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
print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
print unit['density']
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K
volpercell = 7168.*1024*1024/(128*128*896)
vpc = volpercell

# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

Mass = []
Mom_up = []
Mom_p = []
Mom_n = []
Mom_1 = []
Mom_2 = []
Mom_3 = []

stop = 501

for tidx in range(251, stop):  # time step 251, 331, 411, 501

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
    #comp = ds.read_all_data('specific_scalar0') # ism = 0 / icm = 1
    d = ds.read_all_data('density')*unit['density'].value # density
    mass = d # density times volume per cell = mass per cell
    #pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
    vel = ds.read_all_data('velocity')
    vel_z = vel[:,:,:,2]
    #print len(vel_z[vel_z >= 0])
    #print len(vel_z[vel_z < 0])
    #v_z_p = vel_z[vel_z > 0]
    #m_p = mass[vel_z> 0]
    #v_z_n = vel_z[vel_z < 0]
    #m_n = mass[vel_z < 0]

    ###### Mass #######


    ###### momentum #######
    moment_z = mass*vel_z/unit['time'].value #SolMass/pc2/Myr
    #print moment_z[895].shape
    m_up = np.sum(np.sum(moment_z[895]))*512/(8*1024*1024) # if sum, 1kpc^2 occur, so have to devide 1kpc^2 # Momentum at Upper bound
    #m_1kpc = np.sum(np.sum(moment_z[572]))*512/(8*1024*1024) # Momentum at 1kpc from plane
    #m_2kpc = np.sum(np.sum(moment_z[697]))*512/(8*1024*1024) # Momentum at 2kpc from plane
    #m_3kpc = np.sum(np.sum(moment_z[822]))*512/(8*1024*1024) # Momentum at 3kpc from plane

    #mom_p = np.sum(m_up[m_up > 0])
    #mom_n = np.sum(m_up[m_up < 0])

    #m_z_p = m_p*v_z_p
    #m_z_p = moment_z
    #m_z_n = m_n*v_z_n
    #mom_z = np.sum(np.sum(moment_z,axis=1),axis=1) # sum x & y axis / momentum function through z axis
    #mom_p = np.sum(np.sum(m_z_p,axis=1),axis=1)
    #mom_n = np.sum(np.sum(m_z_n, axis=1), axis=1)
    #print mom_z.shape

    Mom_up.append(m_up) ; #Mom_1.append(m_1kpc) ; Mom_2.append(m_2kpc) ; Mom_3.append(m_3kpc) ;

    #Mom_p.append(mom_p)
    #Mom_n.append(mom_n)

    print tidx
    #print len(vel_z)
    #print vel_z.shape
    #nd = ds.read_all_data('number_density') # number density
    #tem = ds.read_all_data('temperature')  # Temperature from data directly
    #coolftn = pa.coolftn()
    #temp = coolftn.get_temp(pre / d)  # Temperature derived from P/d & cooling func
    #magF = ds.read_all_data('magnetic_field')*unit['magnetic_field'].value # magnetic field
    #magP = ds.read_all_data('magnetic_pressure') # magnetic pressure

time = range(251,stop)
plt.plot(time,Mom_up,label='UpperBound') ; #plt.plot(time,Mom_3,label='3kpc') ;  plt.plot(time,Mom_2,label='2kpc'), plt.plot(time,Mom_1,label='1kpc')
plt.title('Mass Flux_nonICM')
plt.xlabel(r'Time [Myr]')
plt.ylabel(r'MassFlux_z $[M_{\odot} pc^{-2}Myr^{-1}]$')
plt.legend(loc=0)
plt.show()

'''
plt.plot(time,Mom_p)
plt.title('Momentum_z_n1_positive')
plt.xlabel(r'Time [Myr]')
plt.ylabel(r'Momentum_z')
plt.show()

plt.plot(time,Mom_n)
plt.title('Momentum_z_n1_negative')
plt.xlabel(r'Time [Myr]')
plt.ylabel(r'Momentum_z')
plt.show()
'''

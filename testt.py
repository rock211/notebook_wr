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
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
#print unit['density']
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K
vpc = 7168.*1024*1024/(128*128*896) # volume per cell


# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

stop = 501

simid_t = ('RPS_8pc_n1e-4_v1414' , 'RPS_8pc_n2e-4_v1414') #'MHD_8pc_new',
labell = ('ICM1','ICM2') #'NonICM',
C = ('dimgrey','coral','royalblue')
S = ('-.','--','-')

z = np.arange(0,7168,1)
g = 4.5181 * 1e-30 # gravitational constant : pc3/solmass*s2
sig_star = 42 # solmass/pc2
z_s = 245 # pc
r_0 = 8000 # pc
rho_dm = 0.0064 # solmass/pc3
g_z = 2.*np.pi*g*sig_star*(z-3584)/((z-3584)**2+z_s**2)**(0.5) + 4.*np.pi*g*rho_dm*((z-3584)/(1+(z-3584)**2/r_0**2)) #
#print g_z
meter = 3.24078*1e-17 # pc
kg = 5.02785*1e-31 # solar mass
# overplot Starformation rate of three different simulations

for j in range(len(simid_t)) :
    basedir = 'D:/yeongu/'
    simid = simid_t[j]
    Mom_up = []

    for tidx in range(281,501):  # time step 251, 331, 411, 501

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
        #gp = ds.read_all_data('gravitational_potential')
        #print unit['gravitational_potential']
        #print gp
        #print gp
        d = ds.read_all_data('density')*unit['density'].value # density : solMass/pc3
        #print d.shape
        d = np.repeat(d,8,axis=0)
        #print d.shape

        restP = d*g_z[:,None,None]
        restP = np.mean(np.mean(restP,axis=1),axis=1)*meter*10/kg/kb

        P = []
        for zz in range(7168) :
            #print len(restP[zz:7196])
            P.append(np.sum(restP[zz:7168]))


        #plt.plot(z,restP)
        plt.semilogy(z, P, label='restoring P', c='firebrick')
        plt.axvline(3584, linestyle='--', c='k', linewidth=0.5)
        plt.xticks([73 * 8, 198 * 8, 323 * 8, 448 * 8, 573 * 8, 698 * 8, 823 * 8],
                   ['-3', '-2', '-1', '0', '1', '2', '3'])
        plt.xlabel('z_axis[kpc]')
        plt.ylabel(r'$Pressure[K/cm^{-3}]$')
        #plt.ylim(1e2,1e6)
        plt.legend(loc=0)
        plt.title('Pressure Variance_%s_%s' % (labell[j],tidx))
        plt.show()
        #plt.savefig('D:/yeongu/plots/restoring/restoring_%s_%s.png' % (labell[j], tidx))
        plt.close()
        print tidx
        '''
        #plt.axhline(0, linestyle='--', c='k', linewidth=0.5)
        plt.axvline(3584, linestyle='--', c='k', linewidth=0.5)
        #plt.plot(g_z, z, c='purple')
        plt.semilogy(z,abs(restP),c='purple')
        # plt.plot(Phi,z)
        plt.xlabel('z_axis[kpc]')
        plt.xticks([73 * 8, 198 * 8, 323 * 8, 448 * 8, 573 * 8, 698 * 8, 823 * 8],
                   ['-3', '-2', '-1', '0', '1', '2', '3'])
        #plt.xlabel(r'gravity $[pc/s^2]$')
        plt.ylabel(r'Restoring Pressure $[Kcm^{-3}]$')
        plt.ylim(1e-2,1e4)
        plt.title('Restoring Pressure_%s_%s' % (labell[j],tidx))
        #plt.title('External Gravity')
        plt.savefig('D:/yeongu/plots/restoring/restoring_%s_%s.png' % (labell[j], tidx))
        #plt.show()
        plt.close()
        print tidx
        '''


#Phi = 2.*np.pi*g*sig_star*z_s*(((1.+((z-3584)**2)/(z_s**2))**(0.5))-1)+2.*np.pi*g*rho_dm*(r_0**2)*np.log(1+((z-3584)**2)/r_0**2)




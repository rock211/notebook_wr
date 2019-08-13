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

# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

start = 250
stop = 500

plt.figure(figsize=(10.5,5)) # original 5.5,5
simid_t = ('R8_8pc_metal','RPS_8pc_ICM1','RPS_8pc_ICM2') #
#simid_t = ('MHD_8pc_new' ,'RPS_8pc_n1e-4_v1414','RPS_8pc_n2e-4_v1414')
#simid_t = ('n1e-4_v1414','n2e-4_v1414')

labell = (r'No ICM','ICM1','ICM2') # r'No ICM',
labelll = ('Cold','Unstable','Warm','Ionized','Hot')
C = ('k','forestgreen','tomato') #('darkkhaki','royalblue','firebrick')
C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')
C3 = ('darkred','salmon','red')
S = (':','-','--')
alpha = (0.8,1,1)

# overplot Starformation rate of three different simulations

for j in range(3) :
    basedir = 'D:/yeongu/'
    simid = simid_t[j]
    Mom_up = []

    Mass_tot = []
    Mass_sub = []
    for tidx in range(start, stop):  # time step 251, 331, 411, 501

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
        #scalar = ds.read_all_data('specific_scalar4') # ism = 0 / icm = 1
        #gp = ds.read_all_data('gravitational_potential')
        #print unit['gravitational_potential']
        #print gp
        d = ds.read_all_data('density')*unit['density'].value # density
        '''
        nd = ds.read_all_data('number_density')*unit['number_density'].value
        icmm = d[0]
        icm_n = nd[0]
        print icmm.shape
        nd_icm = np.sum(np.sum(icm_n))
        d_icm = np.sum(np.sum(icmm))
        print nd_icm*8*8/1024/1024#*unit['time'].value
        print d_icm*8*8/1024/1024*unit['time'].value*1414
        '''
        #den = ds.read_all_data('density')
        #P = ds.read_all_data('pressure')
        #mass = d # density times volume per cell = mass per cell
        #pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
        #vel = ds.read_all_data('velocity')
        #T1 = ds.read_all_data('T1')
        #vel_z = vel[:,:,:,2]
        #coolftn = pa.coolftn()
        #temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func
        if j == 0:
            icm_inflow = 0
        elif j == 1:
            icm_inflow = (tidx-250)*unit['time'].value*0.00051002
        else :
            icm_inflow = (tidx-250)*unit['time'].value*0.00051002*2
        #print icm_inflow
        total_mass = d * 8*8*8
        sub_mass = total_mass-icm_inflow

        surf_mass = total_mass/(1024*1024)
        sub_surf_mass = sub_mass/(1024*1024)

        Mass_tot.append(np.sum(np.sum(np.sum(surf_mass))))
        Mass_sub.append(np.sum(np.sum(np.sum(sub_surf_mass))))
        #print Mass_tot
        print tidx
        #print len(vel_z)

    ut = round(unit['time'].value,4)
    time = np.arange(start, stop)*ut
    #xnew = np.linspace(time.min(),time.max(),5000)
    #smooth = spline(time,Mom_up,xnew)
    #plt.plot(time, Mom_cuw, label='Cold/Warm', color=C2[0])
    #plt.plot(time, Mom_hi, label='Hot', color=C2[3])
    #plt.plot(time, Mom_icm, label='ICM Flux', color='blueviolet')
    #plt.plot(time, Mom_t,label='Total', color='dimgrey',ls=':')

    #plt.plot(time, Mass_tot, label='%s' % labell[j], color=C[j],ls='--') # Hot only
    #plt.plot(time, Mass_sub, color=C[j],ls='--',alpha=alpha[j])
    plt.plot(time, Mass_sub, label='%s' % labell[j], color=C[j],alpha=alpha[j])

#plt.title('Mass Flux_UpperBound')
ml = MultipleLocator(5)
plt.minorticks_on()
plt.axes().tick_params(which='minor', direction='in')
plt.axes().tick_params(which='major', direction='in',labelsize=23)

plt.xlabel(r'Time [Myr]',fontsize=25)
plt.xticks([250,300,350,400,450])
plt.ylabel(r'$\Sigma_{gas}$ [M$_{\odot}$ pc$^{-2}$]',fontsize=25)
plt.legend(loc=0)
plt.xlim(time.min(),time.max())
plt.ylim(0,10)
plt.tight_layout(pad=0.3)
plt.savefig('D:/yeongu/plots/new_data/gas_surf_density_evolve_ONLY_ISM_Trans.png',dpi=500,transparent=True)
plt.show()
plt.close()

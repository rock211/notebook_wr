import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
from mpl_toolkits import axes_grid1
import numpy.polynomial.polynomial as poly
from numpy import inf

from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl

basedir = 'D:/yeongu/'
simid = 'RPS_8pc_n1e-4_v1414'

# ## Unit system
#
# The unit system we choose for this simulation is
# * [length] = pc
# * [velocity] = km/s
# * [density] = 1.4271*m_h/cm^3

unit = pa.set_units(muH=1.4271)
print(unit)
print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K

# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

Total = np.zeros((250,76))
xaxis = np.linspace(0.5,8,76)
print Total.shape
for tidx in range(251, 501):  # time step

    vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
    # read in domain information
    ds = pa.AthenaDataSet(vtkfname)

    # name of original data fields we stored from the simulation
    #print(ds.field_list)

    # It also has predefined data fields can be calculated from the original data.
    #print(ds.derived_field_list)

    # this can be original data fields
    comp = ds.read_all_data('specific_scalar0') # ism = 0 / icm = 1
    d = ds.read_all_data('density')*unit['density'].value # density
    #pre = ds.read_all_data('pressure')*unit['pressure'].value # thermal pressure
    #nd = ds.read_all_data('number_density') # number density
    tem = ds.read_all_data('temperature')  # Temperature from data directly
    #coolftn = pa.coolftn()
    #temp = coolftn.get_temp(P / d)  # Temperature derived from P/d & cooling func
    #magF = ds.read_all_data('magnetic_field')*unit['magnetic_field'].value # magnetic field
    #magP = ds.read_all_data('magnetic_pressure') # magnetic pressure

    #N = tem.shape[0]*tem.shape[1]*tem.shape[2] # total value number

    #nd = np.log10(nd) # Ndensity reshape
    d = (d) # density reshape
    T = np.log10(tem) # temperature reshape
    #P = np.log10(pre) # Pressure reshape

    #print np.max(T), np.min(T), np.median(T)

    Bin = []
    for i in xaxis :
        #print i
        #print np.where((T >= (i-1)/10) & (T <= i/10))
        bins = np.sum(d[(T >= i) & (T < (i+0.1))]) # Temperature bin
        #print bins
        Bin.append(bins)
        #dd =
    #print np.max(Bin)
    Bin = np.log10(Bin / np.sum(Bin))
    #print Bin

    l = tidx - 251
    Total[l,:] = Bin
    #plt.plot(xaxis,Bin)
    #plt.show()
    #print Total
    print tidx

Min = []
Max = []
Median = []
for j in range(76) :
    Min.append(np.min(Total[:,j]))
    #Min[Min < -10] =  -9
    Max.append(np.max(Total[:,j]))
    Median.append(np.median(Total[:,j]))

Min[np.where(Min == -inf)]= -10 # please help................

#np.savetxt('min',Min)

print Min

plt.plot(xaxis,Max,linewidth=0.7,color='red',linestyle='--', label='Max')
plt.plot(xaxis,Median,linewidth=1.7, label='Median')
plt.plot(xaxis,Min,linewidth=0.7,color='lime',linestyle='--', label='Min')

plt.fill_between(xaxis,Max,Min,color='grey',alpha=0.5)
plt.title('MassW_PDFs_n1(250Myr~500Myr)')
plt.legend(loc=1)
plt.xlim(0,8)
plt.ylim(-7,0)
plt.xlabel(r'log T [K]')
plt.ylabel('Fraction of Gas')
plt.show()

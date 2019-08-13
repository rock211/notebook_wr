import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
from scipy.stats import lognorm

from mpl_toolkits import axes_grid1


from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle
import time
# In[7]:
import matplotlib as mpl


# In[8]:

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



simid_t = ('RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1','RPS_8pc_ICM2','RPS_8pc_ICM3','RPS_8pc_ICM4') # 'R8_8pc_metal'
labell = ('ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # 'NonICM',r'No ICM',
C = ('dimgrey','coral','royalblue')
S = ('-.','--','-')

# overplot Starformation rate of three different simulations
start_time = time.time()
for j in range(1,2) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j] # simid_t[j]

    if j == 4:
        stop = 472
    elif j == 5:
        stop = 401
    else:
        stop = 500
    print simid
    Mean = []
    Median = []
    Std = []
    Heightt = []
    for tidx in range(250,stop):  # time step 251, 331, 411, 501

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # name of original data fields we stored from the simulation
        #print(ds.field_list)

        # It also has predefined data fields can be calculated from the original data.
        #print(ds.derived_field_list)

        d = ds.read_all_data('density')*unit['density'].value # density : solMass/pc3
        vel = ds.read_all_data('velocity')
        velz = vel[:, :, :, 2]
        #velx = vel[:, :, :, 0]
        #vely = vel[:, :, :, 1]
        #vel_z = ds.read_all_data('velocity')[:, :, :, 2]
        d_z = np.sum(np.sum(d,axis=1),axis=1)
        d_max = np.argmax(d_z)
        print ('%s sec' %(time.time()-start_time))
        #T1 = ds.read_all_data('T1')
        #coolftn = pa.coolftn()
        #temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func
        #os.mkdir(basedir + '%s/vel_d/%s' % (simid,tidx))
        #print vel_z[438:458].shape
        height = d_max
        #n1, bins1, patches1 = plt.hist(np.ravel(velz[height-5:height+5]), bins=np.arange(velz.min(), velz.max() + 5, 5), log=True,label='total', facecolor='k', alpha=0.5,weights=np.ravel(d[height-5:height+5]))
        vel_sel = velz[height-5:height+5]

        hist, binedges = np.histogram(np.ravel(vel_sel), bins=np.arange(velz.min(), velz.max() + 5, 5),weights=np.ravel(d[height-5:height+5]))#,density=True)
        x_c = (binedges[1::]+binedges[0:-1])/2
        binwidth = binedges[1]-binedges[0]
        plt.bar(x_c, hist/np.sum(hist), align='center', linewidth=0.2, width=binwidth, alpha=0.7,log=1)

        Height = (height-448)*8./1000
        mean = np.mean(vel_sel)
        std = np.std(vel_sel)
        median = np.median(vel_sel)

        plt.text(-1450,0.5,r'Height(@ $\rho_{max}$) = %s kpc $\pm$ 0.04kpc' % Height)
        plt.text(-1450,0.25,'Mean = %s' % round(mean,2))
        plt.text(-1450, 0.14, 'Std. = %s' % round(std,2))
        plt.text(-1450, 0.07, 'Median = %s' % round(median,2))
        plt.title('Weighted_%s_%s' % (labell[j],tidx))
        plt.xlabel(r'Velocity [km/s]')
        plt.ylabel(r'Fraction')
        plt.ylim(1e-7,1)
        plt.xlim(-1500,1500)
        #print n1, np.sum(n1)
        #plt.savefig('D:/yeongu/plots/vel_d2/Vel_hist_%s_%s.png' % (labell[j], tidx), dpi=200)
        #plt.show()
        plt.close()
        print tidx

        Heightt.append(Height)
        Mean.append(mean)
        Std.append(std)
        Median.append(median)

    xx = range(250, stop)
    plt.suptitle('%s' % labell[j])
    plt.subplot(221)
    plt.plot(xx, Heightt)
    plt.title(r'Height(@ $\rho_{max}$)')
    plt.ylabel('kpc')

    plt.subplot(222)
    plt.plot(xx, Mean)
    plt.title('Velocity mean')
    plt.ylabel('km/s')

    plt.subplot(223)
    plt.plot(xx, Median)
    plt.title('Velocity median')
    plt.ylabel('km/s')

    plt.subplot(224)
    plt.plot(xx, Std)
    plt.title('std.')
    plt.ylabel('km/s')

    plt.tight_layout()
    plt.savefig('D:/yeongu/plots/vel_d2/%s.png' % labell[j], dpi=200)
    plt.show()

    '''
        for i in range(323,896) : # -1kpc ~ upper bound
            #print vel_z[i]
            zdata = pd.DataFrame({'velocity_z': [vel_z[i]],'density': [d[i]],'Temperature': [temp[i]]})
            zdata.to_pickle('F:/yeongu/%s/vel_d/%s/vel_distribution_%s_%s' % (simid,tidx,tidx, i))
            # xx = np.linspace(-2,2.5,250)
            #n, bins, patches = plt.hist(s_d,bins=xx,facecolor='blueviolet',alpha=0.7)#,log=True)
            #fit = lognorm.pdf(xx-1,shape,loc=loc,scale=scale)
            #print bins
            #area_hist = .0
            #for ii in range(n.size):
            #    area_hist += (bins[ii + 1] - bins[ii]) *n[ii]
    '''
    '''
        #plt.plot(xx-1,fit*area_hist,c='k')
        #print bins
        #plt.yscale('log')
        plt.xlabel(r'$log$ SurfaceDensity$[M_{\odot}/pc^2]$')
        plt.ylabel('Number')
        #plt.xlim(0,200)
        #plt.ylim(0.5,1e4)
        plt.title('Hist_%s_%s' % (labell[j],tidx))
        ylim = 500
        plt.ylim(0,ylim)
        #plt.xlim(0,10)
        #plt.text(7,5*1e5,r'scalar cut = 0.5')
        plt.text(np.amin(bins)*1.05,ylim*0.9,r'Mean = %s' % mean)
        plt.text(np.amin(bins) *1.05, ylim * 0.85, r'Median = %s' % medi)
        plt.text(np.amin(bins) *1.05, ylim * 0.8, r's_cut = 0.75')
        #plt.text(np.amax(bins) * 0.75, ylim* 0.8, r'bin width = 1')
        plt.savefig('D:/yeongu/plots/surfD_hist/Hist_%s_%s.png' % (labell[j], tidx),dpi=400)

        #plt.show()
        plt.close()
        print tidx
    '''
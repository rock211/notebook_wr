import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
import matplotlib.colors as colors

from mpl_toolkits import axes_grid1


from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:
import matplotlib as mpl


# In[8]:

basedir = 'D:/yeongu/'
#simid = 'MHD_8pc_new' # 'MHD_8pc_new' , 'RPS_8pc_n1e-4_v1414' , 'RPS_8pc_n2e-4_v1414'

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


#simid_t = ('RPS_8pc_noICM_newacc', 'RPS_8pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc')  # 'MHD_8pc_new' ,
#resol='8pc'
simid_t = ('RPS_8pc_noICM_newacc', 'RPS_4pc_ICM1_newacc', 'RPS_4pc_ICM2_newacc')
resol = '4pc'
labell = ('No ICM', 'Weak', 'Strong', 'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
labelll = ('Cold','Unstable','Warm','Ionized','Hot')
C = ('goldenrod','royalblue','firebrick')
C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')
C3 = ('darkred','red','salmon')
S = (':','--','-')

# overplot Starformation rate of three different simulations

for j in (1,3) :
    basedir = 'G:/yeongu/'
    simid = simid_t[j]
    Mom_up = []
    print simid
    #if j==2:
    #    snap=(250,260,300,330,380,410)
    #    lsnap = len(snap)
    #else:
    #    snap = (250,260,300,320,390,450)
    #    lsnap = len(snap)
    #plt.figure(figsize=(lsnap, 5*2))

    cicm = plt.cm.Reds  # RdYlBu_r,Reds
    cicm._init()
    x = np.arange(cicm.N)
    alphas = 0.4 * (np.tanh((x - 100) / 50.) + 1)
    # alphas = np.linspace(0.5, 0.5, cicm.N)
    cicm._lut[:-3, -1] = alphas
    cicm._lut[-3, -1] = alphas.min()
    cicm._lut[-2, -1] = alphas.max()
    if j==2:
        stop=473
    else:
        stop=499
    snap=range(486,stop)
    sl=128
    for k,tidx in enumerate(snap):#range(250, 499):  # time step 251, 331, 411, 501
        plt.figure(figsize=(4.6, 5))
        #plt.figure(figsize=(2.5,7))
        print k
        #surf = ('{}{}/surf/{}.{:04d}.surf.p'.format(basedir,simid,simid,tidx))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # yet, we didn't read data.
        # let's read each data field in a full domain

        # this can be original data fields

        #gp = ds.read_all_data('gravitational_potential')

        dd = ds.read_all_data('density')
        #d = ds.read_all_data('density')*unit['density'].value # density
        #print unit['']
        #pre =ds.read_all_data('pressure')*unit['pressure'].value/kb
        #mass = d # density times volume per cell = mass per cell
        T1 = ds.read_all_data('T1')
        coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func
        #nd = ds.read_all_data('number_density')*unit['number_density'].value
        #pre = pre*unit['pressure'].value/kb # thermal pressure
        #d_mks = d * 6.76745 * 1e-11  # kg/km3
        vel_z = ds.read_all_data('velocity')[:,:,:,2]

        dd = dd[:,sl,:]
        temp = temp[:, sl, :]
        vel_z = vel_z[:,sl,:]
        vel_z[temp>20000]=np.nan

        #dd = np.mean(dd,axis=1)
        #temp = np.mean(temp,axis=1)
        #vel_z = np.mean(vel_z,axis=1)

        #dd = np.sum(dd,axis=1)
        #temp = np.sum(temp,axis=1)
        #vel_z = np.sum(vel_z,axis=1)
        if j!=0:
            scalar = ds.read_all_data('specific_scalar3') # ism = 0 / icm = 1
            #scalar[dd < 1]=0
            #scalar = np.nanmean(np.where(scalar != 0, scalar, np.nan), axis=1)
            scalar = scalar[:,sl,:]
            #scalar = np.mean(scalar,axis=1) # 4pc map
        #vel_z[dd < 1]=0
        #temp[dd < 1]=0
        #dd[dd<1]=0

        #vel_z = np.divide(np.sum(vel_z*dd,axis=1),np.sum(dd,axis=1))
        #temp = np.divide(np.sum(temp * dd, axis=1), np.sum(dd, axis=1))
        #dd = np.nanmean(np.where(dd!=0,dd,np.nan),axis=1)

        #vel_z = np.nanmean(np.where(vel_z!=0,vel_z,np.nan),axis=1)
        #temp = np.nanmean(np.where(temp != 0, temp, np.nan), axis=1)
        #dd = np.nanmean(np.where(dd!=0,dd,np.nan),axis=1)

        plt.subplot(1,4,1)
        im1=plt.imshow(dd,origin='lower',cmap='RdYlBu_r',norm=LogNorm(),aspect='auto',vmin=1e-4,vmax=1e2)#,vmin=1e-1,vmax=1e2)# #cmap='ocean',vmin=1e-4,vmax=20)#
        #im1.set_extent()
        #if j!=0:
            #im4=plt.imshow(scalar,origin='lower',cmap=cicm,aspect='auto',vmin=0,vmax=1)

        if resol=='8pc':
            plt.text(20,800,'%s Myr' %tidx)
            plt.xticks([0,64,128],['-0.5','0','0.5'])
            plt.yticks([64,192,320,448,576,704,832],[])


        plt.ylabel(r'z [kpc]')
        plt.xlabel(r'x [kpc]')
        #print 'number',np.nanmax(dd), np.nanmin(dd)

        ####### 4pc #############
        if resol=='4pc':
            if j!=0:
                plt.text(20, 1650, '%s Myr' % tidx)
                plt.xticks([0,64*2,128*2],['-0.5','0','0.5'])
                plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2],[])
            else:
                plt.text(20, 800, '%s Myr' % tidx)
                plt.xticks([0,64,128],['-0.5','0','0.5'])
                plt.yticks([64,192,320,448,576,704,832],[])

        plt.subplot(1,4,2)
        im2=plt.imshow(scalar,origin='lower',cmap=cicm,aspect='auto',vmin=0,vmax=1)
        #print 'solar',d2.max(), d2.min()
        #plt.title(r'Vel_z [km/s]')
        if resol=='8pc':
            plt.xticks([0,64,128],['-0.5','0','0.5'])
            plt.yticks([64,192,320,448,576,704,832],[])

        ########## 4pc #############
        if resol == '4pc':
            if j!=0:
                plt.xticks([0,64*2,128*2],['-0.5','0','0.5'])
                plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2],[])
            else:
                plt.xticks([0, 64, 128], ['-0.5', '0', '0.5'])
                plt.yticks([64, 192, 320, 448, 576, 704, 832], [])

        #plt.ylabel(r'z [kpc]')
        plt.xlabel(r'x [kpc]')

        plt.subplot(1,4,3)
        im3=plt.imshow(vel_z,origin='lower',cmap='seismic',aspect='auto',vmin=-100,vmax=100) #,vmin=-200,vmax=200)# ,norm=colors.SymLogNorm(linthresh=100,vmin=-200,vmax=200)
        #print 'solar',d2.max(), d2.min()
        #plt.title(r'Vel_z [km/s]')
        if resol=='8pc':
            plt.xticks([0,64,128],['-0.5','0','0.5'])
            plt.yticks([64,192,320,448,576,704,832],[])

        ########## 4pc #############
        if resol == '4pc':
            if j!=0:
                plt.xticks([0,64*2,128*2],['-0.5','0','0.5'])
                plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2],[])
            else:
                plt.xticks([0, 64, 128], ['-0.5', '0', '0.5'])
                plt.yticks([64, 192, 320, 448, 576, 704, 832], [])

        #plt.ylabel(r'z [kpc]')
        plt.xlabel(r'x [kpc]')

        plt.subplot(1,4,4)
        im4=plt.imshow(temp,origin='lower',cmap='gnuplot2_r',norm=LogNorm(),aspect='auto',vmin=1e3,vmax=1e8)#,vmin=1e3,vmax=1e7) #copper 4pc, slice gnuplot2_r
        #print 'solar',d2.max(), d2.min()
        #plt.title(r'Vel_z [km/s]')
        if resol=='8pc':
            plt.xticks([0,64,128],['-0.5','0','0.5'])
            plt.yticks([64,192,320,448,576,704,832],[])

        ########## 4pc #############
        if resol=='4pc':
            if j!=0:
                plt.xticks([0,64*2,128*2],['-0.5','0','0.5'])
                plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2],[])
            else:
                plt.xticks([0, 64, 128], ['-0.5', '0', '0.5'])
                plt.yticks([64, 192, 320, 448, 576, 704, 832], [])


        #plt.ylabel(r'z [kpc]')
        plt.xlabel(r'x [kpc]')


        print tidx
        cbar_1=plt.gcf().add_axes([0.77, 0.78, 0.02, 0.17])
        cbar_2=plt.gcf().add_axes([0.77, 0.56, 0.02, 0.17])
        cbar_3=plt.gcf().add_axes([0.77, 0.34, 0.02, 0.17])
        cbar_4 = plt.gcf().add_axes([0.77, 0.12, 0.02, 0.17])

        cbar1 = plt.colorbar(im1,cax=cbar_1,label='Density [cm$^{-3}$]')
        cbar2 = plt.colorbar(im2, cax=cbar_2, label='Scalar')
        cbar3 = plt.colorbar(im3,cax=cbar_3,label=r'V$_z$ [km/s]')
        cbar4 = plt.colorbar(im4, cax=cbar_4, label=r'Temperature [K]')

        plt.subplots_adjust(bottom=0.1, top=0.99, wspace=.25,right=0.76,left=0.1,hspace=0.11)
        #plt.tight_layout()
        plt.savefig('D:/yeongu/plots/map_new/newmap_slice_y0_%s_%s.png' % (simid,tidx),dpi=200)
        plt.close()
        #plt.show()
    #plt.show()
    #plt.draw()
#plt.show()




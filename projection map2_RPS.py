import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
from pyathena import ath_hst

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

stop = 501
#simid_t = ('RPS_8pc_noICM_newacc', 'RPS_8pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc')  # 'MHD_8pc_new' ,
#resol = '8pc'
simid_t = ('RPS_8pc_noICM_newacc', 'RPS_4pc_ICM1_newacc', 'RPS_4pc_ICM2_newacc')
resol = '4pc'
labell = ('No ICM', 'Weak', 'Strong', 'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',

labelll = ('Cold','Unstable','Warm','Ionized','Hot')
C = ('goldenrod','royalblue','firebrick')
C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')
C3 = ('darkred','red','salmon')
S = (':','--','-')

# overplot Starformation rate of three different simulations

cicm = plt.cm.Reds  # RdYlBu_r,Reds
cicm._init()
x = np.arange(cicm.N)
alphas = 0.4 * (np.tanh((x - 100) / 50.) + 1)
# alphas = np.linspace(0.5, 0.5, cicm.N)
cicm._lut[:-3, -1] = alphas
cicm._lut[-3, -1] = alphas.min()
cicm._lut[-2, -1] = alphas.max()

for j in (0,1,2) :
    basedir = 'G:/yeongu/'
    simid = simid_t[j]
    Mom_up = []
    if j==2:
        snap=(250, 256, 270, 283, 312, 355, 374, 390, 408, 434) # 250, 256, 270, 283, 312, 355, 374, 390, 408, 434 / 250,260,300,330,380,410
        lsnap = len(snap)
    else:
        snap = (250, 256, 270, 283, 335, 374, 400, 421, 455, 486) # 250, 256, 270, 283, 335, 374, 400, 421, 455, 486 / 250,260,300,320,390,450
        lsnap = len(snap)
    plt.figure(figsize=(lsnap, 4.3*3))
    if j==0:
        sl = 64
    else:
        sl = 128

    hstfilename = basedir + simid + '/hst/' + simid + '.hst'
    sn = ath_hst.read_w_pandas(hstfilename.replace('.hst', '.sn'))

    for k,tidx in enumerate(snap):#range(250, 499):  # time step 251, 331, 411, 501
        #plt.figure(figsize=(2.5,7))
        print k
        #surf = ('{}{}/surf/{}.{:04d}.surf.p'.format(basedir,simid,simid,tidx))
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)
        #surff = pickle.load(open(surf,'rb'))
        #print surff

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

        #gp = ds.read_all_data('gravitational_potential')
        #print unit['gravitational_potential']
        #print gp
        #print gp
        #T = ds.read_all_data('temperature')
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

        if j!=0:
            scalar = ds.read_all_data('specific_scalar3') # ism = 0 / icm = 1
            #scalar[dd < 1]=0
            #scalar = np.nanmean(np.where(scalar != 0, scalar, np.nan), axis=1)
            scalar = scalar[:,sl,:]
            #scalar = np.mean(scalar,axis=1)

        #vel_z[dd < 1]=0
        #temp[dd < 1]=0
        #dd[dd<1]=0
        #vel_z = np.nanmean(np.where(vel_z!=0,vel_z,np.nan),axis=1)
        #temp = np.nanmean(np.where(temp != 0, temp, np.nan), axis=1)
        #dd = np.nanmean(np.where(dd!=0,dd,np.nan),axis=1)


        #v_z_p = vel_z[vel_z > 0]
        #m_p = mass[vel_z> 0]
        #v_z_n = vel_z[vel_z < 0]
        #m_n = mass[vel_z < 0]

        #rp = abs(d_mks*gp)*1e-2/kb # restoring force : density x gravitational potantial(=surface density * acceleration)
        #ICM = d_mks*(vel_z**2)*1e-2/kb # Ram Pressure : density x velocity^2, 1e-2 : convert to cgs


        #print rp.shape
        #vz = np.mean(vel_z,axis=1)
        #rp = np.mean(rp,axis=1) # y projection
        #ICM = np.mean(ICM,axis=1) # y projection
        #pre = np.mean(pre,axis=1)
        #Tem = np.mean(temp,axis=1)
        #d1 = np.mean(d,axis=1)
        #dd = np.mean(dd, axis=1)
        #d2 = np.sum(d,axis=1)*8
        #dz = np.sum(d,axis=0)*8
        #print dz.max(),dz.min()
        #print np.log10(np.sum(dz))
        #print np.log10(np.sum(d))
        #print d2
        #nd1 = np.mean(nd,axis=1)
        #nd2 = np.sum(nd,axis=1)/8 * 1.4271
        #nd3 = np.sum(nd,axis=1)
        #print tidx
        #print 'mean max', nd1.max(), nd1.min()
        #print 'from density', d2.max()*6.7694164*1e-23/(1.4271*1.67*1e-24), r'cm^-3'
        #print 'sum/8 max', nd2.max()
        #print 'sum/2 max', nd3.max()/2.
        #print 'sum max', nd3.max()

        plt.subplot(3,lsnap,k+1)
        im1=plt.imshow(dd,origin='lower',cmap='GnBu',norm=LogNorm(),aspect='auto',vmin=1e-4,vmax=1e2) #cmap='ocean',vmin=1e-4,vmax=20)#RdYlBu_r
        if j!=0:
            im4=plt.imshow(scalar,origin='lower',cmap=cicm,aspect='auto',vmin=0,vmax=1)

        if resol=='8pc':
            plt.text(14, 840, '%s Myr' % tidx)
            plt.xticks([0,64,128],[])
            plt.yticks([64,192,320,448,576,704,832],[])


            if k==0:
                plt.ylabel(r'z [kpc]',fontsize=13)
                plt.yticks([64, 192, 320, 448, 576, 704, 832], ['-3', '-2', '-1', '0', '1', '2', '3'])
                #plt.xlabel(r'x [kpc]')

            plt.ylim(260, 896)
        print 'number',np.nanmax(dd), np.nanmin(dd)

        ####### 4pc #############
        if resol == '4pc':
            if j!=0:
                plt.title('%3d Myr' % (tidx*unit['time']).value)
                plt.xticks([0,64*2,128*2],[])
                plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2],[])

                if k == 0:
                    plt.ylabel(r'z [kpc]',fontsize=13)
                    plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2], ['-3', '-2', '-1', '0', '1', '2', '3'])

                plt.ylim(521, 1792)

            else:
                plt.title('%3d Myr' % (tidx*unit['time']).value)
                plt.xticks([0,64,128],[])
                plt.yticks([64,192,320,448,576,704,832],[])

                if k == 0:
                    plt.ylabel(r'z [kpc]',fontsize=13)
                    plt.yticks([64, 192, 320, 448, 576, 704, 832], ['-3', '-2', '-1', '0', '1', '2', '3'])

                plt.ylim(260, 896)


                # plt.xlabel(r'x [kpc]')
        plt.tick_params(which='both', direction='in')

        plt.subplot(3,lsnap,k+1+lsnap)
        im2=plt.imshow(vel_z,origin='lower',cmap='coolwarm',aspect='auto',vmin=-100,vmax=100) #,vmin=-200,vmax=200)#
        #print 'solar',d2.max(), d2.min()
        #plt.title(r'Vel_z [km/s]')
        if resol=='8pc':
            plt.xticks([0,64,128],[])
            plt.yticks([64,192,320,448,576,704,832],['-3','-2','-1','0','1','2','3'])
            plt.ylim(260, 896)

        ########## 4pc #############
        if resol == '4pc':
            if j!=0:
                plt.xticks([0,64*2,128*2],[])
                plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2],[])

                if k == 0:
                    plt.ylabel(r'z [kpc]',fontsize=13)
                    plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2], ['-3', '-2', '-1', '0', '1', '2', '3'])

                plt.ylim(521, 1792)

            else:
                plt.xticks([0, 64, 128], [])
                plt.yticks([64, 192, 320, 448, 576, 704, 832], [])

                if k == 0:
                    plt.ylabel(r'z [kpc]',fontsize=13)
                    plt.yticks([64, 192, 320, 448, 576, 704, 832], ['-3', '-2', '-1', '0', '1', '2', '3'])
                plt.ylim(260, 896)

        if k==0:
            plt.ylabel(r'z [kpc]',fontsize=13)
            #plt.xlabel(r'x [kpc]')
        plt.tick_params(which='both', direction='in')

        plt.subplot(3,lsnap,k+1+lsnap*2)
        im3=plt.imshow(temp,origin='lower',cmap='gnuplot2_r',norm=LogNorm(),aspect='auto',vmin=1e3,vmax=1e8)#gnuplot2_r / BuPu_r
        #print 'solar',d2.max(), d2.min()
        #plt.title(r'Vel_z [km/s]')
        if resol=='8pc':
            plt.xticks([0,64,128],['-0.5','0','0.5'])
            plt.yticks([64,192,320,448,576,704,832],[])
            if k==0:
                plt.ylabel(r'z [kpc]',fontsize=13)
                plt.yticks([64, 192, 320, 448, 576, 704, 832], ['-3', '-2', '-1', '0', '1', '2', '3'])
            plt.ylim(260, 896)
        ########## 4pc #############
        if resol == '4pc':
            if j!=0:
                plt.xticks([0,64*2,128*2],[])
                plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2],[])

                if k == 0:
                    plt.ylabel(r'z [kpc]',fontsize=13)
                    plt.yticks([64*2,192*2,320*2,448*2,576*2,704*2,832*2], ['-3', '-2', '-1', '0', '1', '2', '3'])
                    plt.xticks([0, 64 * 2, 128 * 2], ['-0.5', '0', '0.5'])
                plt.ylim(521, 1792)
            else:
                plt.xticks([0, 64, 128], [])
                plt.yticks([64, 192, 320, 448, 576, 704, 832], [])

                if k == 0:
                    plt.ylabel(r'z [kpc]',fontsize=13)
                    plt.yticks([64, 192, 320, 448, 576, 704, 832], ['-3', '-2', '-1', '0', '1', '2', '3'])
                    plt.xticks([0, 64, 128], ['-0.5', '0', '0.5'])
                plt.ylim(260, 896)

        if k==0:
            plt.ylabel(r'z [kpc]',fontsize=13)
            plt.xlabel(r'x [kpc]',fontsize=13)
        plt.tick_params(which='both', direction='in')

        print tidx
    cbar_1=plt.gcf().add_axes([0.9, 0.75, 0.015, 0.2])
    cbar_2=plt.gcf().add_axes([0.9, 0.41, 0.015, 0.2])
    cbar_3=plt.gcf().add_axes([0.9, 0.1, 0.015, 0.2])
    cbar1 = plt.colorbar(im1,cax=cbar_1)
    cbar1.set_label(label='Density [cm$^{-3}$]',size=12)
    cbar2 = plt.colorbar(im2,cax=cbar_2)
    cbar2.set_label(label=r'V$_z$ [km/s]',size=12)
    cbar3 = plt.colorbar(im3,cax=cbar_3)
    cbar3.set_label(label=r'Temperature [K]',size=12)


    plt.subplots_adjust(bottom=0.05, top=0.975, wspace=.1,hspace=0.06,right=0.89,left=0.07)
    #plt.show()
    #plt.draw()
    plt.savefig('D:/yeongu/plots/map_new/time_new3_temcut_y0_%s.png' % simid,dpi=200)
    plt.savefig('D:/yeongu/plots/map_new/time_new3_temcut_y0_%s.eps' % simid,format='eps')
#plt.show()
#plt.clim(1e1, 1e7)


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
print(unit)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
#print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

# In[10]:

# import pyathena as pa

idn = 'id10'
for tidx in range(251, 501):
    vtkfname = '%s%s/%s/%s.%04d.vtk' % (basedir, simid,idn,idn, tidx)
    fig = plt.subplot()
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
    # print(temp)
    #print(temp[np.argmin(temp)])
    #print(temp[np.argmax(temp)])
    #print((max(temp)[0]),(max(temp)[1]))
    #hh = np.a
    # min(temp,axis=0)
    #
    #print(temp[15])
    #print(tem[15])
    #print(temp[15]/tem[15])

    ##### Phase dividing * Calculate mass ##### Need exact value corresponds to tem. bound
    C=[]
    U=[]
    W=[]
    H=[]
    I=[]
    for z in range(32) :
        d_z = d[z]
        #print(np.amax(temp[z]),np.amin(temp[z]))
        #print(temp[z])
        #temp[z] = temp[z] - 58000.
        #print(temp[z])
        mc = np.sum(d_z[temp[z] < 184])
        mu = np.sum(d_z[(temp[z] > 184) & (temp[z] < 5050)])
        mw = np.sum(d_z[(temp[z] > 5050) & (temp[z] < 20000)])
        mh = np.sum(d_z[(temp[z] > 20000) & (temp[z] < 500000)])
        mi = np.sum(d_z[temp[z] > 500000])

        Mtot = np.sum(mc + mu + mw + mh + mi) #Total Mass at each z
        C.append(mc/Mtot)
        U.append(mu/Mtot)
        W.append(mw/Mtot)
        H.append(mh/Mtot)
        I.append(mi/Mtot)

        #print(z)

    z = range(32)
    plt.plot(z, C,label='C')
    plt.plot(z, U,label='U')
    plt.plot(z, W,label='W')
    plt.plot(z, H,label='H')
    plt.plot(z, I,label='I')
    plt.legend(loc=0)
    print(tidx)
    plt.title('Mass fraction_Phase_time = %s Myr' % tidx)
    plt.xlabel('z_direction')
    plt.ylabel('Mass_fraction')
    plt.savefig('./Tem_den/Tem_den_%s_%s' % (idn, tidx))
    plt.clf()
    #plt.show()
    #print(len(C),len(U),len(W),len(H),len(I))
    #print(Mtot_tidx)


    #print(temp_z.shape)
    #print(np.mean(temp[0]))
    #temp_z1 = temp_z[temp[0] > 58270]
    #temp_z2 = temp_z[temp[0] < np.mean(temp[0])]
    #print(temp_z1)
    #print(temp_z1.shape)
    #print(temp_z2.shape)
    #print(hh[0])
    #print('Temperature from P/d',np.amin(temp),np.amax(temp)) # Temperature from P/d
    #print(np.amax(temp)-np.amin(temp))
    #print('Temperature from P/d', np.amax(temp), np.amin(temp))
    #print (np.amax(d),np.amin(d))
    #print(tem)
    #print(tem.shape)
    #hh = np.amax(tem,axis=0)
    #print(hh)
    #print('Temperature from Data',np.amax(tem,axis=1),np.amin(tem)) # Temperature from P/d
    print('')
    #cold = tem[tem < 10000]
    #print(cold)









'''




    xmin = ds.domain['left_edge']
    xmax = ds.domain['right_edge']
    dx = ds.domain['dx']
    Nx = ds.domain['Nx']

    # set cell centered coordinates
    x = np.arange(xmin[0], xmax[0], dx[0]) + 0.5 * dx[0]
    y = np.arange(xmin[1], xmax[1], dx[1]) + 0.5 * dx[1]
    z = np.arange(xmin[2], xmax[2], dx[2]) + 0.5 * dx[2]

    # calculate background velocity vy_0=-q*Omega*x
    vy0 = -28.e-3 * x


    # find index for maximum density
    imax = np.argmax(d)
    # unravel the index to 3d-form
    imax = np.unravel_index(imax, d.shape)[::-1]
    # print (imax)

    # density projection at density maximum
    # Note that the z-axis comes first and x-axis comes last (C-style array indexing)
    dproj = []
    x_coord = ['x', 'x', 'y']
    i_coord = [0, 0, 1]
    y_coord = ['y', 'z', 'z']
    j_coord = [1, 2, 2]
    max_pos = [x[imax[0]], y[imax[1]], z[imax[2]]]


    # In[36]:

    # tempersure projection with density weights
    # Note that the z-axis comes first and x-axis comes last (C-style array indexing)
    Tslc=[]
    x_coord=['x','x','y']
    i_coord=[0,0,1]
    y_coord=['y','z','z']
    j_coord=[1,2,2]
    max_pos=[x[imax[0]],y[imax[1]],z[imax[2]]]

    fig=plt.figure(figsize=(15,15),dpi=300)

    Tslc.append(temp[imax[2],:,:])
    Tslc.append(temp[:,imax[1],:])
    Tslc.append(temp[:,:,imax[0]])
    #print(Tslc)
    for i in range(3):
        plt.subplot(1,3,i+1)
        im=plt.imshow(Tslc[i],origin='lower',norm=LogNorm())
        im.set_extent([xmin[i_coord[i]],xmax[i_coord[i]],xmin[j_coord[i]],xmax[j_coord[i]]])
        im.set_cmap(plt.cm.Spectral_r)
        #im.set_clim(10,1.e7)  # it makes plot weird(woorak)
        plt.axvline(max_pos[i_coord[i]],ls=':')
        plt.axhline(max_pos[j_coord[i]],ls=':')
        plt.xlabel(x_coord[i])
        plt.ylabel(y_coord[i])
    plt.tight_layout()
    plt.colorbar(im)
    plt.savefig('./temmap/tem_densityweight%s.png' % tidx)


    # # Star Particles

    # In[42]:

    # we have star particles, representing star clusters and runaway OB stars
    # the star particle information is stored at /data-directory/id0/
    # read_starvtk return particle information in "pandas.DataFrame", which provides a good data handler
    # please visit http://pandas.pydata.org/
    print (vtkfname)
    starfname = vtkfname.replace('id0', 'idstarpar').replace('vtk', 'starpar.vtk')
    sp = pa.read_starvtk(starfname)
    # print (sp)
    unit = pa.set_units(muH=1.4271)
    Msun = unit['mass'].to('Msun').value

'''

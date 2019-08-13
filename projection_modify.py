# In[6]:
from IPython import get_ipython

ipython_shell = get_ipython()
# get_ipython().magic(u'matplotlib inline')
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
simid = 'RPS_8pc_n1e-4_v1414'

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
print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
print (unit['mass'], unit['time'], unit['magnetic_field'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.

# In[10]:

# import pyathena as pa

for tidx in range(251, 501):
    vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)

    # read in domain information
    ds = pa.AthenaDataSet(vtkfname)

    # name of original data fields we stored from the simulation
    print(ds.field_list)

    # It also has predefined data fields can be calculated from the original data.
    print(ds.derived_field_list)

    # full domain information
    ds.domain
    print ds.domain

    # information of grid #0
    ds.grids[0]

    # yet, we didn't read data.
    # let's read each data field in a full domain

    # this can be original data fields
    d = ds.read_all_data('density')
    print (d.shape)
    # print (d)
    # note that it adopts C-like indexing, k (z-index) comes first and i (x-index) comes last
    # vector field has 3 component
    v = ds.read_all_data('velocity')
    print (v.shape)
    vx = v[:, :, :, 0]
    vy = v[:, :, :, 1]
    vz = v[:, :, :, 2]

    nd = ds.read_all_data('number_density')
    print (nd.shape)

    tem = ds.read_all_data('temperature')
    print (tem.shape)

    # or you can use derived fields
    v1 = ds.read_all_data('velocity1')
    v2 = ds.read_all_data('velocity2')
    v3 = ds.read_all_data('velocity3')

    print (vx.shape, v1.shape, (vx == v1).all())

    # calculate sound speed
    P = ds.read_all_data('pressure')
    cs = np.sqrt(P / d)

    # calculation of temperature needs additional information about mean molecular weight I adopt
    coolftn = pa.coolftn()
    temp = coolftn.get_temp(P / d)
    # print(temp)


    xmin = ds.domain['left_edge'] # -512x-512x-3584 (8 times larger than data cube
    xmax = ds.domain['right_edge'] # 512x512x3584
    #print xmin, xmax
    dx = ds.domain['dx']
    Nx = ds.domain['Nx']

    # set cell centered coordinates
    x = np.arange(xmin[0], xmax[0], dx[0]) + 0.5 * dx[0]
    y = np.arange(xmin[1], xmax[1], dx[1]) + 0.5 * dx[1]
    z = np.arange(xmin[2], xmax[2], dx[2]) + 0.5 * dx[2]

    # calculate background velocity vy_0=-q*Omega*x
    vy0 = -28.e-3 * x

    # substract it to get a perturbed velocity in y-dir.
    dv2 = v2 - np.tile(vy0.reshape(1, 1, Nx[0]), (Nx[2], Nx[1], 1))

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

    fig = plt.figure(figsize=(15, 15), dpi=300)
    '''
    for i in range(3):
        dproj.append(d.mean(axis=i))

    for i in range(3):
        plt.subplot(1,3,i+1)
        im=plt.imshow(dproj[i],origin='lower',norm=LogNorm())
        im.set_extent([xmin[i_coord[i]],xmax[i_coord[i]],xmin[j_coord[i]],xmax[j_coord[i]]])
        im.set_cmap(plt.cm.cubehelix)
        plt.axvline(max_pos[i_coord[i]],ls=':')
        plt.axhline(max_pos[j_coord[i]],ls=':')
        plt.xlabel(x_coord[i])
        plt.ylabel(y_coord[i])

    plt.tight_layout()
    plt.colorbar(im)
    plt.savefig('./densitymap/densitymap%s.png' % tidx)


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
    '''

    # # Star Particles

    # In[42]:

    # we have star particles, representing star clusters and runaway OB stars
    # the star particle information is stored at /data-directory/id0/
    # read_starvtk return particle information in "pandas.DataFrame", which provides a good data handler
    # please visit http://pandas.pydata.org/
    print (vtkfname)
    starfname=vtkfname.replace('id0','starpar').replace('vtk','starpar.vtk')
    sp=pa.read_starvtk(starfname)
    #print (sp)
    unit = pa.set_units(muH=1.4271)
    Msun = unit['mass'].to('Msun').value

    # In[43]:

    # Let's plot star particles
    plt.subplot(131,aspect='equal')
    plt.plot(sp['x1'],sp['x2'],'o',color = 'r', alpha = 0.3)
    plt.subplot(132,aspect='equal')
    plt.plot(sp['x1'],sp['x3'],'o',color = 'r', alpha = 0.3)
    plt.subplot(133, aspect='equal')
    plt.plot(sp['x2'],sp['x3'],'o',color = 'r', alpha = 0.3)
    plt.tight_layout()
    plt.show()

    # SFR comparison between z > 0 & z < 0
    #star = [sp['mass'] != 0]
    star_upper = sp[sp['x3'] > 0.]
    print(len(sp),len(star_upper))
    star_upper = star_upper[star_upper['mass'] != 0]
    
    print(sum(star_upper['mass'])*Msun)
    star_lower = sp[sp['x3'] < 0.]
    star_lower = star_lower[star_lower['mass'] != 0]
    print(len(star_lower),len(star_upper))

    plt.subplot(121,aspect='equal')
    plt.plot(star_lower['x1'],star_lower['x3'],'o',alpha = 0.3)
    plt.subplot(122, aspect='equal')
    plt.plot(star_upper['x1'], star_upper['x3'], 'o', alpha=0.3)
    plt.tight_layout()
    plt.title('star separate')
    plt.show()


    # select young star particles (age < 40Myr)
    # first, we need to convert age in Myr

    Myr=unit['time'].to('Myr').value

    #print (Myr)
    young_sp=sp[sp['age']*Myr < 40.]
    #print (len(sp),len(young_sp))

    # second, separate clusters (mass != 0) and runaways (mass == 0)
    # mass of runaway OB stars was set to zero
    runaway=young_sp[young_sp['mass'] == 0]
    cluster=young_sp[young_sp['mass'] != 0]
    plt.subplot(111,aspect='equal')
    plt.plot(cluster['x1'],cluster['x2'],'o')
    plt.plot(runaway['x1'],runaway['x2'],'.')
    print (len(runaway),len(cluster))
    plt.show()

    # In[54]:

    # using scatter plot, we can set size and color cluster particles based on their mass and age, respectively.
    # let's convert mass and age in Msun and Myr, respectively.
    Msun=unit['mass'].to('Msun').value
    mass=cluster['mass']*Msun
    age=cluster['age']*Myr
    ax=plt.subplot(111,aspect='equal')
    ax.scatter(cluster['x1'],cluster['x2'],marker='o',s=mass/100.,c=age,
               vmax=40,vmin=0,cmap=plt.cm.cool_r)
    plt.plot(runaway['x1'],runaway['x2'],'.k')
    plt.show()
    #plt.savefig('./2dstar/2d_starpar%s.png' % tidx)
    #from mpl_toolkits.mplot3d import Axes3D
    #ax.scatter(cluster['x1'],cluster['x2'],cluster['x3'],marker='o',s=mass/100.,c=age,vmax=40,vmin=0,cmap=plt.cm.cool_r)
    #ax.scatter(runaway['x1'],runaway['x2'],runaway['x3']'.k')


    # # Surface density map with star particles
    # Let's calculate surface density mape integrated along the z-axis:
    # $$\Sigma \equiv \int \rho dz $$

    # In[55]:

    dz=z[1]-z[0]
    surf=d.sum(axis=0)*unit['density']*dz*unit['length']
    #print (surf.shape)
    #print (surf.unit)
    #print (surf.mean())


    # In[56]:

    #it's always good idea to modulize your code
    def mass_norm(mass):

        #Mass normlization function to determine symbol size
        #This should be called both in sp_plot and sp_legend for the consistent result

        return np.sqrt(mass)
        #return mass/50.

    def sp_plot(ax,sp,plot_runaway=True):

        #This is comment for the function.
        #This will be printed when you call your function in help.
        #You may want to write usage or examples.

        import pyathena as pa
        unit=pa.set_units(muH=1.4271)
        Msun=unit['mass'].to('Msun').value
        Myr=unit['time'].to('Myr').value

        young_sp=sp[sp['age']*Myr < 40.]
        runaway=young_sp[young_sp['mass'] == 0]
        cluster=young_sp[young_sp['mass'] != 0]

        mass=cluster['mass']*Msun
        age=cluster['age']*Myr

        cl=ax.scatter(cluster['x1'],cluster['x2'],marker='o',s=mass_norm(mass),c=age,
                   vmax=40,vmin=0,cmap=plt.cm.cool_r)
        if(plot_runaway): ax.scatter(runaway['x1'],runaway['x2'],marker='.',color='k')

        return cl

    def sp_legend(ax,ref_mass=[1.e3,1.e4,1.e5]):
        ext=ax.images[0].get_extent()

        #plot particle references outside of the domain of interest
        s=[]
        label=[]
        for mass in ref_mass:
            s.append(ax.scatter(ext[1]*2,ext[3]*2,s=mass_norm(mass),color='k',alpha=.5))
            label.append(r'$10^%d M_\odot$' % np.log10(mass))
        ax.set_xlim(ext[0],ext[1])
        ax.set_ylim(ext[2],ext[3])
        legend=ax.legend(s,label,scatterpoints=1,loc=2,ncol=3,fontsize='small',bbox_to_anchor=(0.0, 1.1), frameon=False)

        return legend


    # In[57]:

    import matplotlib as mpl
    # let's combine surface density map with star particles.
    fig=plt.figure(figsize=(12,8))
    plt.rcParams['font.size']=18

    ax=fig.add_subplot(111)
    im=ax.imshow(surf.value,norm=LogNorm(),origin='lower')
    im.set_extent([xmin[0],xmax[0],xmin[1],xmax[1]])

    cl=sp_plot(ax,sp)
    leg=sp_legend(ax)
    im.set_clim(1.e-1,1.e2)
    im.set_cmap(plt.cm.pink_r)
    ax.set_xlabel(r'$x [{\rm pc}]$')
    ax.set_ylabel(r'$y [{\rm pc}]$')

    # Here's one example to draw colorbars
    cax2,kw=mpl.colorbar.make_axes(ax,fraction=0.1,pad=0.06)
    cbar2 = mpl.colorbar.ColorbarBase(cax2, ticks=[0,20,40],
                                      cmap=plt.cm.cool_r, norm=mpl.colors.Normalize(vmin=0,vmax=40),orientation='vertical')
    cbar2.set_label(r'${\rm age [Myr]}$')

    cax,kw=mpl.colorbar.make_axes(ax,fraction=0.1,pad=0.02)
    cbar=plt.colorbar(im,cax=cax)
    cbar.set_label(r'$\Sigma [M_{\odot} {\rm pc}^{-2}]$')

    # You can save figure by commenting out the following command
    plt.show()
    #fig.savefig('./surf/surfmap%s.png' % tidx,bbox_inches='tight')


    # In[58]:

    fig=plt.figure(figsize=(10,10))
    ax=plt.subplot(111)
    i=1
    im=plt.imshow(cs[:,imax[1],:],origin='lower',norm=LogNorm(vmin=1,vmax=1.e3),interpolation='nearest')
    im.set_extent([xmin[0],xmax[0],xmin[2],xmax[2]])
    im.set_cmap(plt.cm.Spectral_r)
    cb=plt.colorbar(im)
    cb.set_label(r'$c_s [{\rm km/s}]$')
    X, Y = np.meshgrid(x,y)
    ct=plt.contour(X,Y,d.mean(axis=0))

    #you can reuse the plotting functions defined before!
    sp_plot(ax,sp,plot_runaway=False)
    sp_legend(ax)

    plt.xlim(im.get_extent()[:2])
    plt.ylim(-1000,1000)
    plt.show()

    #fig.savefig('./contour/contour%s.png' % tidx,bbox_inches='tight')

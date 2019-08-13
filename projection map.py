import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa

from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

import matplotlib as mpl

basedir = 'D:/yeongu/'
simid = 'RPS_8pc_n2e-4_v1414' # 'MHD_8pc_new' , 'RPS_8pc_n1e-4_v1414' , 'RPS_8pc_n2e-4_v1414'

unit = pa.set_units(muH=1.4271)
print(unit)
print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K

# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

# ## Read Full data cube
#
# Original data cube is stored in "vtk" files. For the MPI simulations, Athena dumps the same number of vtk files with the number of proccessors. Each vtk file has data for all physical variables for a part of the simulation "domain" (it is called "grid"). I wrote a reader to read each grid, get data from it, and merge into a big data array for full domain.
tN = 25 # time bin
timestep1 = range(tN-1)
timestep2 = np.arange(251,501,tN)
length = len(timestep2)

fig, ax = plt.subplots(1,length,figsize=(length,7))
xticks = [[14, 39, 64, 89, 114], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
yticks = [[64, 192, 320, 448, 576, 704, 832], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
yticksname = [['-3', '-2', '-1', '0', '1', '2', '3'], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
xticksname = [['-0.4','-0.2', '0','0.2', '0.4'], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], []]
xlabel = [r'x [kpc]', '', '', '', '','' , '', '', '', '', '', '', '', '','' , '', '', '', '', '', '', '', '','' , '', '', '', '', '', '', '', '','' , '', '', '', '']
ylabel = [r'z [kpc]', '', '', '', '','' , '', '', '', '', '', '', '', '','' , '', '', '', '', '', '', '', '','' , '', '', '', '', '', '', '', '','' , '', '', '', '']


for i in range(length):  # time step 251, 331, 411, 501

    tidx = timestep2[i]
    vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
    # read in domain information
    ds = pa.AthenaDataSet(vtkfname)
    starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
    sp = pa.read_starvtk(starfname)
    #print ds.domain
    # name of original data fields we stored from the simulation
    #print(ds.field_list)

    # It also has predefined data fields can be calculated from the original data.
    #print(ds.derived_field_list)
    #print sp
    # full domain information
    #ds.domain

    # this can be original data fields
    #comp = ds.read_all_data('specific_scalar0') # ism = 0 / icm = 1
    d = ds.read_all_data('density') # density
    #pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
    P = ds.read_all_data('pressure') # original pressure
    nd = ds.read_all_data('number_density') # number density
    tem = ds.read_all_data('temperature')  # Temperature from data directly
    coolftn = pa.coolftn()
    temp = coolftn.get_temp(P / d)  # Temperature derived from P/d & cooling func


    #cut = 0.3
    # scalar cut data
    #d = np.log10(d[comp < cut])
    #P = np.log10(pre[comp < cut]/kb) # pressure(P/kb)
    #nd = np.log10(nd[comp < cut]) # number density
    #T = np.log10(tem[comp < cut]) # temperature

    nd_proj = np.mean(nd, axis=1)
    #nd_proj = nd[:,64,:] # y = 0 slice
    #tmep = np.mean(temp, axis=1)
    temp = temp[:,64,:]
    #temp = np.mean(temp, axis=1)
    #temp = tem[:,64,:]
    #scalar_proj = np.mean(comp, axis=1)
    #pre_proj = np.mean(pre , axis=1)
    #scalar_proj = np.sum(comp/len(comp[1]),axis=1)
    #os.mkdir('D:/yeongu/plots/8pc_n0_slice/%s_slice' % tidx)

###### projection image ######

    #plt.suptitle(r'Number density $[cm^{-3}]$', fontsize=13) # title change
    #plt.suptitle(r'Temperature $[K]$', fontsize=13)

    unit = pa.set_units(muH=1.4271)
    Msun = unit['mass'].to('Msun').value
    Myr = unit['time'].to('Myr').value
    young_sp = sp[sp['age'] * Myr < 40.] # choose young star
    star = young_sp[young_sp['mass'] != 0] # choose not runaway star
    mass=star['mass']*Msun
    age=star['age']*Myr
    star_x = star['x1']/8+64 # x-axis calibration
    star_z = star['x3']/8+448 # z-axis calibration

    time = round(tidx*unit['time'].value,1)

    im=ax[i].imshow(nd_proj, origin='lower',aspect='auto',norm=LogNorm(),cmap='RdYlBu_r')
    #im = ax[i].imshow(temp, origin='lower', aspect='auto', norm=LogNorm(), cmap='RdYlBu_r')

    ax[i].set_xlabel(xlabel[i])
    ax[i].set_xticks([14, 39, 64, 89, 114])
    ax[i].set_xticklabels(xticksname[i], rotation=-45, fontsize=8)
    ax[i].set_ylabel(ylabel[i])
    ax[i].set_yticks([64, 192, 320, 448, 576, 704, 832])
    ax[i].set_yticklabels(yticksname[i])
    ax[i].set_xlim(0,128)
    ax[i].set_title(r'%s Myr' % time, fontsize=8)
    #im.set_clim(5*1e-5, 5*1e+1) # number density slice colorlim
    im.set_clim(1e-5, 1e+1) # number density mean colorlim
    #im.set_clim(1e+2,1e+9) # Temperature(from data) colorlim
    #im.set_clim(1e+1, 1e+6) # Temperature(from function) colorlim(mean, slice both)

    im1 = ax[i].scatter(star_x, star_z,marker='o',lw=0.5,edgecolor='k',s=np.sqrt(mass)/8,c=age,vmax=40,vmin=0,cmap='coolwarm')
    #im2 = ax[i].imshow(scalar_proj, origin='lower',aspect='auto',norm=LogNorm(),cmap='binary',alpha=0.5)

    print tidx

    plt.subplots_adjust(wspace=0.01)

################# colorbar ##########################################

cbar_ax1 = fig.add_axes([0.91,0.575,0.015,0.3])
cbar1 = plt.colorbar(im,cax=cbar_ax1)
cbar1.ax.set_ylabel(r'$n_H$ [$cm^{-3}$]',fontsize=7) # Number density
#cbar1.ax.set_ylabel(r'$T$ [$K$]',fontsize=7) # Temperature
cbar1.ax.tick_params(labelsize=7)
cbar_ax2 = fig.add_axes([0.91,0.11,0.015,0.3])
cbar2 = plt.colorbar(im1,cax=cbar_ax2)
cbar2.ax.set_ylabel(r'Age [Myr]',fontsize=7)
cbar2.ax.tick_params(labelsize=7)

labell = [r'$10^3M_{\odot}$',r'$10^4M_{\odot}$',r'$10^5M_{\odot}$']
area = [np.sqrt(1000)/8,np.sqrt(10000)/8,np.sqrt(100000)/8]
for j in range(3):
    ax=plt.scatter([],[],c='k', s=area[j], label=labell[j])
plt.legend(loc='center left',bbox_to_anchor=(-0.5,1.3),scatterpoints=1, frameon=False,labelspacing=0.5,fontsize=8)

###################################################################################

#plt.savefig('D:/yeongu/plots/new/n2_bin%s_mean_nd.pdf' % tN,bbox_inches='tight',format='pdf')
plt.show()

'''
###### slice image #######
for m in range(128) :
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3) # originally 1 x 3 plot
    plt.suptitle(r'8pc_n0_%sMyr (y=%s)' % (tidx,m))
    nd_slice = nd[:,m,:] # choose y axis data
    pre_slice = pre[:,m,:]
    tem_slice = tem[:,m,:]
    #scalar_slice = comp[:,m,:]

    ax1.set_title(r'$N density [cm^{-3}]$',fontsize=10)
    ax1.set_xlabel('x[kpc]')
    ax1.set_ylabel('z[kpc]')
    ax1.set_xticks([0,64,128])
    ax1.set_xticklabels(['-0.5','0','0.5'])
    ax1.set_yticks([64,192,320,448,576,704,832])
    ax1.set_yticklabels(['-3','-2','-1','0','1','2','3'])
    im1 = ax1.imshow(np.log10(nd_slice), origin='lower',aspect='equal')
    divider1 = axes_grid1.make_axes_locatable(ax1)
    cax1 = divider1.append_axes('right',size = '20%',pad = 0.05)
    cbar1 = plt.colorbar(im1, cax = cax1)
    im1.set_clim(-3.5,2)

    ax2.set_title(r'$P/K_b [Kcm^{-3}]$',fontsize=10)
    ax2.set_xlabel('x[kpc]')
    #ax2.set_ylabel('z[kpc]')
    ax2.set_xticks([0,64,128])
    ax2.set_xticklabels(['-0.5','0','0.5'])
    ax2.set_yticks([64,192,320,448,576,704,832])
    #ax2.set_yticklabels(['-3','-2','-1','0','1','2','3'])
    ax2.set_yticklabels([])
    im2 = ax2.imshow(np.log10(pre_slice), origin='lower',cmap='copper',aspect='equal')
    divider2 = axes_grid1.make_axes_locatable(ax2)
    cax2 = divider2.append_axes('right',size = '20%',pad = 0.05)
    cbar2 = plt.colorbar(im2, cax = cax2)
    im2.set_clim(1, 6)

    #ax3.set_title(r'Scalar N',fontsize=10)
    ax3.set_title(r'Temperature [K]', fontsize=10)
    ax3.set_xlabel('x[kpc]')
    #ax3.set_ylabel('z[kpc]')
    ax3.set_xticks([0,64,128])
    ax3.set_xticklabels(['-0.5','0','0.5'])
    ax3.set_yticks([64,192,320,448,576,704,832])
    #ax3.set_yticklabels(['-3','-2','-1','0','1','2','3'])
    ax3.set_yticklabels([])
    im3 = ax3.imshow(np.log10(tem_slice), origin='lower',cmap='cool',aspect='equal') # originally scalar slice
    divider3 = axes_grid1.make_axes_locatable(ax3)
    cax3 = divider3.append_axes('right',size = '20%',pad = 0.05)
    cbar3 = plt.colorbar(im3, cax = cax3)
    im3.set_clim(2, 7.5)

    plt.subplots_adjust(wspace=0.1)
    plt.savefig('D:/yeongu/plots/8pc_n0_slice/%s_slice/n0_slice_%s_%s.png' % (tidx,tidx, m),bbox_inches='tight')

    #plt.show()
    plt.close()
print tidx
'''

'''
ax2.set_title(r'$P/K_b [Kcm^{-3}]$', fontsize=10)
ax2.set_xlabel('x[kpc]')
#ax2.set_ylabel('z[kpc]')
ax2.set_xticks([0,64,128])
ax2.set_xticklabels(['-0.5','0','0.5'])
ax2.set_yticks([64,192,320,448,576,704,832])
#ax2.set_yticklabels(['-3','-2','-1','0','1','2','3'])
ax2.set_yticklabels([])
im2 = ax2.imshow(np.log10(pre_proj), origin='lower',cmap='copper',aspect='equal')
ax2.scatter(star_x, star_z,marker='o',s=mass/5000.,c=age,vmax=40,vmin=0,cmap='autumn')
divider2 = axes_grid1.make_axes_locatable(ax2)
cax2 = divider2.append_axes('right',size = '20%',pad = 0.05)
cbar2 = plt.colorbar(im2, cax = cax2)
im2.set_clim(1, 6)

ax3.set_title(r'Scalar N', fontsize=10)
ax3.set_xlabel('x[kpc]')
#ax3.set_ylabel('z[kpc]')
ax3.set_xticks([0,64,128])
ax3.set_xticklabels(['-0.5','0','0.5'])
ax3.set_yticks([64,192,320,448,576,704,832])
#ax3.set_yticklabels(['-3','-2','-1','0','1','2','3'])
ax3.set_yticklabels([])
im3 = ax3.imshow(scalar_proj, origin='lower',cmap='cool',aspect='equal')
im4 = ax3.scatter(star_x, star_z, marker='o',s=mass/8000.,c=age,vmax=40,vmin=0,cmap='autumn')
divider3 = axes_grid1.make_axes_locatable(ax3)
cax3 = divider3.append_axes('right',size = '20%',pad = 0.05)
cbar3 = plt.colorbar(im3, cax = cax3)
im3.set_clim(0, 1)
cbar4 = plt.colorbar(im4,ticks=np.linspace(0,40,5))
cbar4.set_label(r'Age [Myr]')
'''



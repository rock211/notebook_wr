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
from astropy.io import fits

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

simid_t = ('RPS_8pc_noICM_newacc', 'RPS_8pc_ICM0_newacc', 'RPS_8pc_ICM1_newacc','RPS_4pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc','RPS_4pc_ICM2_newacc', 'RPS_8pc_ICM3_newacc')  # 'MHD_8pc_new' ,
#labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # r'No ICM',
labell = ('No ICM','P1', 'P3','P3h', 'P7','P7h','P14')
labelll = ('Cold','Unstable','Warm','Ionized','Hot')

C2 = ('darkblue','orange','goldenrod','red','firebrick')

S = (':','--','-')
hh = [0.0062,0.0062,0.0062,0.02,0.152,0.152,0.25]
ylim = [0.0062,0.0062,0.0062,0.02,0.152,0.152,0.25]
ylim_c = [0.15,0.15,0.15,0.35,10,10,10]
inflow = [0,0.005002/4,0.005002,0.005002,0.005002*2,0.005002*2,0.005002*4] # icm inflow rate
k = 0


icm_c = '#03d803ff'
jj = (1,6) # (0,2,4,1,6) #
#fig, axs = plt.subplots(len(jj),1, figsize=(10,12), sharex = True)

for j in jj: #range(1,7) :
    basedir = 'G:/yeongu/'
    #basedir='/media/woorak/data2/yeongu/' # for ubuntu
    simid = simid_t[j]
    Mom_up = []

    if j == 5 or j==6:
        stop = 473
    else:
        stop = 499

    Mom_1 = []
    Mom_2 = []
    Mom_3 = []
    Mom_c = []
    Mom_u = []
    Mom_w = []
    Mom_i = []
    Mom_h = []
    Mom_cuw = []
    Mom_hi = []
    Mom_t = []
    Mom_icm = []
    Mass = []
    #plt.figure(figsize=(8, 5)) # shrink 12,5

    dcut = 10  # average H2/H1 ratio is around 0.15 which is corresponding to observation box gas ratio / 4pc

    for tidx in range(250, stop):  # time step 251, 331, 411, 501

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # name of original data fields we stored from the simulation
        #print(ds.field_list)

        # It also has predefined data fields can be calculated from the original data.
        #print(ds.derived_field_list)

        #gp = ds.read_all_data('gravitational_potential')
        #print unit['gravitational_potential']
        #print gp
        if j==3 or j==5:

            d = ds.read_all_data('density')*unit['density'].value*4*4*4 # density in solarmass per pc3
        else:
            d = ds.read_all_data('density') * unit['density'].value * 8 * 8 * 8  # density in solarmass per pc3

        #print d.shape
        nd = ds.read_all_data('density')

        #d_cc = np.sum(d[nd >=dcut])
        #d_ww = np.sum(d[(nd <dcut)])# & (nd >0.001)])
        #ratio = d_cc/d_ww
        #print 'H2/H1 ratio', ratio

        #hdu = fits.PrimaryHDU(d)
        #hdu.writeto('%s_%s.fits' % (labell[j],tidx))
        #nd = ds.read_all_data('density')*8*8*8 # number density
        #den = ds.read_all_data('density')
        #P = ds.read_all_data('pressure')
        #mass = d # density times volume per cell = mass per cell
        #pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
        #vel = ds.read_all_data('velocity')
        T1 = ds.read_all_data('T1')
        #vel_z = vel[:,:,:,2]
        coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func
        #temp = ds.read_all_data('temperature')
        #v_z_p = vel_z[vel_z > 0]
        #m_p = mass[vel_z> 0]
        #v_z_n = vel_z[vel_z < 0]
        #m_n = mass[vel_z < 0]

        ###### momentum #######
        #moment_z = d*vel_z/unit['time'].value #SolMass/pc2/Myr
        #m_c = np.sum(np.sum(d_c * v_c / unit['time'].value))*512/(8*1024*1024)
        #m_u = np.sum(np.sum(d_u * v_u / unit['time'].value))*512/(8*1024*1024)
        #m_w = np.sum(np.sum(d_w * v_w / unit['time'].value))*512/(8*1024*1024)
        #m_h = np.sum(np.sum(d_h * v_h / unit['time'].value))*512/(8*1024*1024)
        #m_i = np.sum(np.sum(d_i * v_i / unit['time'].value))*512/(8*1024*1024)

        #m_cuw = np.sum(np.sum(d_cuw[v_cuw > 67.1] * v_cuw[v_cuw > 67.1] / unit['time'].value)) * 512 / (8 * 1024 * 1024) # velocity cut
        #m_hi = np.sum(np.sum(d_hi[v_hi > 67.1] * v_hi[v_hi > 67.1] / unit['time'].value)) * 512 / (8 * 1024 * 1024) # velocity cut

        #tot = m_cuw+m_hi
        #Mom_cuw.append(m_cuw); Mom_hi.append(m_hi)
        #d_c = d[temp < 184];# v_d_cuw = vel_z_down[temp_down < 20000]
        #d_u = d[(temp > 184) & (temp < 5050)];#  v_d_i = vel_z_down[(temp_down > 20000) & (temp_down < 500000)]
        #d_w = d[(temp > 5050) & (temp < 20000)];#  v_d_h = vel_z_down[temp_down > 500000]
        #print np.log10(np.sum(d))
        #d[temp > 20000] = 0
        #nd[temp > 20000]=0
        #d_c = d[temp < 200]
        #d_c2 = d[(temp > 200) & (temp < 5050)]
        #print np.log10(np.sum(d_c))
        #print np.log10(np.sum(d_c2))


        if j==3 or j==5:
            ### 4pc ###
            d_top = d[1046:1646,:,:]
            nd_top = nd[1046:1646, :, :]
            temp_top = temp[1046:1646,:,:]# 1046-1792 : 0.6~3 kpc

            d_bot = d[746:1046,:,:] # 746-1046 : -0.6~0.6 kpc
            nd_bot = nd[746:1046, :, :]
            temp_bot = temp[746:1046, :, :]
        else:
            ### 8pc ###
            d_top = d[523:823, :, :] # 523-896 : 0.6~3
            nd_top = nd[523:823, :, :]
            temp_top = temp[523:823, :, :] # 583-896 : 1~3.5kpc

            d_bot = d[373:523, :, :] # 373-523 : -0.6~0.6 kpc
            nd_bot = nd[373:523, :, :] # 271-583 : -1~1 kpc
            temp_bot = temp[373:523, :, :]

        #print 'whole',np.log10(np.sum(d))
        #print 'top',np.log10(np.sum(d_top))
        #print 'bot',np.log10(np.sum(d_bot))
        #print 'top_cold', np.log10(np.sum(d_top[temp_top<5050]))
        #print 'top_warm', np.log10(np.sum(d_top[temp_top>5050]))
        #print 'bot_cold', np.log10(np.sum(d_bot[temp_bot<5050]))
        #print 'bot_warm', np.log10(np.sum(d_bot[temp_bot>5050]))
        top = np.log10(np.sum(d_top))
        bot = np.log10(np.sum(d_bot))
        print 'top', top
        print 'bot', bot
        top_c = np.log10(np.sum(d_top[nd_top > dcut]))
        top_w = np.log10(np.sum(d_top[nd_top < dcut]))
        bot_c = np.log10(np.sum(d_bot[nd_bot > dcut]))
        bot_w = np.log10(np.sum(d_bot[nd_bot < dcut]))
        #top_c = np.log10(np.sum(d_top[temp_top < 5050]))
        #top_w = np.log10(np.sum(d_top[temp_top > 5050]))
        #bot_c = np.log10(np.sum(d_bot[temp_bot < 5050]))
        #bot_w = np.log10(np.sum(d_bot[temp_bot > 5050]))
        #
        # top_c = np.log10(np.sum(d_top[temp_top < 200]))
        #top_w = np.log10(np.sum(d_top[(temp_top > 200)& (temp_top < 5050)]))
        #bot_c = np.log10(np.sum(d_bot[temp_bot < 200]))
        #bot_w = np.log10(np.sum(d_bot[(temp_top > 200)& (temp_top < 5050)]))
        Mass.append([top,top_c,top_w,bot,bot_c,bot_w])
        
        # mass_proj=np.sum(d,axis=1)
        #plt.imshow(np.sum(d,axis=0)/8/8,origin='lower')
        #plt.subplot(1,2,1)
        #plt.imshow(mass_proj,origin='lower',norm=LogNorm(),cmap='bone_r')
        #plt.colorbar()
        #plt.subplot(1,2,2)
        #plt.imshow(mass_proj/8/8, origin='lower', norm=LogNorm(), cmap='bone_r')
        #plt.colorbar()
        #plt.show()
        #nd_c = nd[temp < 184];# v_d_cuw = vel_z_down[temp_down < 20000]
        #nd_u = nd[(temp > 184) & (temp < 5050)];#  v_d_i = vel_z_down[(temp_down > 20000) & (temp_down < 500000)]
        #nd_w = nd[(temp > 5050) & (temp < 20000)];#  v_d_h = vel_z_down[temp_down > 500000]
        #print np.sum(nd_c), np.sum(d_c)
        #print np.sum(nd_u), np.sum(d_u)
        #print np.sum(nd_w), np.sum(d_w)
        #print Mom_cuw
         #; Mom_t.append(tot)
        #Mom_c.append(m_c) ; Mom_u.append(m_u) ; Mom_w.append(m_w) ; Mom_h.append(m_h) ; Mom_i.append(m_i)
        #print moment_z[895].shape
        #m_up = np.sum(np.sum(moment_z[895]))*512/(8*1024*1024) # if sum, 1kpc^2 occur, so have to devide 1.024*1.024kpc^2 # Momentum at Upper bound
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

        #Mom_up.append(m_up) ; #Mom_1.append(m_1kpc) ; Mom_2.append(m_2kpc) ; Mom_3.append(m_3kpc) ;


        print tidx
    Mass = np.array(Mass)
    np.savetxt('Mass_%s_dcut_%s_newrange.txt' % (labell[j],dcut),Mass)
        #print len(vel_z)
    '''
    ##### data save #########
    Mom_cuw = np.array(Mom_cuw)
    Mom_i = np.array(Mom_i)
    Mom_h = np.array(Mom_h)
    np.savetxt('Mom_cuw_%s.txt' % labell[j],Mom_cuw)
    np.savetxt('Mom_i_%s.txt' % labell[j], Mom_i)
    np.savetxt('Mom_h_%s.txt' % labell[j], Mom_h)
    if j!=0:
        Mom_icm = np.array(Mom_icm)
        np.savetxt('Mom_icm_%s.txt' % labell[j], Mom_icm)
    #########################

    time = np.arange(250*unit['time'].value,stop*unit['time'].value,unit['time'].value)

    #### cumsum overlap ####
    ax1 = axs[k].twinx()
    ax1.plot(time,np.cumsum(Mom_cuw), label = '',c = C2[0], lw = 1, ls = '--')
    ax1.plot(time,np.cumsum(Mom_i), label='', c='#ef6fffff', lw=1, ls='--')
    ax1.plot(time,np.cumsum(Mom_h), label='', c=C2[3], lw=1, ls='--')
    ax1.set_ylim(-ylim_c[j]*0.02,ylim_c[j])
    ax1.tick_params(which='major',direction='in',labelsize=20)
    #if k == 1:
    #    ax1.set_ylabel('Cumulative mass [M$_{\odot}$ pc$^{-2}$]',fontsize=24)
    ########################

    #### Mass flux #########
    axs[k].plot(time, Mom_cuw, label='Cold/Warm', color=C2[0],lw=2.5)
    axs[k].plot(time, Mom_i, label='Ionized', color='#ef6fffff',lw=2.5)
    axs[k].plot(time, Mom_h, label='Hot', color=C2[3],lw=2.5)
    if j!=0:
        axs[k].axhline(inflow[j], label='ICM inflow', ls=':', color='k')
        axs[k].plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')

    axs[k].text(249,hh[j]*0.9,'%s' % labell[j],fontsize=20)
    axs[k].set_ylim(-ylim[j]*0.02,ylim[j])
    axs[k].tick_params(which='major', direction='in', labelsize=20)
    #########################
    plt.tight_layout()
    k = k+1

#axs[1].set_ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]',fontsize=24)
axs[k].set_xlabel('Time [Myr]',fontsize=20)
plt.xlim(250*unit['time'].value,499*unit['time'].value)

######### y label ###########
ax=fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
ax.tick_params(which='both',labelcolor='none', top=False, bottom=False, left=False, right=False)
ax.set_ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]',fontsize=24,labelpad=45)

ax2 = fig.add_subplot(111, sharex=ax, frameon=False)
#ax2.set_frame_on()
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.tick_params(which='both',labelcolor='none', top=False, bottom=False, left=False, right=False)
ax2.set_ylabel('Cumulative mass [M$_{\odot}$ pc$^{-2}$]',fontsize=24,labelpad=35)
#############################
plt.tight_layout()
#plt.savefig('D:/yeongu/plots/paperplot/new/massflux_No12_neww.png',dpi=300) # ,transparent=True
plt.show()
'''

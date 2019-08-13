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
jj = (5,10)
fig, axs = plt.subplots(len(jj),1, figsize=(10,12), sharex = True)

for j in jj: #range(1,7) :
    basedir = 'G:/yeongu/'
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
    #plt.figure(figsize=(8, 5)) # shrink 12,5
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
        d = ds.read_all_data('density')*unit['density'].value # density
        #den = ds.read_all_data('density')
        #P = ds.read_all_data('pressure')
        #mass = d # density times volume per cell = mass per cell
        #pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
        vel = ds.read_all_data('velocity')
        T1 = ds.read_all_data('T1')
        vel_z = vel[:,:,:,2]
        coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func
        #temp = ds.read_all_data('temperature')
        #v_z_p = vel_z[vel_z > 0]
        #m_p = mass[vel_z> 0]
        #v_z_n = vel_z[vel_z < 0]
        #m_n = mass[vel_z < 0]

        ###### Mass #######
        if j==3 or j==5:
            height=-1
        else:
            height = 895 #698 823 #895 #data set
        #print d.shape
        d_up = d[height]; d_down = d[0]
        #print d_up.shape
        vel_z_up = vel_z[height]; vel_z_down = vel_z[0]
        temp_up = temp[height]; temp_down = temp[0]

        #print np.median(temp_down), np.mean(temp[50])
        #d_c = d_up[temp_up < 184] ; v_c = vel_z_up[temp_up < 184]
        #d_u = d_up[(temp_up > 184) & (temp_up < 5050)] ; v_u = vel_z_up[(temp_up > 184) & (temp_up < 5050)]
        #d_w = d_up[(temp_up > 5050) & (temp_up < 20000)] ; v_w = vel_z_up[(temp_up > 5050) & (temp_up < 20000)]
        #d_h = d_up[(temp_up > 20000) & (temp_up < 500000)] ; v_h = vel_z_up[(temp_up > 20000) & (temp_up < 500000)]
        #d_i = d_up[temp_up > 500000] ; v_i = vel_z_up[temp_up > 500000]
        d_cuw = d_up[temp_up < 20000] ; v_cuw = vel_z_up[temp_up < 20000]
        d_i = d_up[(temp_up > 20000) & (temp_up < 500000)] ; v_i = vel_z_up[(temp_up > 20000) & (temp_up < 500000)]
        d_h = d_up[temp_up > 500000];   v_h = vel_z_up[temp_up > 500000]
        if j==3 or j==5:
            Area = 64. / (4 * 1024 * 1024)
        else:
            Area = 512. / (8 * 1024 * 1024)

        m_cuw = np.sum(np.sum(d_cuw * v_cuw / unit['time'].value)) * Area
        m_i = np.sum(np.sum(d_i * v_i / unit['time'].value)) * Area
        m_h = np.sum(np.sum(d_h * v_h / unit['time'].value)) * Area

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
        if j==0:

            d_d_cuw = d_down[temp_down < 20000]; v_d_cuw = vel_z_down[temp_down < 20000]
            d_d_i = d_down[(temp_down > 20000) & (temp_down < 500000)];  v_d_i = vel_z_down[(temp_down > 20000) & (temp_down < 500000)]
            d_d_h = d_down[temp_down > 500000];  v_d_h = vel_z_down[temp_down > 500000]

            m_d_cuw = np.sum(np.sum(d_d_cuw * v_d_cuw / unit['time'].value)) * Area
            m_d_i = np.sum(np.sum(d_d_i * v_d_i / unit['time'].value)) * Area
            m_d_h = np.sum(np.sum(d_d_h * v_d_h / unit['time'].value)) * Area
            #print m_d_cuw+m_d_hi
            Mom_cuw.append(m_cuw - m_d_cuw); Mom_i.append(m_i - m_d_i) ;Mom_h.append(m_h - m_d_h)
        else:

            scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1
            scalar_up = scalar[height]
            m_icm = np.sum(np.sum(scalar_up * d_up * vel_z_up / unit['time'].value)) * Area
            Mom_cuw.append(m_cuw); Mom_i.append(m_i); Mom_h.append(m_h); Mom_icm.append(m_icm)
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
        #print len(vel_z)

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


#fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axis
#plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
#plt.xlabel("common X")
#plt.ylabel("common Y")
#plt.tight_layout()
'''
    icm_c = '#03d803ff'
    if j == 0:
        plt.ylim(-0.00031,0.0062)
        #plt.ylim(-0.00445, 0.0062)
    if j == 4 : # ICM 00
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')
        icm_inflow = 0.0051002 / 16 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
        plt.ylim(-0.000275, 0.0019)
    elif j== 3: # ICM 0
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')
        icm_inflow = 0.0051002 / 4 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
        plt.ylim(-0.000275, 0.0019)
    elif j== 1: # ICM 1
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='-.')
        icm_inflow = 0.005002  # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k',lw=2.5)

        plt.ylim(-0.00031, 0.0062) # for 3 kpc
        #plt.ylim(-0.0325, 0.027) # for 2 kpc
    elif j== 2: # ICM 2
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='-.')
        icm_inflow = 0.005002 * 2 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k',lw=2.5)

        plt.ylim(-0.00655, 0.151) # for 3 kpc
        #plt.ylim(-0.033, 0.131) # for 2 kpc
    elif j== 5: # ICM 3
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')
        icm_inflow = 0.0051002 * 4 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
        plt.ylim(-0.0175, 0.39)
    elif j== 6: # ICM 4
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')
        icm_inflow = 0.0051002 * 10 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
        plt.ylim(-0.0175, 0.39)

    #plt.plot(time, Mom_t,label='Total', color='dimgrey',ls=':')
    #plt.plot(time, Mom_hi, label='Hot_%s' % labell[j], color=C3[j],ls=S[j]) # Hot only
    #plt.plot(time, Mom_c, label='%s' % labelll[0], color=C2[0])
    #plt.plot(time, Mom_u, label='%s' % labelll[1], color=C2[1])
    #plt.plot(time, Mom_w, label='%s' % labelll[2], color=C2[2])
    #plt.plot(time, Mom_h, label='%s' % labelll[3], color=C2[3])
    #plt.plot(time, Mom_i, label='%s' % labelll[4], color=C2[4])
    #plt.plot(xnew,smooth,label='%s' % labell[j],color=C[j]) ; #plt.plot(time,Mom_3,label='3kpc') ;  plt.plot(time,Mom_2,label='2kpc'), plt.plot(time,Mom_1,label='1kpc')

#plt.title('Mass Flux_UpperBound')
    #ml = MultipleLocator(5)
    #plt.minorticks_on()

    #plt.axes().tick_params(which='minor', direction='in')
    #plt.axes().tick_params(which='major', direction='in',labelsize=22)
    #if j== 1 or j ==2 or j ==3 :
    #    plt.xticks([])
    #if j == 4 or j==5 or j ==6 :
    #    plt.xlabel(r'time [Myr]',fontsize=14)
    #    plt.xticks([250,300,350,400,450])
    #if j == 1 or j == 4 :
    #    plt.ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]',fontsize=15)
    if j==1:
        plt.legend(loc=0,fontsize=19)
        plt.ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]', fontsize=24)
    if j==2:
        plt.xlabel(r'time [Myr]',fontsize=20)
        plt.xticks([250,300,350,400,450],fontsize=20)
    else :
        plt.xticks([250, 300, 350, 400, 450],[])
    #if j==3:
    #    plt.ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]', fontsize=19)
    plt.yticks(fontsize=19)
    plt.xlim(time.min(),time.max())
    plt.tick_params(direction='in')
    plt.locator_params(axis='y', nbins=6)
    #plt.title('%s' % labell[j])
    #if j == 1 :
    #    plt.legend(loc=0,fontsize=18)
    #if j != 2 :
    #    plt.ylim(-0.0002,0.0065)

    #plt.tight_layout(pad=0.3)
    #plt.savefig('D:/yeongu/plots/massflux_upper_%s.png' % labell[j], dpi=500)  # ,transparent=True
    k = k + 1
plt.subplots_adjust(bottom=0.06, top=0.99, hspace=.01, left=0.17, right=0.99)
#ax = fig.add_subplot(111,frameon=False)
#ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
#ax.spines['top'].set_color('none')
#ax.spines['bottom'].set_color('none')
#ax.spines['left'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.patch.set_facecolor('None')
#ax.set_ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]', fontsize=24,labelpad=45)
#plt.minorticks_on()
#plt.axes().tick_params(which='minor', direction='in')
#plt.axes().tick_params(which='major', direction='in',labelsize=22)
#plt.tight_layout(pad=0.3)
plt.savefig('D:/yeongu/plots/paperplot/new/massflux_No12_new2.png',dpi=300) # ,transparent=True
plt.savefig('D:/yeongu/plots/paperplot/new/massflux_No12_new2.eps',format='eps',dpi=300)
plt.show()
#plt.close()
'''


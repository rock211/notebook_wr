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

simid_t = ('R8_8pc_metal','RPS_8pc_noICM_newacc','RPS_8pc_ICM1','RPS_8pc_ICM1_newacc','RPS_8pc_ICM2','RPS_8pc_ICM2_newacc') # 'MHD_8pc_new' ,
#labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4')
labell = ('No ICM','No ICM_new','ICM1','ICM1_new','ICM2','ICM2_new')
labelll = ('Cold','Unstable','Warm','Ionized','Hot')
C = ('goldenrod','royalblue','firebrick')
#C2 = ('darkblue','deepskyblue','goldenrod','red','firebrick')
C2 = ('salmon','r','limegreen','g','deepskyblue','b')
C3 = ('darkred','red','salmon')
S = (':','--','-')


# overplot Starformation rate of three different simulations
#plt.figure(figsize=(14,8))
fig =plt.figure(figsize=(8.5,12))
for j in range(0,6): #range(1,7) :
    basedir = 'F:/yeongu/'
    if j==5:
        basedir = 'D:/yeongu/'
    simid = simid_t[j]
    Mom_up = []
    stop = 499
    #if j == 5:
    #    stop = 472
    #elif j == 6:
    #    stop = 401
    #else:
    #    stop = 500
    Mom_1 = []
    Mom_2 = []
    Mom_3 = []
    Mom_c = []
    Mom_u = []
    Mom_w = []
    Mom_h = []
    Mom_i = []
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

        # full domain information
        #print ds.domain

        # information of grid #0
        #print ds.grids[0]

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
        height = 895 #698 823 #895 #data set

        d_up = d[height]; d_down = d[0]
        vel_z_up = vel_z[height]; vel_z_down = vel_z[0]
        temp_up = temp[height]; temp_down = temp[0]

        #print np.median(temp_down), np.mean(temp[50])
        #d_c = d_up[temp_up < 184] ; v_c = vel_z_up[temp_up < 184]
        #d_u = d_up[(temp_up > 184) & (temp_up < 5050)] ; v_u = vel_z_up[(temp_up > 184) & (temp_up < 5050)]
        #d_w = d_up[(temp_up > 5050) & (temp_up < 20000)] ; v_w = vel_z_up[(temp_up > 5050) & (temp_up < 20000)]
        #d_h = d_up[(temp_up > 20000) & (temp_up < 500000)] ; v_h = vel_z_up[(temp_up > 20000) & (temp_up < 500000)]
        #d_i = d_up[temp_up > 500000] ; v_i = vel_z_up[temp_up > 500000]
        d_cuw = d_up[temp_up < 20000] ; v_cuw = vel_z_up[temp_up < 20000]
        d_hi = d_up[temp_up > 20000];   v_hi = vel_z_up[temp_up > 20000]

        d_d_cuw = d_down[temp_up < 20000] ; v_d_cuw = vel_z_down[temp_up < 20000]
        d_d_hi = d_down[temp_up > 20000]; v_d_hi = vel_z_down[temp_up > 20000]
        ###### momentum #######
        #moment_z = d*vel_z/unit['time'].value #SolMass/pc2/Myr
        #m_c = np.sum(np.sum(d_c * v_c / unit['time'].value))*512/(8*1024*1024)
        #m_u = np.sum(np.sum(d_u * v_u / unit['time'].value))*512/(8*1024*1024)
        #m_w = np.sum(np.sum(d_w * v_w / unit['time'].value))*512/(8*1024*1024)
        #m_h = np.sum(np.sum(d_h * v_h / unit['time'].value))*512/(8*1024*1024)
        #m_i = np.sum(np.sum(d_i * v_i / unit['time'].value))*512/(8*1024*1024)
        m_cuw = np.sum(np.sum(d_cuw * v_cuw / unit['time'].value)) * 512 / (8 * 1024 * 1024)
        m_hi = np.sum(np.sum(d_hi * v_hi / unit['time'].value)) * 512 / (8 * 1024 * 1024)


        #tot = m_cuw+m_hi
        #Mom_cuw.append(m_cuw); Mom_hi.append(m_hi)
        if j==0 or j==1:
            m_d_cuw = np.sum(np.sum(d_d_cuw * v_d_cuw / unit['time'].value)) * 512 / (8 * 1024 * 1024)
            m_d_hi = np.sum(np.sum(d_d_hi * v_d_hi / unit['time'].value)) * 512 / (8 * 1024 * 1024)
            #print m_d_cuw+m_d_hi
            Mom_cuw.append(m_cuw - m_d_cuw); Mom_hi.append(m_hi - m_d_hi)
        elif j==2 or j==4:
            scalar = ds.read_all_data('specific_scalar4')  # ism = 0 / icm = 1
            scalar_up = scalar[height]
            m_icm = np.sum(np.sum(scalar_up * d_up * vel_z_up / unit['time'].value)) * 512 / (8 * 1024 * 1024)
            Mom_cuw.append(m_cuw); Mom_hi.append(m_hi); Mom_icm.append(m_icm)
        elif j == 3 or j == 5:
            scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1
            scalar_up = scalar[height]
            m_icm = np.sum(np.sum(scalar_up * d_up * vel_z_up / unit['time'].value)) * 512 / (8 * 1024 * 1024)
            Mom_cuw.append(m_cuw);
            Mom_hi.append(m_hi);
            Mom_icm.append(m_icm)

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

    time = np.arange(250*unit['time'].value,stop*unit['time'].value,unit['time'].value)
    #xnew = np.linspace(time.min(),time.max(),5000)
    #smooth = spline(time,Mom_up,xnew)

    if j==0 or j==1:
        k=1
    elif j==2 or j ==3 :
        k=2
    elif j==4 or j ==5:
        k=3

    plt.subplot(3,1,k)
    plt.plot(time, Mom_cuw, label=labell[j]+'cold', color=C2[j],ls='--')
    plt.plot(time, Mom_hi, label=labell[j]+'hot', color=C2[j])
    icm_c = 'green'
    '''
    if j == 0:
        #plt.ylim(-0.00031,0.0062)
        plt.ylim(-0.00445, 0.0062)
    if j == 1 : # ICM 00
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')
        icm_inflow = 0.0051002 / 16 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
        plt.ylim(-0.000275, 0.0019)
    elif j== 2: # ICM 0
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')
        icm_inflow = 0.0051002 / 4 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
        plt.ylim(-0.000275, 0.0019)
    elif j== 3: # ICM 1
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')
        icm_inflow = 0.0051002  # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
        #plt.ylim(-0.000275, 0.0055)
        #plt.ylim(-0.00031, 0.0062) # for 3 kpc
        plt.ylim(-0.0325, 0.027) # for 2 kpc
    elif j== 4: # ICM 2
        plt.plot(time, Mom_icm, label='ICM outflow', color=icm_c,ls='--')
        icm_inflow = 0.0051002 * 2 # m_d_cuw + m_d_hi
        plt.axhline(icm_inflow,label='ICM inflow',ls=':',color='k')
        #plt.ylim(-0.0175, 0.39)
        #plt.ylim(-0.00655, 0.131) # for 3 kpc
        plt.ylim(-0.033, 0.131) # for 2 kpc
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
    '''
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
    if j==3:
        plt.legend(loc=0,fontsize=19)
    elif j==1:
        plt.legend(loc=0, fontsize=19)
    elif j==5:
        plt.legend(loc=0, fontsize=19)

    if j==5:
        plt.xlabel(r'time [Myr]',fontsize=20)
        plt.xticks([250,300,350,400,450],fontsize=20)
    else :
        plt.xticks([250, 300, 350, 400, 450],[])
    #if j==3:
    #    plt.ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]', fontsize=19)
    plt.yticks(fontsize=19)
    plt.xlim(time.min(),time.max())
    plt.tick_params(direction='in')
    #plt.title('%s' % labell[j])
    #if j == 1 :
    #    plt.legend(loc=0,fontsize=18)
    #if j != 2 :
    #    plt.ylim(-0.0002,0.0065)

    #plt.tight_layout(pad=0.3)
    #plt.savefig('D:/yeongu/plots/massflux_upper_%s.png' % labell[
    # j], dpi=500)  # ,transparent=True
ax = fig.add_subplot(111,frameon=False)
ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.patch.set_facecolor('None')
ax.set_ylabel(r'Mass flux [M$_{\odot}$ pc$^{-2}$ Myr$^{-1}$]', fontsize=24,labelpad=45)
#plt.minorticks_on()
#plt.axes().tick_params(which='minor', direction='in')
#plt.axes().tick_params(which='major', direction='in',labelsize=22)
plt.tight_layout(pad=0.3)
plt.savefig('D:/yeongu/plots/massflux_compare.png',dpi=400) # ,transparent=True
#plt.savefig('D:/yeongu/plots/massflux_No12_2kpc.eps',format='eps',dpi=400)
plt.show()
#plt.close()



import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
from scipy.interpolate import spline
from matplotlib.ticker import MultipleLocator
import cPickle



sys.path.insert(0, '../')
from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

import pyathena as pa


unit = pa.set_units(muH=1.4271)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
#print (unit['mass'], unit['time'], unit['magnetic_field'])

Msun = unit['mass'].to('Msun').value
print Msun
Myr=unit['time'].to('Myr').value

agebin = 10 # unit : Myr, 10 : H-alpha like, 40 : cluster lifetime

M_T_0 = []
M_T_1 = []
M_T_2 = []

simid_t = ('R8_8pc_metal','RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1','RPS_8pc_ICM2','RPS_8pc_ICM3','RPS_8pc_ICM4') # 'MHD_8pc_new' ,
labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # r'No ICM',
#C = ('k',(1,0.796,0.224),(1,0.506,0.035),(0.871,0.0,0.102),'deepskyblue','blue','navy')#('k','plum','orchid','purple','deepskyblue','royalblue','navy')#'lightgreen','forestgreen','darkgreen') #'darkkhaki','royalblue','firebrick'
C = ('k','lightsalmon','skyblue','darkmagenta','goldenrod','royalblue','crimson')
S = ('--','-','-','-','-','-','-')
#L = (1.5,1.5,1.5)
lw = (4,3,2.5,2,3,2.5,2)
alpha = (0.5,1,1,1,0.8,1,1)
# overplot Starformation rate of three different simulations
start = 250
stop = 500
plt.figure(figsize=(11,10))
for j in (0,1,2,3,4,5,6) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j]
    Mj = []
    M_star_c = 0.0
    Mass_tot = []
    Mc = []
    SFE = []
    if j == 5:
        stop = 472
    elif j == 6:
        stop = 401
    else:
        stop = 500
    for tidx in range(250, stop):
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir,simid, simid, tidx)

        # read in domain information
        #ds = pa.AthenaDataSet(vtkfname)

        # full domain information
        #ds.domain

        # information of grid #0
        #ds.grids[0]

        # yet, we didn't read data.
        # let's read each data field in a full domain

        # this can be original data fields
        #d = ds.read_all_data('density')*unit['density'].value
        #nd = ds.read_all_data('number_density')
        #tem = ds.read_all_data('temperature')
        #coolftn = pa.coolftn()
        #temp = coolftn.get_temp(P / d)  # Temperature derived from P/d & cooling func
        #d = ds.read_all_data('density')*unit['density'].value
        #T1 = ds.read_all_data('T1')
        #coolftn = pa.coolftn()
        #temp = coolftn.get_temp(T1)

        starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)
        #print(sp.shape)
        # print (sp)

        print(tidx)
        star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
        star_clu2 = star_clu[star_clu['age'] * Myr < agebin]
        #print(star_clu['x1'])
        #M_star = []
        #for age in range(0,250):
        #    star_clu2 = star_clu[star_clu['age'] * Myr < age] # time scale cut
            #print star_clu2['mass']
        #    M_star.append(np.sum(star_clu2['mass']))
            #print M_star

        #total_mass = d * 8*8*8
        #cold_warm_tot = total_mass[temp < 20000]
        #surf_mass = np.sum(np.sum(np.sum(cold_warm_tot/(1024*1024))))
        #Mass_tot.append(np.sum(np.sum(np.sum(surf_mass))))

        M_star_c += sum(star_clu2['mass'])*Msun # mass sum in time scale / cumulative
        M_star = sum(star_clu2['mass']) * Msun # SFR time variation

        Mj.append(M_star* unit['time'].value / (1e+6*agebin*1.024*1.024)) #* unit['time'].value/ mass divided by time scale & area
        Mc.append(M_star_c* unit['time'].value / (1e+6*agebin*1.024*1.024)) #* unit['time'].value/ mass divided by time scale & area

        #sfe = M_star* unit['time'].value / (1e+6*agebin*1.024*1.024)/surf_mass
        #SFE.append(sfe)

    #print M_star
    #np.array([M_star],dtype=object)
    #agee = range(250, 250+len(M_star))
    #plt.plot(agee,M_star)
    #plt.show()
        print j
        ut = round(unit['time'].value,2)
        xx = np.arange(start, stop)*ut
        ml = MultipleLocator(5)


        if j==0 or j==1 or j==2 or j==3:
            ax1=plt.subplot(2,2,1)
            #plt.minorticks_on()
            ax1.tick_params(which='minor', direction='in')
            ax1.tick_params(which='major', direction='in', labelsize=16)
            plt.minorticks_on()
            plt.plot(xx, Mj, label='%s' % labell[j], color=C[j],linewidth=lw[j],alpha=alpha[j],ls=S[j])
            plt.xticks([250, 300, 350, 400, 450], [])
            #plt.plot(xx, SFE, label='%s' % labell[j], color=C[j], linestyle=S[j], linewidth=L[j], alpha=alpha[j])
            plt.xlim(xx.min(), 499*ut)
            plt.ylim(0, 0.019)
            #plt.xlabel(r'time [Myr]', fontsize=14)
            # plt.ylabel(r'SFE',fontsize=17) # For SFE
            plt.ylabel(r'$\Sigma_{SFR}$ [M$_{\odot}$ kpc$^{-2}$ yr$^{-1}$]',fontsize=19) # For SFR
            # plt.title(r'SFR variation $(Agebin=%s)$' % agebin)
            # plt.yticks([0,0.005,0.01,0.015,0.02]) # For cumulative
            plt.yticks([0,0.004,0.008,0.012,0.016]) # For SFR

            ax2=plt.subplot(2,2,2)
            plt.minorticks_on()
            ax2.tick_params(which='minor', direction='in')
            ax2.tick_params(which='major', direction='in', labelsize=16)
            plt.plot(xx,Mc,label='%s' % labell[j], color=C[j],linewidth=lw[j],alpha=alpha[j],ls=S[j])
            #plt.xlabel(r'time [Myr]', fontsize=14)
            plt.xticks([250,300,350,400,450],[])
            plt.ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ pc$^{-2}$]', fontsize=19)  # For cumulative
            plt.xlim(xx.min(), 499*ut)
            plt.ylim(0, 1.1)  # For cumulative
            plt.legend(loc=0, fontsize=17)
            plt.tight_layout(pad=0.3)
            plt.savefig('D:/yeongu/plots/paperplot/sfr/SFR_%s_%s.png' % (j, tidx), dpi=100)

        if j==0 or j==4 or j==5 or j==6:
            ax1 = plt.subplot(2, 2, 3)
            # plt.minorticks_on()
            ax1.tick_params(which='minor', direction='in')
            ax1.tick_params(which='major', direction='in', labelsize=16)
            plt.minorticks_on()
            plt.plot(xx, Mj, label='%s' % labell[j], color=C[j], linewidth=lw[j], alpha=alpha[j], ls=S[j])
            # plt.plot(xx, SFE, label='%s' % labell[j], color=C[j], linestyle=S[j], linewidth=L[j], alpha=alpha[j])
            plt.xlim(xx.min(), 499 * ut)
            plt.ylim(0, 0.019)
            plt.xlabel(r'time [Myr]', fontsize=16)
            # plt.ylabel(r'SFE',fontsize=17) # For SFE
            plt.ylabel(r'$\Sigma_{SFR}$ [M$_{\odot}$ kpc$^{-2}$ yr$^{-1}$]', fontsize=19)  # For SFR
            # plt.title(r'SFR variation $(Agebin=%s)$' % agebin)
            # plt.yticks([0,0.005,0.01,0.015,0.02]) # For cumulative
            plt.yticks([0, 0.004, 0.008, 0.012, 0.016])  # For SFR

            ax2 = plt.subplot(2, 2, 4)
            plt.minorticks_on()
            ax2.tick_params(which='minor', direction='in')
            ax2.tick_params(which='major', direction='in', labelsize=16)
            plt.plot(xx, Mc, label='%s' % labell[j], color=C[j], linewidth=lw[j], alpha=alpha[j], ls=S[j])
            plt.xlabel(r'time [Myr]', fontsize=16)
            plt.ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ pc$^{-2}$]', fontsize=19)  # For cumulative
            plt.xlim(xx.min(), 499 * ut)
            plt.ylim(0, 1.1)  # For cumulative
            plt.legend(loc=0, fontsize=17)
            plt.tight_layout(pad=0.3)
            plt.show()
            #plt.savefig('D:/yeongu/plots/paperplot/sfr/SFR_%s_%s.png' % (j, tidx), dpi=100)

    #plt.savefig('D:/yeongu/plots/paperplot/sfr/SFR_%s_%s.png' % (j,tidx), dpi=100)
    #plt.savefig('D:/yeongu/plots/paperplot/SFR_%s.eps' % j,format='eps', dpi=300)
    #plt.show()
    #plt.close()

'''
    ax1=plt.subplot(1,2,1)
    #plt.minorticks_on()
    ax1.tick_params(which='minor', direction='in')
    ax1.tick_params(which='major', direction='in', labelsize=12)
    plt.minorticks_on()
    plt.plot(xx, Mj, label='%s' % labell[j], color=C[j],linewidth=lw[j],alpha=alpha[j],ls=S[j])
    #plt.plot(xx, SFE, label='%s' % labell[j], color=C[j], linestyle=S[j], linewidth=L[j], alpha=alpha[j])
    plt.xlim(xx.min(), 499*ut)
    plt.ylim(0, 0.019)
    plt.xlabel(r'time [Myr]', fontsize=14)
    # plt.ylabel(r'SFE',fontsize=17) # For SFE
    plt.ylabel(r'$\Sigma_{SFR}$ [M$_{\odot}$ kpc$^{-2}$ yr$^{-1}$]',fontsize=17) # For SFR
    # plt.title(r'SFR variation $(Agebin=%s)$' % agebin)
    # plt.yticks([0,0.005,0.01,0.015,0.02]) # For cumulative
    plt.yticks([0,0.004,0.008,0.012,0.016]) # For SFR
    
    ax2=plt.subplot(1,2,2)
    plt.minorticks_on()
    ax2.tick_params(which='minor', direction='in')
    ax2.tick_params(which='major', direction='in', labelsize=12)
    plt.plot(xx,Mc,label='%s' % labell[j], color=C[j],linewidth=lw[j],alpha=alpha[j],ls=S[j])
    plt.xlabel(r'time [Myr]', fontsize=14)
    plt.ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ pc$^{-2}$]', fontsize=17)  # For cumulative
    plt.xlim(xx.min(), 499*ut)
    plt.ylim(0, 1.1)  # For cumulative
    plt.legend(loc=0, fontsize=15)
    plt.tight_layout(pad=0.3)
#plt.savefig('D:/yeongu/plots/paperplot/SFR_%s.png' % j, dpi=300)
plt.show()
plt.close()
'''
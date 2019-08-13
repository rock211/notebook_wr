import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
from scipy.interpolate import spline
from matplotlib.ticker import MultipleLocator
import cPickle
import seaborn


sys.path.insert(0, '../')
from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

import pyathena as pa

unit = pa.set_units(muH=1.4271)
# print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'])

Msun = unit['mass'].to('Msun').value
print Msun
Myr = unit['time'].to('Myr').value

agebin = 10  # unit : Myr, 10 : H-alpha like, 40 : cluster lifetime

M_T_0 = []
M_T_1 = []
M_T_2 = []

simid_t = ('RPS_8pc_noICM_newacc', 'RPS_8pc_ICM0_newacc', 'RPS_8pc_ICM1_newacc','RPS_4pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc','RPS_4pc_ICM2_newacc', 'RPS_8pc_ICM3_newacc')  # 'MHD_8pc_new' ,
labell = ('No ICM','Very Weak', 'Weak', 'Strong','Very Strong' ,'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
labell = ('No ICM','P1', 'P3','P3h', 'P7','P7h','P14' ,'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
# C = ('k',(1,0.796,0.224),(1,0.506,0.035),(0.871,0.0,0.102),'deepskyblue','blue','navy')#('k','plum','orchid','purple','deepskyblue','royalblue','navy')#'lightgreen','forestgreen','darkgreen') #'darkkhaki','royalblue','firebrick'
C = ('k', 'lightskyblue', 'dodgerblue','mediumblue' ,'goldenrod','salmon', 'firebrick')
#C = ('dimgray', 'lightskyblue','mediumblue' ,'','goldenrod','', 'firebrick','darkmagenta','goldenrod','royalblue','crimson') # 'plum','orchid','purple'
S = ('--', '-', '-', '-', '-', '-', '-')
# L = (1.5,1.5,1.5)
lw = (2.7, 1.3, 2.4, 3.5, 2.4, 3.5, 1.3)
alpha = (1, 1, 1, 1, 1, 1, 1)
# overplot Starformation rate of three different simulations
start = 250
stop = 499
plt.figure(figsize=(12, 8))
SFR =[]
ut = round(unit['time'].value, 2)
xx = np.arange(start, stop) * ut
#SFR.append(xx)
rev = (6,5,4,3,2,1,0)
nor = (0,1,2,3,4,5,6)
k=0
for j in (6,1,2,4,3,5,0):

    basedir = 'G:/yeongu/'

    simid = simid_t[j]
    Mj = []
    M_star_c = 0.0
    Mass_tot = []
    Mc = []
    SFE = []

    if j == 6 or j==5:
        stop = 474
    else:
        stop = 499
    #stop = 499
    for tidx in range(250, stop):
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)

        # read in domain information
        # ds = pa.AthenaDataSet(vtkfname)

        # full domain information
        # ds.domain

        # information of grid #0
        # ds.grids[0]

        # yet, we didn't read data.
        # let's read each data field in a full domain

        # this can be original data fields
        # d = ds.read_all_data('density')*unit['density'].value

        if j==6 or j==5:
            if tidx < 474:
                starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
                sp = pa.read_starvtk(starfname)
                # print(sp.shape)
                # print (sp)


                star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
                star_clu2 = star_clu[star_clu['age'] * Myr < agebin]
                M_star_c += sum(star_clu2['mass']) * Msun  # mass sum in time scale / cumulative
                M_star = sum(star_clu2['mass']) * Msun  # SFR time variation
            else:
                M_star = 0

        else :
            starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
            sp = pa.read_starvtk(starfname)
            # print(sp.shape)
            # print (sp)

            star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
            star_clu2 = star_clu[star_clu['age'] * Myr < agebin]
            M_star_c += sum(star_clu2['mass']) * Msun  # mass sum in time scale / cumulative
            M_star = sum(star_clu2['mass']) * Msun  # SFR time variation

        Mj.append(M_star * unit['time'].value / (
                    1e+6 * agebin * 1.024 * 1.024))  # * unit['time'].value/ mass divided by time scale & area
        Mc.append(M_star_c * unit['time'].value / (
                    1e+6 * agebin * 1.024 * 1.024))  # * unit['time'].value/ mass divided by time scale & area

        # sfe = M_star* unit['time'].value / (1e+6*agebin*1.024*1.024)/surf_mass
        # SFE.append(sfe)

        # print M_star
        # np.array([M_star],dtype=object)
        # agee = range(250, 250+len(M_star))
        # plt.plot(agee,M_star)
        # plt.show()
        #print j

        ml = MultipleLocator(5)
        print(tidx)
    print len(Mj)

    ax1 = plt.subplot(2, 1, 1)
    # plt.minorticks_on()
    ax1.tick_params(which='minor', direction='in',length=5,width=1.5)
    ax1.tick_params(which='major', direction='in', labelsize=21,length=9,width=1.5)
    plt.minorticks_on()
    plt.plot(xx, Mj, color=C[j], linewidth=lw[j], alpha=alpha[j], ls=S[j])

    plt.xticks([250, 300, 350, 400, 450], [])
    # plt.plot(xx, SFE, label='%s' % labell[j], color=C[j], linestyle=S[j], linewidth=L[j], alpha=alpha[j])
    plt.xlim(xx.min(), 499 * ut)
    plt.ylim(0, 0.019)
    # plt.xlabel(r'time [Myr]', fontsize=14)
    # plt.ylabel(r'SFE',fontsize=17) # For SFE
    plt.ylabel(r'$\Sigma_{SFR}$ [M$_{\odot}$ kpc$^{-2}$ yr$^{-1}$]', fontsize=23)  # For SFR
    # plt.title(r'SFR variation $(Agebin=%s)$' % agebin)
    # plt.yticks([0,0.005,0.01,0.015,0.02]) # For cumulative
    plt.yticks([0, 0.004, 0.008, 0.012, 0.016])  # For SFR

    ax2 = plt.subplot(2, 1, 2)
    plt.minorticks_on()
    ax2.tick_params(which='minor', direction='in',length=5,width=1.5)
    ax2.tick_params(which='major', direction='in', labelsize=22,length=9,width=1.5)
    plt.plot(xx, Mc, color=C[j], linewidth=lw[j], alpha=alpha[j], ls=S[j])
    plt.plot([],[],label='%s' % labell[k], color=C[k], linewidth=lw[k], ls=S[k])
    k = k+1
    plt.xlabel(r'Time [Myr]', fontsize=22)
    #plt.xticks([250, 300, 350, 400, 450], [])
    plt.ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ pc$^{-2}$]', fontsize=23)  # For cumulative
    plt.xlim(xx.min(), 499 * ut)
    plt.ylim(0, 1.2)  # For cumulative
    plt.tight_layout(pad=0.3)
    #SFR.append(Mj)
#SFR=np.array(SFR)
#np.savetxt('SFR.txt',SFR,fmt='%s')

#handles, labels = ax2.get_legend_handles_labels()
#print handles[0], labels[0]
#ax2.legend(handles[::-1],labels[::-1],loc=0, fontsize=15)
ax2.legend(loc=0,fontsize=16)
plt.savefig('D:/yeongu/plots/paperplot/new/SFR_test_w_4.png' , dpi=300)
plt.savefig('D:/yeongu/plots/paperplot/new/SFR_test_w_4.eps' ,format='eps', dpi=300,rasterized=True)
plt.show()



    # plt.savefig('D:/yeongu/plots/paperplot/sfr/SFR_%s_%s.png' % (j,tidx), dpi=100)
    # plt.savefig('D:/yeongu/plots/paperplot/SFR_%s.eps' % j,format='eps', dpi=300)
    # plt.show()
    # plt.close()

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
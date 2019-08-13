import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator
import copy
from mpl_toolkits import axes_grid1

from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

labell = ('No ICM', 'P1', 'P3', 'P3h', 'P7', 'P7h', 'P14')  # r'No ICM',
# C = ('gray', 'mediumturquoise', 'dodgerblue','mediumblue' ,'goldenrod','salmon', 'firebrick','darkmagenta','goldenrod','royalblue','crimson') # 'plum','orchid','purple'
C = ('k', 'lightskyblue', 'dodgerblue', 'mediumblue', 'salmon', 'crimson', 'goldenrod')
Model = [0, 8.63 * 1e3, 3.46 * 1e4, 3.46 * 1e4, 6.92 * 1e4, 6.92 * 1e4, 1.38 * 1e5]
crit = 94
k = 1

plt.figure(figsize=(10,6))

#plt.figure(figsize=(5, 5))
for j in (2, 4):
    # plt.figure(figsize=(6,4))
    variable = 'density'
    if j==2 or j==4:
        z = range(895)
    else:
        z = range(1791)
    P = np.genfromtxt('./proj/Pism_%s.txt' % labell[j])
    SFR = np.genfromtxt('./proj/SFR_%s.txt' % labell[j])
    W = np.genfromtxt('./proj/Weight_%s.txt' % labell[j])
    D = np.genfromtxt('./proj/Dc_%s.txt' % labell[j])
    if j != 0:
        icm_tot = np.genfromtxt('./proj/Picm_tot_%s.txt' % labell[j])
        icm_th = np.genfromtxt('./proj/Picm_ther_%s.txt' % labell[j])
        icm_tu = np.genfromtxt('./proj/Picm_turb_%s.txt' % labell[j])
    # print D.shape
    # for t in range(50):
    #    plt.semilogy(P[t,:],c='k')
    # plt.axhline(np.max(Tot),ls='--')
    #    print np.max(P[t,:])
    # plt.ylim(1e2,1e5)
    # plt.show()
    # plt.close()
    # print P[crit::,:].shape
    P_e = P[0:crit, :];
    P_l = P[crit::, :]  # +icm_tot[0:crit,0:-1] +icm_tot[crit::,0:-1]
    SFR_e = SFR[0:crit];
    SFR_l = SFR[crit::]
    W_e = W[0:crit, :];
    W_l = W[crit::, :]
    D_e = D[0:crit, :];
    D_l = D[crit::, :]

    P_medi_e = np.nanmedian(P_e, axis=0);
    P_medi_l = np.nanmedian(P_l, axis=0)
    P_min_e = np.nanpercentile(P_e, 25, axis=0);
    P_max_e = np.nanpercentile(P_e, 75, axis=0)
    P_min_l = np.nanpercentile(P_l, 25, axis=0);
    P_max_l = np.nanpercentile(P_l, 75, axis=0)

    SFR_medi_e = np.nanmedian(SFR_e);
    SFR_medi_l = np.nanmedian(SFR_l)
    SFR_min_e = np.nanpercentile(SFR_e, 25);
    SFR_min_l = np.nanpercentile(SFR_l, 25)
    SFR_max_e = np.nanpercentile(SFR_e, 75);
    SFR_max_l = np.nanpercentile(SFR_l, 75)

    W_medi_e = np.nanmean(W_e, axis=0);
    W_medi_l = np.nanmean(W_l, axis=0)
    W_min_e = np.nanpercentile(W_e, 10, axis=0);
    W_min_l = np.nanpercentile(W_l, 10, axis=0)
    W_max_e = np.nanpercentile(W_e, 90, axis=0);
    W_max_l = np.nanpercentile(W_l, 90, axis=0)
    print np.nanmax(W_e)
    print np.nanmax(W_l)
    D_medi_e = np.nanmedian(D_e, axis=0) * 1e6;
    D_medi_l = np.nanmedian(D_l, axis=0) * 1e6
    D_min_e = np.nanpercentile(D_e, 25, axis=0);
    D_min_l = np.nanpercentile(D_l, 25, axis=0)
    D_max_e = np.nanpercentile(D_e, 75, axis=0);
    D_max_l = np.nanpercentile(D_l, 75, axis=0)

    # ICM calculation below

    icmtot_e = icm_tot[0:crit, :]; icmtot_l = icm_tot[crit::, :]
    icmth_e = icm_th[0:crit, :]; icmth_l = icm_th[crit::, :]
    icmtu_e = icm_tu[0:crit, :]; icmtu_l = icm_tu[crit::, :]

    icmtot_medi_e = np.nanmean(icmtot_e,axis=0); icmtot_medi_l = np.nanmean(icmtot_l,axis=0)
    icmtot_min_e = np.nanpercentile(icmtot_e,10,axis=0) ; icmtot_min_l = np.nanpercentile(icmtot_l,10,axis=0)
    icmtot_max_e = np.nanpercentile(icmtot_e, 90, axis=0) ; icmtot_max_l = np.nanpercentile(icmtot_l,90,axis=0)

    icmth_medi_e = np.nanmean(icmth_e,axis=0); icmth_medi_l = np.nanmean(icmth_l,axis=0)
    icmth_min_e = np.nanpercentile(icmth_e,10,axis=0) ; icmth_min_l = np.nanpercentile(icmth_l,10,axis=0)
    icmth_max_e = np.nanpercentile(icmth_e, 90, axis=0) ; icmth_max_l = np.nanpercentile(icmth_l,90,axis=0)

    icmtu_medi_e = np.nanmean(icmtu_e,axis=0); icmtu_medi_l = np.nanmean(icmtu_l,axis=0)
    icmtu_min_e = np.nanpercentile(icmtu_e,10,axis=0) ; icmtu_min_l = np.nanpercentile(icmtu_l,10,axis=0)
    icmtu_max_e = np.nanpercentile(icmtu_e, 90, axis=0) ; icmtu_max_l = np.nanpercentile(icmtu_l,90,axis=0)


    d_e_p = np.argmax(D_medi_e);
    d_l_p = np.argmax(D_medi_l)  # Density profile peak position

    print labell[j], d_e_p, d_l_p
    # print P_medi_e.shape
    # print SFR_medi_e
    # print W_medi_e.shape
    # print D_medi_e.shape

    # print D_medi_e[0:-1]

    plt.subplot(2,2,2*k-1)
    plt.plot(z,icmtot_medi_e[0:-1],'m',label='P$_{ICM}$')
    plt.fill_between(z,icmtot_min_e[0:-1],icmtot_max_e[0:-1],facecolor='m',alpha=0.3)
    plt.plot(z,icmth_medi_e[0:-1],c='coral',ls='-',label='P$_{ICM,ther}$')
    plt.fill_between(z, icmth_min_e[0:-1], icmth_max_e[0:-1], facecolor='coral', alpha=0.3)
    plt.plot(z,icmtu_medi_e[0:-1],c='gold',ls='-',label='P$_{ICM,turb}$')
    plt.fill_between(z, icmtu_min_e[0:-1], icmtu_max_e[0:-1], facecolor='gold', alpha=0.3)
    plt.plot(z,W_medi_e,'c-',label='W$_{ISM}$')
    plt.fill_between(z,W_min_e,W_max_e,facecolor='c',alpha=0.3)
    plt.axhline(Model[j], ls='-.', c='k', label='P$_{ICM, Input}$')

    plt.yscale('log', nonposy='clip')
    plt.ylim(5e2, 1e5)
    if j==3 or j==5:
        plt.xlim(0, 1792)
    else:
        plt.xlim(0, 896)
    plt.ylabel('Pressure/$k_B$, Weight [K cm$^{-3}$]',fontsize=12)
    if k==2:
        plt.xlabel('z [kpc]',fontsize=13)
        if j==3 or j==5:
            plt.xticks([73*2, 198*2, 323*2, 448*2, 573*2, 698*2, 823*2], ['-3', '-2', '-1', '0', '1', '2', '3'])
        else:
            plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    else:
        if j==3 or j==5:
            plt.xticks([73*2, 198*2, 323*2, 448*2, 573*2, 698*2, 823*2], [])
        else:
            plt.xticks([73, 198, 323, 448, 573, 698, 823], [])

    plt.tick_params(which='both', direction='in')
    if j==3 and k==1:
        plt.legend(loc='upper right')

    plt.subplot(2,2,2*k)
    plt.plot(z,icmtot_medi_l[0:-1],'m',label='P$_{ICM}$_%s' % labell[j])
    plt.fill_between(z,icmtot_min_l[0:-1],icmtot_max_l[0:-1],facecolor='m',alpha=0.3)
    plt.plot(z,icmth_medi_l[0:-1],c='coral',label='P$_{ICM,ther}$_%s' % labell[j])
    plt.fill_between(z, icmth_min_l[0:-1], icmth_max_l[0:-1], facecolor='coral', alpha=0.3)
    plt.plot(z,icmtu_medi_l[0:-1],c='gold',label='P$_{ICM,turb}$_%s' % labell[j])
    plt.fill_between(z, icmtu_min_l[0:-1], icmtu_max_l[0:-1], facecolor='gold', alpha=0.3)
    plt.plot(z,W_medi_l,'c-',label='W$_{ISM}$_%s' % labell[j])
    plt.fill_between(z,W_min_l,W_max_l,facecolor='c',alpha=0.3)

    #plt.plot([],[],label='Early',ls='-',c='k')
    plt.axhline(Model[j], ls='-.', c='k', label='P$_{ICM, Input}$')
    #plt.plot([],[],label='Late',ls='-',c='k')


    #h, l = plt.gca().get_legend_handles_labels()
    #legend1 = plt.legend(h[:4], l[:4], loc='upper right')  # ,fontsize=14,framealpha=0.3)
    #if k==1:
    #    plt.legend(h[4:], l[4:], loc='upper left',ncol=2 ,fontsize=11.5)
    #    plt.gca().add_artist(legend1)


    #plt.fill_between(z,P_min_e,P_max_e,facecolor=C[j],alpha=0.5)
    #plt.fill_between(z,P_min_l,P_max_l,facecolor=C[j],alpha=0.5)

    #plt.plot(z, P_medi_e, c=C[j], ls='-', label='P_Early_%s' % labell[j])
    #plt.plot(z, P_medi_l, c=C[j], ls='--', label='P_Late_%s' % labell[j])

    #plt.plot(z,W_medi_e,'m-',label='W_Early_%s' % labell[j])
    #plt.fill_between(z,W_min_e,W_max_e,facecolor=C[j],alpha=0.5)
    #plt.plot(z,W_medi_l,'m--',label='W_Late_%s' % labell[j])
    #print np.max(D_medi_e)
    #print D_medi_e[0:-1].shape
    #print len(z)
    #plt.plot(z,D_medi_e[:-1],'y-',label='D_Early_%s' % labell[j])
    #plt.fill_between(z,W_min_e,W_max_e,facecolor=C[j],alpha=0.5)
    #plt.plot(z,D_medi_l[:-1],'y--',label='D_Late_%s' % labell[j])

    #e= np.argmax(W_medi_e)
    #l = np.argmax(W_medi_l)
    #print e, l
    #plt.axvline(e,ls='-',c='r')
    #plt.axvline(l,ls='--',c='b')
    #plt.fill_between(z,W_min_l,W_max_l,facecolor=C[j],alpha=0.5)
    plt.yscale('log', nonposy='clip')
    plt.ylim(5e2, 1e5)
    if j == 3 or j == 5:
        plt.xlim(0, 1792)
    else:
        plt.xlim(0, 896)
    #plt.ylabel('Pressure [P/$k_b$]')
    if k==2:
        plt.xlabel('z [kpc]',fontsize=13)
        if j==3 or j==5:
            plt.xticks([73*2, 198*2, 323*2, 448*2, 573*2, 698*2, 823*2], ['-3', '-2', '-1', '0', '1', '2', '3'])
        else:
            plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    else:
        if j==3 or j==5:
            plt.xticks([73*2, 198*2, 323*2, 448*2, 573*2, 698*2, 823*2], [])
        else:
            plt.xticks([73, 198, 323, 448, 573, 698, 823], [])
    #else:
    #    plt.xticks([73, 198, 323, 448, 573, 698, 823], [])
    plt.tick_params(which='both', direction='in',labelleft=False)
    #plt.legend(loc='upper right')

    k = k+1
    plt.tight_layout()
    plt.savefig('D:/yeongu/plots/paperplot/new/prof_new_8pc.png',dpi=300)
    plt.savefig('D:/yeongu/plots/paperplot/new/prof_new_8pc.eps',format='eps',dpi=300)
plt.show()
#plt.close()

'''
    # print np.argmax(W_medi_e)
    # print np.argmax(W_medi_l)
    # if j!=1:
    pe = P_medi_e[d_e_p];
    pl = P_medi_l[d_l_p]
    we = W_medi_e[d_e_p];
    wl = W_medi_l[d_l_p]
    # else:
    #    w_e_p = np.argmax(W_medi_e) ; w_l_p = np.argmax(W_medi_l)
    #    pe = P_medi_e[w_e_p] ; pl = P_medi_l[w_l_p]
    #    we = W_medi_e[w_e_p] ; wl = W_medi_l[w_l_p]

    SFR_e_errors = [SFR_medi_e - SFR_min_e, SFR_max_e - SFR_medi_e]  # [np.nonzero(SFR_e)] nonzero selection
    SFR_l_errors = [SFR_medi_l - SFR_min_l, SFR_max_l - SFR_medi_l]
    P_e_errors = [pe - P_min_e[d_e_p], P_max_e[d_e_p] - pe]
    P_l_errors = [pl - P_min_l[d_l_p], P_max_l[d_l_p] - pl]
    W_e_errors = [we - W_min_e[d_e_p], W_max_e[d_e_p] - we]
    W_l_errors = [wl - W_min_l[d_l_p], W_max_l[d_l_p] - wl]
    ############################# Plot ################################

    plt.errorbar(wl, SFR_medi_l, xerr=np.array([W_l_errors]).T, yerr=np.array([SFR_l_errors]).T, color='none',
                 ecolor=C[j], marker='o', mfc=C[j], mec=C[j], capsize=3.5, label=labell[j])
    plt.errorbar(we, SFR_medi_e, xerr=np.array([W_e_errors]).T, yerr=np.array([SFR_e_errors]).T, color='none',
                 ecolor=C[j], marker='x', mfc='none', mec=C[j], capsize=3.5, ms=10)

    # plt.errorbar(we, pe,  xerr=np.array([W_e_errors]).T, yerr=np.array([P_e_errors]).T,color='none',marker='o', mfc=C[j], mec=C[j], capsize=3.5, label=labell[j],ecolor=C[j])
    # plt.errorbar(wl, pl,  xerr=np.array([W_l_errors]).T, yerr=np.array([P_l_errors]).T,marker='s', mfc='none', mec=C[j], capsize=3.5,ecolor=C[j])

##### legend #####
s1 = plt.scatter([], [], marker='o', edgecolors='k', facecolors='k', label='Early')
s2 = plt.scatter([], [], marker='x', edgecolors='k', facecolors='k', label='Late')
##################

##### legend modify #####
plt.rcParams['legend.numpoints'] = 1
h, l = plt.gca().get_legend_handles_labels()
legend1 = plt.legend(h[2:], l[2:], loc='lower right', fontsize=10)  # ,framealpha=0.3)
plt.legend(h[:2], l[:2], loc='center right', fontsize=10)
plt.rcParams['legend.numpoints'] = 1
plt.gca().add_artist(legend1)
#########################

pressure = np.arange(1e2, 7e5)

############ Pressure reference
# sfr_l1 = 2.1*1e-3*(pressure/1e4)**1.18
# plt.plot((pressure),(sfr_l1),ls='--',c='k')

############ Weight reference
sfr_l2 = 1.8 * 1e-3 * (pressure / 1e4) ** 1.13
plt.plot((pressure), (sfr_l2), ls='--', c='k')
##############################

############# 1 to 1 correlation
# y = pressure
# plt.plot(pressure,y,'k--')
##################################

plt.xscale('log', nonposx='clip')
plt.yscale('log', nonposy='clip')
plt.xlim(1e3, 1e5)
# plt.ylim(1e3,1e5) # W vs P
plt.ylim(10 ** (-4), 10 ** (-2))

plt.ylabel('log($\Sigma$$_{SFR}$) [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]', fontsize=14)
# plt.xlabel('log(P/k$_b$) [K cm$^{-3}$]', fontsize=14)
plt.xlabel(r'log(W) [K cm$^{-3}$]', fontsize=14)

# plt.ylabel('log(P/k$_B$) [K cm$^{-3}$]', fontsize=14) # W vs P
# plt.ylabel(r'log(W) [K cm$^{-3}$]', fontsize=14) # W vs P

plt.tight_layout()

# plt.savefig('D:/yeongu/plots/paperplot/new/W-SFR_new.png',dpi=400)
# plt.savefig('D:/yeongu/plots/paperplot/new/W-SRF_new.eps',format='eps',dpi=400)
plt.savefig('W-SFR_new.png', dpi=400)
plt.savefig('W-SRF_new.eps', format='eps', dpi=400)
plt.show()
'''


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
C = ('gray', 'lightskyblue', 'dodgerblue', 'mediumblue', 'goldenrod', 'salmon', 'firebrick')
Model = [0, 8.63 * 1e3, 3.46 * 1e4, 3.46 * 1e4, 6.92 * 1e4, 6.92 * 1e4, 1.38 * 1e5]
crit = 94
k = 1

plt.figure(figsize=(6,16))

for j in (0,1,2,3,4,5,6):
    #plt.figure(figsize=(6, 8))
    if j==3 or j==5:
        z=range(1791)
    else:
        z = range(895)
    variable = 'Density'
    P_j = np.genfromtxt('./proj/Pism_%s.txt' % labell[j])
    P = np.genfromtxt('./proj/Pism_%s.txt' % labell[j])
    SFR = np.genfromtxt('./proj/SFR_%s.txt' % labell[j])
    W = np.genfromtxt('./proj/Weight_%s.txt' % labell[j])
    Wext = np.genfromtxt('./proj/WExt_%s.txt' % labell[j])
    Wsg = np.genfromtxt('./proj/WSg_%s.txt' % labell[j])
    # Density is calculated by np.mean(d*unit['density']) twice which is corresponding to np.sum(np.sum(d))/Area (Area=1024*1024/8/8)
    # Unit: solar mass per cubic parsec
    Dc = np.genfromtxt('./proj/Dc_%s.txt' % labell[j])
    Dh = np.genfromtxt('./proj/Dh_%s.txt' % labell[j])
    #Dc=[]
    #Dh=[]
    Dism = Dc+Dh

    # print P[crit::,:].shape
    P_e = P[0:crit, :]; P_l = P[crit::, :]
    SFR_e = SFR[0:crit]; SFR_l = SFR[crit::]
    W_e = W[0:crit, :]; W_l = W[crit::, :]
    D_e = Dc[0:crit, :]; D_l = Dc[crit::, :]
    #D_e=[]; D_l=[]
    low= 10; high = 90
    if j!=0:
        #Dicm = np.genfromtxt('./proj/Dicm_%s.txt' % labell[j])

        icm_tot = np.genfromtxt('./proj/Picm_tot_%s.txt' % labell[j])
        icm_th = np.genfromtxt('./proj/Picm_ther_%s.txt' % labell[j])
        icm_tu = np.genfromtxt('./proj/Picm_turb_%s.txt' % labell[j])
        # print D.shape
        P_tot = P+icm_tot[:,0:-1]
        #D_icm_e = Dicm[0:crit, :]; D_icm_l = Dicm[crit::,:]
        icmtot_e = icm_tot[0:crit, :]; icmtot_l = icm_tot[crit::, :]
        icmth_e = icm_th[0:crit, :]; icmth_l = icm_th[crit::, :]
        icmtu_e = icm_tu[0:crit, :]; icmtu_l = icm_tu[crit::, :]

        #### ICM properties ####
        #Dicm_medi_e = np.nanmedian(D_icm_e, axis=0) ; Dicm_medi_l = np.nanmedian(D_icm_l, axis=0)
        #Dicm_min_e = np.nanpercentile(D_icm_e, low, axis=0); Dicm_min_l = np.nanpercentile(D_icm_l, low, axis=0)
        #Dicm_max_e = np.nanpercentile(D_icm_e, high, axis=0); Dicm_max_l = np.nanpercentile(D_icm_l, high, axis=0)

        icmtot_medi_e = np.nanmedian(icmtot_e,axis=0); icmtot_medi_l = np.nanmedian(icmtot_l,axis=0)
        icmtot_min_e = np.nanpercentile(icmtot_e,low,axis=0) ; icmtot_min_l = np.nanpercentile(icmtot_l,low,axis=0)
        icmtot_max_e = np.nanpercentile(icmtot_e, high, axis=0) ; icmtot_max_l = np.nanpercentile(icmtot_l,high,axis=0)

    #### ISM properties #####
    P_medi_e = np.nanmedian(P_e, axis=0); P_medi_l = np.nanmedian(P_l, axis=0)
    P_min_e = np.nanpercentile(P_e, low, axis=0); P_min_l = np.nanpercentile(P_l, low, axis=0)
    P_max_e = np.nanpercentile(P_e, high,axis=0); P_max_l = np.nanpercentile(P_l, high, axis=0)

    W_medi_e = np.nanmedian(W_e, axis=0); W_medi_l = np.nanmedian(W_l, axis=0)
    W_min_e = np.nanpercentile(W_e, low, axis=0); W_min_l = np.nanpercentile(W_l, low, axis=0)
    W_max_e = np.nanpercentile(W_e, high,axis=0); W_max_l = np.nanpercentile(W_l, high, axis=0)
    W_ext = np.nanmedian(Wext, axis=0) ; W_sg = np.nanmedian(Wsg, axis=0) ;

    D_medi_e = np.nanmedian(D_e, axis=0) ; D_medi_l = np.nanmedian(D_l, axis=0)
    D_min_e = np.nanpercentile(D_e, low, axis=0); D_min_l = np.nanpercentile(D_l, low, axis=0)
    D_max_e = np.nanpercentile(D_e, high,axis=0); D_max_l = np.nanpercentile(D_l, high, axis=0)

    plt.subplot(7,1,k)

    #plt.plot(z, (W_medi_e+W_medi_l)/2., 'm-', label='W_all_%s' % labell[j])
    #plt.plot(z, W_ext, 'r-', label='W_ext_%s' % labell[j])
    #plt.plot(z, W_sg, 'b-', label='W_sg_%s' % labell[j])

    plt.plot(z, P_medi_e, c='c', ls='-', label='P_Early_%s' % labell[j])
    plt.plot(z, P_medi_l, c='c', ls='--', label='P_Late_%s' % labell[j])
    plt.fill_between(z,P_min_e,P_max_e,facecolor='c',alpha=0.3)
    plt.fill_between(z,P_min_l,P_max_l,facecolor='c',alpha=0.3)
    print 'Pressure',np.max(P_medi_e), np.max(P_medi_l)

    plt.plot(z,W_medi_e,'m-',label='W_Early_%s' % labell[j])
    plt.plot(z,W_medi_l,'m--',label='W_Late_%s' % labell[j])
    plt.fill_between(z,W_min_e,W_max_e,facecolor='m',alpha=0.3)
    plt.fill_between(z,W_min_l,W_max_l,facecolor='m',alpha=0.3)
    print 'Weight',np.max(W_medi_e), np.max(W_medi_l)

    plt.plot(z,D_medi_e[:-1]*1e5,'y-',label='D_Early_%s' % labell[j])
    plt.plot(z,D_medi_l[:-1]*1e5,'y--',label='D_Late_%s' % labell[j])
    #plt.fill_between(z,D_min_e[:-1]*1e5,D_max_e[:-1]*1e5,facecolor='y',alpha=0.3)
    #plt.fill_between(z,D_min_l[:-1]*1e5,D_max_l[:-1]*1e5,facecolor='y',alpha=0.3)
    if j!=0:
        plt.plot(z,icmtot_medi_e[0:-1],'k-',label='I_Early_%s' % labell[j])
        plt.plot(z,icmtot_medi_l[0:-1],'k--',label='I_Late_%s' % labell[j])

    plt.ylim(1e2,1e5)
    if j==3 or j==5:
        plt.xlim(0, 1791)
        plt.xticks([73*2, 198*2, 323*2, 448*2, 573*2, 698*2, 823*2], ['-3', '-2', '-1', '0', '1', '2', '3'])
    else:
        plt.xlim(0,895)
        plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])

    plt.yscale('log')
    plt.legend(loc='lower left')
    print j
    k = k+1
    if j==0:
        plt.title('median')
plt.xlabel('z [kpc]')
plt.tight_layout()
#plt.savefig('D:/yeongu/plots/paperplot/new/Allprof_new2.png', dpi=300)
plt.show()

'''
    ##### Early ####
    plt.subplot(3,1,1)
    #plt.plot(np.nanmedian(D_e + icmtot_e[:, 0:-1], axis=0), label='%s_Total_Early' % labell[j], c='k', ls='-',lw=3)
    plt.plot(np.nanmedian(D_e + D_icm_e, axis=0), label='%s_Total_Early' % labell[j], c='k', ls='-',lw=3)
    plt.plot(D_medi_e,label='%s_ISM_Early' % labell[j],c=C[j])
    plt.plot(Dicm_medi_e, label='%s_ICM_Early' % labell[j], c=C[j],ls='--')
    plt.yscale('log', nonposy='clip')
    plt.xlim(0, 896)
    #plt.ylim(5e2, 1e5)  # Pressure
    #plt.ylabel('Pressure [P/$k_b$]') # Pressure
    #plt.ylabel('Weight')
    plt.ylim(1e-8, 5e-2)
    plt.ylabel('Density')
    plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    plt.xlabel('z [kpc]')
    plt.legend(loc='upper right')

    ##### Late ####
    plt.subplot(3,1,2)
    #plt.plot(np.nanmedian(P_l + icmtot_l[:, 0:-1], axis=0), label='%s_Total_Late' % labell[j], c='k', ls='-',lw=3)
    plt.plot(np.nanmedian(D_l + D_icm_l, axis=0), label='%s_Total_Late' % labell[j], c='k', ls='-',lw=3)
    plt.plot(D_medi_l,label='%s_ISM_Late' % labell[j],c=C[j])
    plt.plot(Dicm_medi_l, label='%s_ICM_Late' % labell[j], c=C[j],ls='--')
    plt.yscale('log', nonposy='clip')
    plt.xlim(0, 896)
    #plt.ylim(5e2, 1e5)  # Pressure
    #plt.ylabel('Pressure [P/$k_b$]') # Pressure
    #plt.ylabel('Weight')
    plt.ylim(1e-8, 5e-2)
    plt.ylabel('Density')
    plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    plt.xlabel('z [kpc]')
    plt.legend(loc='upper right')

    ##### All time ####
    plt.subplot(3,1,3)
    plt.plot(np.nanmedian(Dc+Dicm, axis=0), label='%s_Total_All' % labell[j], c='k', ls='-',lw=3)
    plt.plot(np.nanmedian(Dc, axis=0),label='%s_ISM_All' % labell[j],c=C[j])
    plt.plot(np.nanmedian(Dicm, axis=0), label='%s_ICM_All' % labell[j], c=C[j],ls='--')
    plt.yscale('log', nonposy='clip')
    plt.xlim(0, 896)
    #plt.ylim(5e2, 1e5)  # Pressure
    #plt.ylabel('Pressure [P/$k_b$]')# Pressure
    #plt.ylabel('Weight')
    plt.ylim(1e-8,5e-2)
    plt.ylabel('Density')
    plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    plt.xlabel('z [kpc]')
    plt.legend(loc='upper right')
    plt.tight_layout()
    #plt.savefig('D:/yeongu/plots/paperplot/new/%s-%s-zprof.png' % (labell[j], variable), dpi=300)
    plt.show()
'''
'''
    icmtot_medi_e = np.nanmedian(icmtot_e,axis=0); icmtot_medi_l = np.nanmedian(icmtot_l,axis=0)
    icmtot_min_e = np.nanpercentile(icmtot_e,low,axis=0) ; icmtot_min_l = np.nanpercentile(icmtot_l,low,axis=0)
    icmtot_max_e = np.nanpercentile(icmtot_e, high, axis=0) ; icmtot_max_l = np.nanpercentile(icmtot_l,high,axis=0)

    icmth_medi_e = np.nanmedian(icmth_e,axis=0); icmth_medi_l = np.nanmedian(icmth_l,axis=0)
    icmth_min_e = np.nanpercentile(icmth_e,low,axis=0) ; icmth_min_l = np.nanpercentile(icmth_l,low,axis=0)
    icmth_max_e = np.nanpercentile(icmth_e, high, axis=0) ; icmth_max_l = np.nanpercentile(icmth_l,high,axis=0)

    icmtu_medi_e = np.nanmedian(icmtu_e,axis=0); icmtu_medi_l = np.nanmedian(icmtu_l,axis=0)
    icmtu_min_e = np.nanpercentile(icmtu_e,low,axis=0) ; icmtu_min_l = np.nanpercentile(icmtu_l,low,axis=0)
    icmtu_max_e = np.nanpercentile(icmtu_e, high, axis=0) ; icmtu_max_l = np.nanpercentile(icmtu_l,high,axis=0)
'''

#d_e_p = np.argmax(D_medi_e);
#d_l_p = np.argmax(D_medi_l)  # Density profile peak position

#print labell[j], d_e_p, d_l_p




'''
    #plt.subplot(4,1,k)
    plt.plot(z,icmtot_medi_e[0:-1],'r-',label='P$_{ICM}$_%s' % labell[j])
    #plt.fill_between(z,P_min_e,P_max_e,facecolor=C[j],alpha=0.5)
    plt.plot(z,icmtot_medi_l[0:-1],'r--')
    #plt.fill_between(z,P_min_l,P_max_l,facecolor=C[j],alpha=0.5)

    plt.plot(z,icmth_medi_e[0:-1],c='coral',ls='-',label='P$_{ICM,ther}$_%s' % labell[j])
    #plt.fill_between(z,P_min_e,P_max_e,facecolor=C[j],alpha=0.5)
    plt.plot(z,icmth_medi_l[0:-1],c='coral',ls='--')
    #plt.fill_between(z,P_min_l,P_max_l,facecolor=C[j],alpha=0.5)

    plt.plot(z,icmtu_medi_e[0:-1],c='gold',ls='-',label='P$_{ICM,turb}$_%s' % labell[j])
    #plt.fill_between(z,P_min_e,P_max_e,facecolor=C[j],alpha=0.5)
    plt.plot(z,icmtu_medi_l[0:-1],c='gold',ls='--')
    #plt.fill_between(z,P_min_l,P_max_l,facecolor=C[j],alpha=0.5)

    plt.plot(z,W_medi_e,'c-',label='W$_{ISM}$_%s' % labell[j])
    #plt.fill_between(z,W_min_e,W_max_e,facecolor=C[j],alpha=0.5)
    plt.plot(z,W_medi_l,'c--')



    plt.plot([],[],label='Early',ls='-',c='k')
    plt.axhline(Model[j], ls='-.', c='k', label='Input ICM pressure')
    plt.plot([],[],label='Late',ls='--',c='k')


    h, l = plt.gca().get_legend_handles_labels()
    legend1 = plt.legend(h[:4], l[:4], loc='upper right')  # ,fontsize=14,framealpha=0.3)
    if k==1:
        plt.legend(h[4:], l[4:], loc='upper left',ncol=2 ,fontsize=11.5)
        plt.gca().add_artist(legend1)
'''
'''
    plt.fill_between(z, P_min_e, P_max_e, facecolor=C[j], alpha=0.5)
    plt.fill_between(z, P_min_l, P_max_l, facecolor=C[j], alpha=0.5)

    plt.plot(z, P_medi_e, c=C[j], ls='-', label='P_Early_%s' % labell[j])
    plt.plot(z, P_medi_l, c=C[j], ls='--', label='P_Late_%s' % labell[j])

    # plt.plot(z,W_medi_e,'m-',label='W_Early_%s' % labell[j])
    # plt.fill_between(z,W_min_e,W_max_e,facecolor=C[j],alpha=0.5)
    # plt.plot(z,W_medi_l,'m--',label='W_Late_%s' % labell[j])
    # print np.max(D_medi_e)
    # print D_medi_e[0:-1].shape
    # print len(z)
    # plt.plot(z,D_medi_e[:-1],'y-',label='D_Early_%s' % labell[j])
    # plt.fill_between(z,W_min_e,W_max_e,facecolor=C[j],alpha=0.5)
    # plt.plot(z,D_medi_l[:-1],'y--',label='D_Late_%s' % labell[j])

    # e= np.argmax(W_medi_e)
    # l = np.argmax(W_medi_l)
    # print e, l
    # plt.axvline(e,ls='-',c='r')
    # plt.axvline(l,ls='--',c='b')
    # plt.fill_between(z,W_min_l,W_max_l,facecolor=C[j],alpha=0.5)
    plt.yscale('log', nonposy='clip')
    plt.ylim(5e2, 1e5)
    plt.xlim(0, 896)
    plt.ylabel('Pressure [P/$k_b$]')
    # if k==4:
    plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    plt.xlabel('z [kpc]')
    # else:
    #    plt.xticks([73, 198, 323, 448, 573, 698, 823], [])
    plt.tick_params(which='both', direction='in')
    plt.legend(loc='upper right')

    k = k + 1
    plt.tight_layout()
    plt.savefig('D:/yeongu/plots/paperplot/new/%s-%s-zprof.png' % (labell[j], variable), dpi=300)
    # plt.savefig('D:/yeongu/plots/paperplot/new/WismPicm.eps',format='eps',dpi=300)
    plt.show()
    plt.close()
'''
'''
    #print np.argmax(W_medi_e)
    #print np.argmax(W_medi_l)
    #if j!=1:
    pe = P_medi_e[d_e_p]; pl = P_medi_l[d_l_p]
    we = W_medi_e[d_e_p]; wl = W_medi_l[d_l_p]
    #else:
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

    plt.errorbar(we, SFR_medi_e, yerr=np.array([SFR_e_errors]).T, xerr=np.array([W_e_errors]).T,
                 marker='o', mfc=C[j], mec=C[j], capsize=3.5, label=labell[j],ecolor=C[j])
    plt.errorbar(wl, SFR_medi_l, yerr=np.array([SFR_l_errors]).T, xerr=np.array([W_l_errors]).T,
                 marker='s', mfc='none', mec=C[j], capsize=3.5,ecolor=C[j])


    plt.errorbar(we, pe, yerr=np.array([P_e_errors]).T, xerr=np.array([W_e_errors]).T,
                 marker='o', mfc=C[j], mec=C[j], capsize=3.5, label=labell[j],ecolor=C[j])
    plt.errorbar(wl, pl, yerr=np.array([P_l_errors]).T, xerr=np.array([W_l_errors]).T,
                 marker='s', mfc='none', mec=C[j], capsize=3.5,ecolor=C[j])
    '''
'''
##### legend #####
s1 = plt.scatter([], [], marker='o',edgecolors='k',facecolors='k', alpha=0.7,label='Early')
s2 = plt.scatter([], [], marker='s',edgecolors='k',facecolors='none', alpha=0.7,label='Late')
##################

##### legend modify #####
h, l =plt.gca().get_legend_handles_labels()
legend1=plt.legend(h[2:],l[2:],loc='upper left')#,fontsize=14,framealpha=0.3)
plt.legend(h[:2],l[:2],loc='center left')#,fontsize=12)
plt.gca().add_artist(legend1)
#########################

pressure=np.arange(1e2,7e5)

############ Pressure reference
#sfr_l1 = 2.1*1e-3*(pressure/1e4)**1.18
#plt.plot((pressure),(sfr_l1),ls='--',c='k')

############ Weight reference
sfr_l2 = 1.8*1e-3*(pressure/1e4)**1.13
plt.plot((pressure),(sfr_l2),ls='--',c='k')
##############################
y = pressure
#plt.plot(pressure,y,'k--')
plt.xscale('log', nonposx='clip')
plt.yscale('log', nonposy='clip')
plt.xlim(3e2,3e5)
#plt.ylim(3e2,3e5) # W vs P
plt.ylim(10**(-4.5),10**(-1.5))
plt.ylabel('log(SFR)')
#plt.ylabel('log(P/kb)') # W vs P
#plt.xlabel('log(P/kb)')
plt.xlabel('log(Weight)')
plt.show()
'''
'''
#print P_medi_e.shape
#plt.figure(figsize=(6,3))
#plt.plot(P_medi_e,'c-')
#plt.plot(P_medi_l,'m--')
#plt.yscale('log', nonposy='clip')
#plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
#plt.ylim(1e2,1e5)
#plt.xlim(0,895)
#plt.show()
'''
'''
if min(SFR_e)==0 or min(SFR_l)==0 or min(P_e)==0 or min(P_l)==0:
    sfr_e_min=0
    sfr_l_min=0
    p_e_min=0
    p_l_min=0
else:
    sfr_e_min=np.log10(np.min(SFR_e))
    sfr_l_min=np.log10(np.min(SFR_l))
    p_e_min=np.log10(np.min(P_e))
    p_l_min=np.log10(np.min(P_l))

############################### Upper and Lower #########################################
SFR_e_errors = [np.mean(SFR_e)-np.min(SFR_e),np.max(SFR_e)-np.mean(SFR_e)] # [np.nonzero(SFR_e)] nonzero selection
print SFR_e_errors
SFR_l_errors = [np.mean(SFR_l)-np.min(SFR_l),np.max(SFR_l)-np.mean(SFR_l)]
print SFR_l_errors
P_e_errors = [np.mean(P_e)-np.min(P_e),np.max(P_e)-np.mean(P_e)]
print np.array([P_e_errors]).T
P_l_errors = [np.mean(P_l)-np.min(P_l),np.max(P_l)-np.mean(P_l)]
print np.array([P_l_errors]).T
#########################################################################################

plt.errorbar(np.mean(P_e),np.mean(SFR_e),yerr=np.array([SFR_e_errors]).T,xerr=np.array([P_e_errors]).T,ecolor=C[j],marker='o',mfc=C[j],mec=C[j],capsize=3.5,label=labell[j])
plt.errorbar(np.mean(P_l),np.mean(SFR_l),yerr=np.array([SFR_l_errors]).T,xerr=np.array([P_l_errors]).T,ecolor=C[j],marker='s',mfc='none',mec=C[j],capsize=3.5)

##### legend #####
s1 = plt.scatter([], [], marker='o',edgecolors='k',facecolors='k', alpha=0.7,label='Early')
s2 = plt.scatter([], [], marker='s',edgecolors='k',facecolors='none', alpha=0.7,label='Late')
##################

##### theoretical line #####
pressure=np.arange(1e2,7e5)
sfr_l = 2.1*1e-3*(pressure/1e4)**1.18
plt.plot((pressure),(sfr_l),ls='--',c='k')
############################

##### legend modify #####
h, l =plt.gca().get_legend_handles_labels()
legend1=plt.legend(h[2:],l[2:],loc='upper left')#,fontsize=14,framealpha=0.3)
plt.legend(h[:2],l[:2],loc='center left')#,fontsize=12)
plt.gca().add_artist(legend1)
#########################

plt.xscale('log', nonposx='clip')
plt.yscale('log', nonposy='clip')
#plt.xlim(3e2,3e5)
plt.ylim(10**(-4.5),10**(-1.5))
plt.ylabel('log(SFR)')
plt.xlabel('log(Weight)')
#plt.legend(loc=0)
plt.tight_layout()
#plt.savefig('D:/yeongu/plots/paperplot/new/SFR_P_wz.png',dpi=300)
#plt.savefig('D:/yeongu/plots/paperplot/new/SFR_P_wz.eps',format='eps',dpi=300)
#plt.show()
'''
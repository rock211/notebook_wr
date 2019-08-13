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
labell = ('No ICM','P1', 'P3','P3h', 'P7','P7h','P14')  # r'No ICM',
#C = ('gray', 'mediumturquoise', 'dodgerblue','mediumblue' ,'goldenrod','salmon', 'firebrick','darkmagenta','goldenrod','royalblue','crimson') # 'plum','orchid','purple'
C = ('k', 'powderblue', 'dodgerblue','mediumblue' ,'salmon', 'crimson','goldenrod')
Model = [0,8.63*1e3,3.46*1e4,3.46*1e4,6.92*1e4,6.92*1e4,1.38*1e5]
crit = 94
k=1
#plt.figure(figsize=(6,10))
plt.figure(figsize=(5,5))
sf_ind = ([446, 453],[456, 445],[460, 456],[915, 914],[469, 553],[965, 1052],[506, 684])
sf_ind = np.array(sf_ind)
print sf_ind.shape
for j in (0,1,2,3):
    #plt.figure(figsize=(6,4))
    variable = 'density'

    P = np.genfromtxt('./proj/Pism_%s.txt' % labell[j])
    SFR = np.genfromtxt('./proj/SFR_%s.txt'% labell[j])
    W = np.genfromtxt('./proj/Weight_%s.txt'% labell[j])
    D = np.genfromtxt('./proj/Dc_%s.txt' % labell[j])
    if j!=0:
        icm_tot = np.genfromtxt('./proj/Picm_tot_%s.txt' % labell[j])
        icm_th = np.genfromtxt('./proj/Picm_ther_%s.txt' % labell[j])
        icm_tu = np.genfromtxt('./proj/Picm_turb_%s.txt' % labell[j])
    #print D.shape
    #for t in range(50):
    #    plt.semilogy(P[t,:],c='k')
        #plt.axhline(np.max(Tot),ls='--')
    #    print np.max(P[t,:])
        #plt.ylim(1e2,1e5)
        #plt.show()
        #plt.close()
    #print P[crit::,:].shape
    P_e = P[0:crit,:]; P_l = P[crit::,:] # +icm_tot[0:crit,0:-1] +icm_tot[crit::,0:-1]
    SFR_e = SFR[0:crit]; SFR_l = SFR[crit::]
    W_e = W[0:crit,:]; W_l = W[crit::,:]
    D_e = D[0:crit, :]; D_l = D[crit::, :]

    P_tot_medi = np.nanmean(P,axis=0)

    low= 25; high = 75

    P_medi_e = np.nanmean(P_e,axis=0); P_medi_l = np.nanmean(P_l,axis=0)
    P_min_e = np.nanpercentile(P_e,low,axis=0) ; P_min_l = np.nanpercentile(P_l,low,axis=0)
    P_max_e = np.nanpercentile(P_e,high,axis=0) ; P_max_l = np.nanpercentile(P_l,high,axis=0)

    SFR_medi_e = np.nanmean(SFR_e); SFR_medi_l = np.nanmean(SFR_l)
    SFR_min_e = np.nanpercentile(SFR_e,low) ; SFR_min_l = np.nanpercentile(SFR_l,low)
    SFR_max_e = np.nanpercentile(SFR_e, high) ; SFR_max_l = np.nanpercentile(SFR_l,high)

    W_medi_e = np.nanmean(W_e,axis=0); W_medi_l = np.nanmean(W_l,axis=0)
    W_min_e = np.nanpercentile(W_e,low,axis=0) ; W_min_l = np.nanpercentile(W_l,low,axis=0)
    W_max_e = np.nanpercentile(W_e, high, axis=0) ; W_max_l = np.nanpercentile(W_l,high,axis=0)
    print np.nanmax(W_e)
    print np.nanmax(W_l)
    D_medi_e = np.nanmean(D_e,axis=0)*1e6; D_medi_l = np.nanmean(D_l,axis=0)*1e6
    D_min_e = np.nanpercentile(D_e,low,axis=0) ; D_min_l = np.nanpercentile(D_l,low,axis=0)
    D_max_e = np.nanpercentile(D_e, high, axis=0) ; D_max_l = np.nanpercentile(D_l,high,axis=0)

    # ICM calculation below
    '''
    icmtot_e = icm_tot[0:crit, :]; icmtot_l = icm_tot[crit::, :]
    icmth_e = icm_th[0:crit, :]; icmth_l = icm_th[crit::, :]
    icmtu_e = icm_tu[0:crit, :]; icmtu_l = icm_tu[crit::, :]
    
    icmtot_medi_e = np.nanmean(icmtot_e,axis=0); icmtot_medi_l = np.nanmean(icmtot_l,axis=0)
    icmtot_min_e = np.nanpercentile(icmtot_e,25,axis=0) ; icmtot_min_l = np.nanpercentile(icmtot_l,25,axis=0)
    icmtot_max_e = np.nanpercentile(icmtot_e, 75, axis=0) ; icmtot_max_l = np.nanpercentile(icmtot_l,75,axis=0)

    icmth_medi_e = np.nanmean(icmth_e,axis=0); icmth_medi_l = np.nanmean(icmth_l,axis=0)
    icmth_min_e = np.nanpercentile(icmth_e,25,axis=0) ; icmth_min_l = np.nanpercentile(icmth_l,25,axis=0)
    icmth_max_e = np.nanpercentile(icmth_e, 75, axis=0) ; icmth_max_l = np.nanpercentile(icmth_l,75,axis=0)

    icmtu_medi_e = np.nanmean(icmtu_e,axis=0); icmtu_medi_l = np.nanmean(icmtu_l,axis=0)
    icmtu_min_e = np.nanpercentile(icmtu_e,25,axis=0) ; icmtu_min_l = np.nanpercentile(icmtu_l,25,axis=0)
    icmtu_max_e = np.nanpercentile(icmtu_e, 75, axis=0) ; icmtu_max_l = np.nanpercentile(icmtu_l,75,axis=0)
    '''
    # Slab peak choose
    '''
    if j == 3 or j == 5:
        l_b = 75  # box length corresponding to 300pc
        cen = 895
        z = range(1792)
    else:
        l_b = 38
        cen = 447
        z = range(896)
    # print D_medi_e,D_medi_e.shape
    box_max_e = []
    box_max_l = []
    for i in range(len(D_medi_e) - l_b):
        a = D_medi_e[i:i + l_b]
        b = D_medi_l[i:i + l_b]
        box_max_e.append(np.mean(a))
        box_max_l.append(np.mean(b))

    d_peak_e_ind = np.argmax(D_medi_e)
    box_p_e = np.argmax(box_max_e) + l_b / 2

    d_peak_l_ind = np.argmax(D_medi_l)
    box_p_l = np.argmax(box_max_l) + l_b / 2
    '''

    #d_e_p = np.argmax(D_medi_e); d_l_p = np.argmax(D_medi_l) # Density profile peak position
    #d_e_p = np.argmax(W_medi_e); d_l_p = np.argmax(W_medi_l)  # Weight profile peak position
    #d_e_p = np.argmax(P_medi_e); d_l_p = np.argmax(P_medi_l)  # Pressure profile peak position
    d_e_p = np.argmax(P_tot_medi); d_l_p = np.argmax(P_tot_medi)  # Pressure profile peak position
    #d_e_p = box_p_e ; d_l_p = box_p_l # Box peak
    #d_e_p = sf_ind[j,0]; d_l_p = sf_ind[j,1] # SF peak
    #if j==3 or j==5:
    #    d_e_p = 896; d_l_p = 896
    #else:
    #    d_e_p = 448; d_l_p = 448
    print labell[j], d_e_p, d_l_p
    #print P_medi_e.shape
    #print SFR_medi_e
    #print W_medi_e.shape
    #print D_medi_e.shape

    #print D_medi_e[0:-1]
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
    

    plt.fill_between(z,P_min_e,P_max_e,facecolor=C[j],alpha=0.5)
    plt.fill_between(z,P_min_l,P_max_l,facecolor=C[j],alpha=0.5)

    plt.plot(z, P_medi_e, c=C[j], ls='-', label='P_Early_%s' % labell[j])
    plt.plot(z, P_medi_l, c=C[j], ls='--', label='P_Late_%s' % labell[j])

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
    plt.xlim(0, 896)
    plt.ylabel('Pressure [P/$k_b$]')
    #if k==4:
    plt.xticks([73, 198, 323, 448, 573, 698, 823], ['-3', '-2', '-1', '0', '1', '2', '3'])
    plt.xlabel('z [kpc]')
    #else:
    #    plt.xticks([73, 198, 323, 448, 573, 698, 823], [])
    plt.tick_params(which='both', direction='in')
    plt.legend(loc='upper right')

    k = k+1
    plt.tight_layout()
    #plt.savefig('D:/yeongu/plots/paperplot/new/%s-%s-zprof.png' % (labell[j],variable),dpi=300)
    #plt.savefig('D:/yeongu/plots/paperplot/new/WismPicm.eps',format='eps',dpi=300)
    plt.show()
    plt.close()
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
    #print P_e_errors,P_l_errors
    plt.errorbar(we, SFR_medi_e,  xerr=np.array([W_e_errors]).T, yerr=np.array([SFR_e_errors]).T,color='none',ecolor=C[j], marker='o', mfc=C[j],mec=C[j], capsize=3.5, label=labell[j])
    if j != 5 and j != 4 and j != 6:
        plt.errorbar(wl, SFR_medi_l, xerr=np.array([W_l_errors]).T, yerr=np.array([SFR_l_errors]).T, color='none',ecolor=C[j], marker='x', mfc='none', mec=C[j], capsize=3.5,ms=10)


    
    #plt.errorbar(we, pe,  xerr=np.array([W_e_errors]).T, yerr=np.array([P_e_errors]).T,color='none',ecolor=C[j], marker='o', mfc=C[j],mec=C[j], capsize=3.5, label=labell[j])
    #if j != 5 and j!=4 and j!=6:
    #    plt.errorbar(wl, pl,  xerr=np.array([W_l_errors]).T, yerr=np.array([P_l_errors]).T, color='none',ecolor=C[j], marker='x', mfc='none', mec=C[j], capsize=3.5,ms=10)


##### legend #####
s1 = plt.scatter([], [], marker='o',edgecolors='k',facecolors='k',label='Early')
s2 = plt.scatter([], [], marker='x',edgecolors='k',facecolors='k',label='Late')
##################

##### legend modify #####
plt.rcParams['legend.numpoints']=1
h, l =plt.gca().get_legend_handles_labels()
legend1=plt.legend(h[2:],l[2:],loc='lower right',fontsize=10)#,framealpha=0.3)
plt.legend(h[:2],l[:2],loc='lower center',fontsize=10)
plt.rcParams['legend.numpoints']=1
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

############# 1 to 1 correlation
#y = pressure
#plt.plot(pressure,y,'k--')
##################################

plt.xscale('log', nonposx='clip')
plt.yscale('log', nonposy='clip')

plt.xlim(3e3,7e4) # Weight xaxis
plt.xticks([3e3,1e4,3e4,7e4],[r'3$\times$$10^3$',r'$10^4$',r'3$\times$$10^4$',r'7$\times$$10^4$'],fontsize=12)
plt.xlabel(r'log(W) [K cm$^{-3}$]', fontsize=14)
#plt.xlabel('log(P/k$_b$) [K cm$^{-3}$]', fontsize=14)

#plt.ylim(3e3,7e4) # W vs P
#plt.yticks([3e3,1e4,3e4,7e4],[r'3$\times$$10^3$',r'$10^4$',r'3$\times$$10^4$',r'7$\times$$10^4$'],fontsize=12)
#plt.ylabel('log(P/k$_B$) [K cm$^{-3}$]', fontsize=14) # W vs P

plt.ylim(4e-4,1e-2) # W vs SFR
plt.yticks([1e-3,1e-2],[r'$10^{-3}$',r'$10^{-2}$'],fontsize=12)
plt.ylabel('log($\Sigma$$_{SFR}$) [M$_\odot$ yr$^{-1}$ kpc$^{-2}$]', fontsize=14)

#plt.title('Pressure peak, mean')
plt.tight_layout()

#plt.savefig('D:/yeongu/plots/paperplot/new/W-SFR_new.png',dpi=400)
#plt.savefig('D:/yeongu/plots/paperplot/new/W-SRF_new.eps',format='eps',dpi=400)
plt.savefig('W-P_new_mean_origin.png',dpi=300)
plt.savefig('W-P_new_mean_origin.eps',format='eps',dpi=300)
plt.show()


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

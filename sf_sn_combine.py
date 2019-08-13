import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
from scipy.interpolate import spline
from matplotlib.ticker import MultipleLocator
from pyathena import ath_hst
import pyathena as pa

unit = pa.set_units(muH=1.4271)
Msun = unit['mass'].to('Msun').value
print Msun
Myr=unit['time'].to('Myr').value
agebin = 10

simid_t = ('R8_8pc_metal','RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1','RPS_8pc_ICM2','RPS_8pc_ICM3','RPS_8pc_ICM4') # 'MHD_8pc_new' ,
labell = ('No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # r'No ICM',
C = ('k','r','g','b','cyan','magenta','y') #'darkkhaki','royalblue','firebrick'

S = ('-','-','-')
L = (1.5,1.5,1.5)
alpha = (0.8,1,1)
# overplot Starformation rate of three different simulations
start = 250
stop = 500
plt.figure(figsize=(5.5,5))

for j in range(len(simid_t)) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j]
    Mj = []
    M_star = 0.0
    Mass_tot = []
    SFE = []
    if j == 5:
        stop = 472
    elif j == 6:
        stop = 401
    else:
        stop = 500
    id = simid_t[j]
    ########### SN data #############
    hstfilename = basedir + id + '/hst/' + id + '.hst'
    hstp = ath_hst.read_w_pandas(hstfilename)
    sn=ath_hst.read_w_pandas(hstfilename.replace('.hst','.sn'))
    #print sn['time']

    for tidx in range(250, stop):
        id = simid_t[j]  # 'MHD_8pc_new' , 'RPS_8pc_n1e-4_v1414' , 'RPS_8pc_n2e-4_v1414'
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)

        ########## star data ###########
        starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)
        #print(sp.shape)
        # print (sp)

        star_clu = sp[sp['mass'] != 0]  # choose cluster, if choose 'mass' = 0, it means runaway star
        star_clu2 = star_clu[star_clu['age'] * Myr < agebin]
        M_star = sum(star_clu2['mass']) * Msun # SFR time variation

        Mj.append(M_star* unit['time'].value / (1e+6*agebin*1.024*1.024))

    time = np.arange(250 * unit['time'].value, stop * unit['time'].value, unit['time'].value)

    plt.subplot(211)
    plt.plot(time, Mj, label='%s' % labell[j], color=C[j])
    plt.xticks([])
    # plt.ylabel(r'SFE',fontsize=17) # For SFE
    plt.ylabel(r'$\Sigma_{SFR}$ [M$_{\odot}$ kpc$^{-2}$ yr$^{-1}$]',fontsize=10) # For SFR
    plt.xlim(time.min(), time.max())
    plt.ylim(0,0.019)
    #plt.ylim(0, 1.1)  # For cumulative
    plt.legend(loc=0, fontsize=10)

    time2 = np.arange(250 * unit['time'].value, stop * unit['time'].value, 10*unit['time'].value)
    plt.subplot(212)
    hist, binedge = np.histogram(sn['time'],bins=time2)
    h = hist/float(np.sum(hist))
    x_c = (binedge[1::] + binedge[0:-1]) / 2
    plt.plot(x_c,h,color=C[j])
    plt.xlim(time.min(), time.max())
    plt.ylim(0, 0.22)
    #plt.ylim(0,0.023)
    plt.xlabel(r'Time [Myr]', fontsize=13)
    plt.ylabel('SN fraction')
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('D:/yeongu/plots/SFR_SN_10_%s.png' % labell[j], dpi=200)
    #plt.show()
    plt.close()
    print j
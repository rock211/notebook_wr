import matplotlib.pyplot as plt
import os
from shutil import copyfile
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
from scipy.stats import lognorm
import cPickle
from mpl_toolkits import axes_grid1


from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle
import time
# In[7]:
import matplotlib as mpl


# In[8]:

# ## Unit system
#
# The unit system we choose for this simulation is
# * [length] = pc
# * [velocity] = km/s
# * [density] = 1.4271*m_h/cm^3

unit = pa.set_units(muH=1.4271)
print(unit)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
#print unit['density']
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K
#vpc = 7168.*1024*1024/(128*128*896) # volume per cell

# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])

simid_t = ('R8_8pc_metal','RPS_8pc_ICM00','RPS_8pc_ICM0','RPS_8pc_ICM1','RPS_8pc_ICM2','RPS_8pc_ICM3','RPS_8pc_ICM4') # 'R8_8pc_metal'
labell = (r'No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # 'NonICM',r'No ICM',
C = ('dimgrey','coral','royalblue')
S = ('-.','--','-')

#start_time = time.time()
for j in (1,2,6) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j] # simid_t[j]

    if j == 5:
        stop = 472
    elif j == 6:
        stop = 401
    else:
        stop = 500
    print simid

    H = [] ; C = [] ; I = []

    for tidx in range(250,stop):

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        ds = pa.AthenaDataSet(vtkfname)

        #print(ds.field_list) ; print(ds.derived_field_list)

        d = ds.read_all_data('density')*unit['density'].value # density : solMass/pc3
        vel = ds.read_all_data('velocity')
        velz = vel[:, :, :, 2]
        T1 = ds.read_all_data('T1') ; coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling function

        if j!=0:
            scalar = ds.read_all_data('specific_scalar4')  # ism = 0 / icm = 1
        #surf_d = np.sum(d[scalar < 0.7] * (8 * 8 * 8) / (1024 * 1024.))
        #surf_d = np.sum(d* (8 * 8 * 8) / (1024 * 1024.))
        #print ('%s sec' %(time.time()-start_time))

        bin = np.arange(-1000, 1000, 5)
        #total = []
        cold = [] ; hot = [] ; ion = []
        ############ make histogram ##############
        for height in range(0,896):
            dh = d[height] ; velzh = velz[height]; temph = temp[height];
            if j!=0:
                scalarh = scalar[height]; dh = dh[scalarh < 0.7] ; velzh = velzh[scalarh < 0.7]; temph=temph[scalarh < 0.7] # choose only ism

            hist_c, binedges_c = np.histogram(np.ravel(velzh[temph < 20000]), bins=bin,weights=np.ravel(dh[temph < 20000]))  # ,density=True)
            hist_h, binedges_h = np.histogram(np.ravel(velzh[(temph > 20000) & (temph < 500000)]), bins=bin,weights=np.ravel(dh[(temph > 20000) & (temph < 500000)]))
            hist_i, binedges_i = np.histogram(np.ravel(velzh[temph > 500000]), bins=bin,weights=np.ravel(dh[temph > 500000]))
            #x_c0 = (binedges[1::] + binedges[0:-1]) / 2
            cold.append(hist_c)
            hot.append(hist_h)
            ion.append(hist_i)
        H.append(hot)
        C.append(cold)
        I.append(ion)

        print len(H), len(C)


    ########## data save ##############
        #a =np.stack([hist0,hist0_c,hist0_h,histh,histh_c,histh_h,hist1,hist1_c,hist1_h,hist2,hist2_c,hist2_h,hist3,hist3_c,hist3_h])

    C = np.array(C)
    H = np.array(H)
    I = np.array(I)
    print C.shape, H.shape
    #print b
    cPickle.dump(C,open('surf_ism_cold_%s.pkl' % labell[j],'wb'))
    cPickle.dump(H, open('surf_ism_hot_%s.pkl' % labell[j], 'wb'))
    cPickle.dump(I, open('surf_ism_ion_%s.pkl' % labell[j], 'wb'))
    #dd = cPickle.load(open('surf_%s.pkl' % labell[j],'rb'))
    #print dd
    ####################################
    '''
        s0 = float(np.sum(hist0)) ; sh = float(np.sum(histh)) ; s1 = float(np.sum(hist1)) ; s2 = float(np.sum(hist2)) ; s3 = float(np.sum(hist3))

        ##### histogram plot ##########

        if s0!=0:
            plt.semilogy(x_c0_c, hist0_c / float(np.sum(hist0)), label='Cold/Warm',c='royalblue')
            plt.semilogy(x_c0_h, hist0_h / float(np.sum(hist0)), label='Hot',c='crimson')
            plt.bar(x_c0, hist0 / float(np.sum(hist0)), align='center', linewidth=0.2, width=5, alpha=0.7, log=1,facecolor='k',label='total')
            plt.title('%s_0 kpc_%s' % (labell[j],round(tidx*unit['time'].value)))
            plt.legend(loc=1)
            plt.ylim(1e-5,1)
            #plt.savefig('D:/yeongu/plots/surf_d/Hist_%s_0kpc_%s.png' % (labell[j], tidx), dpi=100)
            #plt.show()
            plt.close()

        if sh!=0:
            plt.semilogy(x_ch_c, histh_c / float(np.sum(histh)), label='Cold/Warm',c='royalblue')
            plt.semilogy(x_ch_h, histh_h / float(np.sum(histh)), label='Hot',c='crimson')
            plt.bar(x_ch, histh / float(np.sum(histh)), align='center', linewidth=0.2, width=5, alpha=0.7, log=1,facecolor='k',label='total')
            plt.title('%s_0.5 kpc_%s' % (labell[j],round(tidx*unit['time'].value)))
            plt.legend(loc=1)
            plt.ylim(1e-5,1)
            #plt.savefig('D:/yeongu/plots/surf_d/Hist_%s_0.5kpc_%s.png' % (labell[j], tidx), dpi=100)
            #plt.show()
            plt.close()
        if s1!=0:
            plt.semilogy(x_c1_c, hist1_c / float(np.sum(hist1)), label='Cold/Warm',c='royalblue')
            plt.semilogy(x_c1_h, hist1_h / float(np.sum(hist1)), label='Hot',c='crimson')
            plt.bar(x_c1, hist1 / float(np.sum(hist1)), align='center', linewidth=0.2, width=5, alpha=0.7, log=1,facecolor='k',label='total')
            plt.title('%s_1 kpc_%s' % (labell[j],round(tidx*unit['time'].value)))
            plt.legend(loc=1)
            plt.ylim(1e-5,1)
            #plt.savefig('D:/yeongu/plots/surf_d/Hist_%s_1kpc_%s.png' % (labell[j], tidx), dpi=100)
            #plt.show()
            plt.close()

        plt.semilogy(x_c2_c, hist2_c / float(np.sum(hist2)), label='Cold/Warm',c='royalblue')
        plt.semilogy(x_c2_h, hist2_h / float(np.sum(hist2)), label='Hot',c='crimson')
        plt.bar(x_c2, hist2 / float(np.sum(hist2)), align='center', linewidth=0.2, width=5, alpha=0.7, log=1,facecolor='k',label='total')
        plt.title('%s_2 kpc_%s' % (labell[j],round(tidx*unit['time'].value)))
        plt.legend(loc=1)
        plt.ylim(1e-5,1)
        #plt.savefig('D:/yeongu/plots/surf_d/Hist_%s_2kpc_%s.png' % (labell[j], tidx), dpi=100)
        #plt.show()
        plt.close()

        plt.semilogy(x_c3_c, hist3_c / float(np.sum(hist3)), label='Cold/Warm',c='royalblue')
        plt.semilogy(x_c3_h, hist3_h / float(np.sum(hist3)), label='Hot',c='crimson')
        plt.bar(x_c3, hist3 / float(np.sum(hist3)), align='center', linewidth=0.2, width=5, alpha=0.7, log=1,facecolor='k',label='total')
        plt.title('%s_3 kpc_%s' % (labell[j],round(tidx*unit['time'].value)))
        plt.legend(loc=1)
        plt.ylim(1e-5,1)
        #plt.savefig('D:/yeongu/plots/surf_d/Hist_%s_3kpc_%s.png' % (labell[j], tidx), dpi=100)
        #plt.show()
        plt.close()

        ######### total height ########### hist over surface density is not clear

        plt.title('Total_%s_%s' % (labell[j],tidx))
        plt.semilogy(x_c0, hist0 / surf_d, label='0 kpc', c='k')
        plt.semilogy(x_ch, histh / surf_d, label='0.5 kpc', c='navy',alpha=0.8)
        plt.semilogy(x_c1, hist1 / surf_d, label='1 kpc', c='blue',alpha=0.9)
        plt.semilogy(x_c2, hist2 / surf_d, label='2 kpc', c='deepskyblue')
        plt.semilogy(x_c3, hist3 / surf_d, label='3 kpc', c='skyblue')
        plt.legend(loc=1)
        plt.ylim(1e-8,1e-2)
        plt.savefig('D:/yeongu/plots/surf_d_test/Hist_%s_total_%s.png' % (labell[j], tidx), dpi=100)
        #plt.show()
        plt.close()

        print tidx
    '''
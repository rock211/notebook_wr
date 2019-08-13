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
for j in range(1,3) :
    basedir = 'F:/yeongu/'
    simid = simid_t[j] # simid_t[j]

    if j == 5:
        stop = 472
    elif j == 6:
        stop = 401
    else:
        stop = 500
    print simid

    b = []
    for tidx in range(250,stop):

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        ds = pa.AthenaDataSet(vtkfname)

        #print(ds.field_list) ; print(ds.derived_field_list)

        d = ds.read_all_data('density')*unit['density'].value # density : solMass/pc3
        vel = ds.read_all_data('velocity')
        velz = vel[:, :, :, 2]
        T1 = ds.read_all_data('T1') ; coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)  # Temperature derived from P/d & cooling func
        #scalar = ds.read_all_data('specific_scalar4')  # ism = 0 / icm = 1
        #surf_d = np.sum(d[scalar < 0.7] * (8 * 8 * 8) / (1024 * 1024.))
        surf_d = np.sum(d* (8 * 8 * 8) / (1024 * 1024.))
        #print ('%s sec' %(time.time()-start_time))

        h0 = 448 ; hh = 510 ; h1 = 573 ; h2 = 698 ; h3 = 823 # 0kpc 0.5kpc 1kpc 2kpc 3pkc
        hhn = 386 ; h1n = 323 ; h2n = 198 ; h3n = 73 # for no icm
        ####### all data at each height ###########
        d0 = d[h0 - 5:h0 + 5]; v0 = velz[h0 - 5:h0 + 5]; t0 = temp[h0 - 5:h0 + 5]; #s0 = scalar[h0 - 5:h0 + 5]

        dh = d[hh - 5:hh + 5]; vh = velz[hh - 5:hh + 5]; th = temp[hh - 5:hh + 5]; #sh = scalar[hh - 5:hh + 5]
        d1 = d[h1 - 5:h1 + 5]; v1 = velz[h1 - 5:h1 + 5]; t1 = temp[h1 - 5:h1 + 5]; #s1 = scalar[h1 - 5:h1 + 5]
        d2 = d[h2 - 5:h2 + 5]; v2 = velz[h2 - 5:h2 + 5]; t2 = temp[h2 - 5:h2 + 5]; #s2 = scalar[h2 - 5:h2 + 5]
        d3 = d[h3 - 5:h3 + 5]; v3 = velz[h3 - 5:h3 + 5]; t3 = temp[h3 - 5:h3 + 5]; #s3 = scalar[h3 - 5:h3 + 5]

        dhn = d[hhn - 5:hhn + 5]; vhn = velz[hhn - 5:hhn + 5]; thn = temp[hhn - 5:hhn + 5]; #sh = scalar[hh - 5:hh + 5]
        d1n = d[h1n - 5:h1n + 5]; v1n = velz[h1n - 5:h1n + 5]; t1n = temp[h1n - 5:h1n + 5]; #s1 = scalar[h1 - 5:h1 + 5]
        d2n = d[h2n - 5:h2n + 5]; v2n = velz[h2n - 5:h2n + 5]; t2n = temp[h2n - 5:h2n + 5]; #s2 = scalar[h2 - 5:h2 + 5]
        d3n = d[h3n - 5:h3n + 5]; v3n = velz[h3n - 5:h3n + 5]; t3n = temp[h3n - 5:h3n + 5]; #s3 = scalar[h3 - 5:h3 + 5]
        '''
        ########## Only ISM at each height ############
        d0 = d0[s0 < 0.7]; v0 = v0[s0 < 0.7]; t0 = t0[s0 < 0.7]
        dh = dh[sh < 0.7]; vh = vh[sh < 0.7]; th = th[sh < 0.7]
        d1 = d1[s1 < 0.7]; v1 = v1[s1 < 0.7]; t1 = t1[s1 < 0.7]
        d2 = d2[s2 < 0.7]; v2 = v2[s2 < 0.7]; t2 = t2[s2 < 0.7]
        d3 = d3[s3 < 0.7]; v3 = v3[s3 < 0.7]; t3 = t3[s3 < 0.7]
        '''
        ########## Devide cold/warm & Hot ###########
        d0_c = d0[t0 < 20000]; v0_c = v0[t0 < 20000]; d0_h = d0[t0 > 20000]; v0_h = v0[t0 > 20000]
        dh_c = dh[th < 20000]; vh_c = vh[th < 20000]; dh_h = dh[th > 20000]; vh_h = vh[th > 20000]
        d1_c = d1[t1 < 20000]; v1_c = v1[t1 < 20000]; d1_h = d1[t1 > 20000]; v1_h = v1[t1 > 20000]
        d2_c = d2[t2 < 20000]; v2_c = v2[t2 < 20000]; d2_h = d2[t2 > 20000]; v2_h = v2[t2 > 20000]
        d3_c = d3[t3 < 20000]; v3_c = v3[t3 < 20000]; d3_h = d3[t3 > 20000]; v3_h = v3[t3 > 20000]

        dh_cn = dhn[thn < 20000]; vh_cn = vhn[thn < 20000]; dh_hn = dhn[thn > 20000]; vh_hn = vhn[thn > 20000]
        d1_cn = d1n[t1n < 20000]; v1_cn = v1n[t1n < 20000]; d1_hn = d1n[t1n > 20000]; v1_hn = v1n[t1n > 20000]
        d2_cn = d2n[t2n < 20000]; v2_cn = v2n[t2n < 20000]; d2_hn = d2n[t2n > 20000]; v2_hn = v2n[t2n > 20000]
        d3_cn = d3n[t3n < 20000]; v3_cn = v3n[t3n < 20000]; d3_hn = d3n[t3n > 20000]; v3_hn = v3n[t3n > 20000]

        if j==0 :
            np.append(dh, dhn); np.append(d1, d1n) ; np.append(d2, d2n) ; np.append(d3, d3n) ;
            np.append(vh, vhn); np.append(v1, v1n) ; np.append(v2, v2n) ; np.append(v3, v3n) ;

            np.append(dh_c,dh_cn); np.append(dh_h,dh_hn)
            np.append(d1_c,d1_cn) ; np.append(d1_h, d1_hn)
            np.append(d2_c, d2_cn);  np.append(d2_h, d2_hn)
            np.append(d3_c, d3_cn);  np.append(d3_h, d3_hn)

            np.append(vh_c, vh_cn);  np.append(vh_h, vh_hn)
            np.append(v1_c, v1_cn ); np.append(v1_h, v1_hn)
            np.append(v2_c, v2_cn);  np.append(v2_h, v2_hn)
            np.append(v3_c, v3_cn);  np.append(v3_h, v3_hn)

        bin = np.arange(-1000, 1000, 5)


        ############ make histogram ############## % density over surf_d is not clear

        hist0, binedges0 = np.histogram(np.ravel(v0), bins=bin,weights=np.ravel(d0/surf_d * (8 * 8 * 8 / (1024 * 1024.))))  # ,density=True)
        x_c0 = (binedges0[1::] + binedges0[0:-1]) / 2
        hist0_c, binedges0_c = np.histogram(np.ravel(v0_c), bins=bin,weights=np.ravel(d0_c/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c0_c = (binedges0_c[1::]+binedges0_c[0:-1])/2
        hist0_h, binedges0_h = np.histogram(np.ravel(v0_h), bins=bin,weights=np.ravel(d0_h/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c0_h = (binedges0_h[1::]+binedges0_h[0:-1])/2

        histh, binedgesh = np.histogram(np.ravel(vh), bins=bin,weights=np.ravel(dh/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_ch = (binedgesh[1::]+binedgesh[0:-1])/2
        histh_c, binedgesh_c = np.histogram(np.ravel(vh_c), bins=bin,weights=np.ravel(dh_c/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_ch_c = (binedgesh_c[1::]+binedgesh_c[0:-1])/2
        histh_h, binedgesh_h = np.histogram(np.ravel(vh_h), bins=bin,weights=np.ravel(dh_h/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_ch_h = (binedgesh_h[1::]+binedgesh_h[0:-1])/2

        hist1, binedges1 = np.histogram(np.ravel(v1), bins=bin,weights=np.ravel(d1/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c1 = (binedges1[1::]+binedges1[0:-1])/2
        hist1_c, binedges1_c = np.histogram(np.ravel(v1_c), bins=bin,weights=np.ravel(d1_c/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c1_c = (binedges1_c[1::]+binedges1_c[0:-1])/2
        hist1_h, binedges1_h = np.histogram(np.ravel(v1_h), bins=bin,weights=np.ravel(d1_h/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c1_h = (binedges1_h[1::]+binedges1_h[0:-1])/2

        hist2, binedges2 = np.histogram(np.ravel(v2), bins=bin,weights=np.ravel(d2/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c2 = (binedges2[1::]+binedges2[0:-1])/2
        hist2_c, binedges2_c = np.histogram(np.ravel(v2_c), bins=bin,weights=np.ravel(d2_c/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c2_c = (binedges2_c[1::]+binedges2_c[0:-1])/2
        hist2_h, binedges2_h = np.histogram(np.ravel(v2_h), bins=bin,weights=np.ravel(d2_h/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c2_h = (binedges2_h[1::]+binedges2_h[0:-1])/2

        hist3, binedges3 = np.histogram(np.ravel(v3), bins=bin,weights=np.ravel(d3/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c3 = (binedges3[1::]+binedges3[0:-1])/2
        hist3_c, binedges3_c = np.histogram(np.ravel(v3_c), bins=bin,weights=np.ravel(d3_c/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c3_c = (binedges3_c[1::]+binedges3_c[0:-1])/2
        hist3_h, binedges3_h = np.histogram(np.ravel(v3_h), bins=bin,weights=np.ravel(d3_h/surf_d*(8*8*8/(1024*1024.))))#,density=True)
        x_c3_h = (binedges3_h[1::]+binedges3_h[0:-1])/2


    ########## data save ##############
        a =np.stack([hist0,hist0_c,hist0_h,histh,histh_c,histh_h,hist1,hist1_c,hist1_h,hist2,hist2_c,hist2_h,hist3,hist3_c,hist3_h])
        b.append(a)
        print tidx

    b = np.array(b)
    #print b.shape
    cPickle.dump(b,open('surf_%s.pkl' % labell[j],'wb'))
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
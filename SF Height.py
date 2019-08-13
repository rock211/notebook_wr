import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
from scipy.interpolate import spline
from matplotlib.ticker import MultipleLocator

sys.path.insert(0, '../')
from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

# In[7]:

import pyathena as pa


unit = pa.set_units(muH=1.4271)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])

# other units can be easily obtained
#print (unit['mass'], unit['time'], unit['magnetic_field'])

Msun = unit['mass'].to('Msun').value
Myr=unit['time'].to('Myr').value

agebin = 10 # unit : Myr, 10 : H-alpha like, 40 : cluster lifetime

M_T_0 = []
M_T_1 = []
M_T_2 = []


simid_t = ('RPS_8pc_noICM_newacc','RPS_8pc_ICM0_newacc','RPS_8pc_ICM1_newacc','RPS_4pc_ICM1_newacc','RPS_8pc_ICM2_newacc','RPS_4pc_ICM2_newacc','RPS_8pc_ICM3_newacc')
#labell = ('No ICM',r'Very Weak','Weak','Weak_4pc','Strong','Strong_4pc',r'Very Strong')
labell = ('P3','P3h', 'P7','P7h','No ICM', 'P14' ,'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
#C = ('k', 'salmon', 'mediumblue','deepskyblue' ,'darkgreen','lime', 'magenta','darkmagenta','goldenrod','royalblue','crimson') # 'plum','orchid','purple'
C = ('dimgray', 'lightskyblue', 'dodgerblue','mediumblue' ,'salmon','crimson', 'firebrick')
C2 = ( 'dodgerblue','mediumblue','salmon','crimson','dimgray', 'firebrick')
#C = ('k', 'salmon', 'deepskyblue', 'green', 'magenta','goldenrod','royalblue','crimson')
M = ((5,1),'o','o','o','o','o','o')
S = ('-','-','-')
L = (1.5,1.5,1.5)
a = (0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8)
a = (1,1,1,1,1,1,1,1,1,1,1)

# overplot Starformation rate of three different simulations
leg = []
norm = 2
fig,ax = plt.subplots(figsize=(10,5))
k=0
model = (1,2,4,3,5,6,0)
for j in model :
    basedir = 'G:/yeongu/'
    simid = simid_t[j]

    if j == 5 or j==6:
        stop = 474
    #elif j == 6:
    #    stop = 401
    else:
        stop = 499
    print simid

    Mj = []
    M_star = 0.0
    o_z = []
    y_z = []
    y_z_t = []
    y_z_b = []
    z = []
    t = []
    m = []
    SFR=[]

    for tidx in range(250, stop):
        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir,simid, simid, tidx)

        # read in domain information
        #ds = pa.AthenaDataSet(vtkfname)

        # full domain information
        #print ds.domain

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

        starfname = vtkfname.replace('id0', 'starpar').replace('vtk', 'starpar.vtk')
        sp = pa.read_starvtk(starfname)
        #print sp

        #mass = sp['mass']*Msun
        #age = sp['age']
        star_clu=sp[sp['mass']!=0]
        new_star = star_clu[star_clu['age']<=1]
        new_z = new_star['x3'].values
        new_t = new_star['time'].values
        #print new_t
        new_m = np.sqrt(new_star['mass'].values * Msun)/norm
        #print list(sp)
        #print np.sum(mass)
        #print len(mass)
        #print len(mass[mass==0]) #runaway star
        #print np.sum(mass[mass!=0]) #cluster
        #print age[mass!=0].min() , age[mass!=0].max()
        #print np.argmin(sp['age']), np.argmax(sp['age'])
        #print sp['age'][52], sp['age'][151]
        #print mass[52],mass[151]
        if len(new_star)!=0 :
            z.append(new_z[0])
            t.append(new_t[0])
            m.append(new_m[0])
            if len(new_z)==2 :
                z.append(new_z[1])
                t.append(new_t[1])
                m.append(new_m[1])
        #print tidx
    #print len(z)
    z = np.array(z)
    m = np.array(m)
    t = np.array(t)
    #print z, z.shape
    #print np.mean(z[t<344]), np.mean(z[t>=344])
    if j==3 or j==5:
        print 896+int(round(np.average(z[t<344],weights=m[t<344])/4)), 896+int(round(np.average(z[t>=344],weights=m[t>=344])/4))
    else:
        print 448+int(round(np.average(z[t < 344], weights=m[t < 344]) / 8)), 448+int(round(np.average(z[t >= 344], weights=m[t >= 344]) / 8))
    t = t * Myr
    #print t


    #print max(m),min(m)

    #sfr=np.append([t,z],[m],axis=0)

    #np.vstack((SFR,z))

    #print sfr.shape

    #np.savetxt('SF_height_%s.txt' % labell[j],sfr,fmt='%s')

    #plt.show()
    #plt.minorticks_on()
    #plt.axes().tick_params(which='minor', direction='in')
    plt.axes().tick_params(which='major', direction='in')
    plt.xlim(250*Myr,500*Myr)
    plt.ylim(-1.1,3.1)
    plt.xticks(fontsize=21)
    plt.yticks([-1,0,1,2,3],fontsize=21)
    plt.xlabel('time [Myr]',fontsize=21)
    plt.ylabel('z [kpc]',fontsize=21)
    plt.scatter(t, z, s=m, c=C[j], edgecolors='k', linewidths=0.65, alpha=a[j], marker=M[j])
    plt.scatter([], [], marker=M[j], s=np.sqrt(3e4)/norm,c=C2[k],label=labell[k],edgecolors='k',linewidths=0.65,alpha=a[k])
    k = k+1
    #plt.legend(loc='upper left',fontsize=14,framealpha=0.1)

starlabels=(r'$10^3 M_\odot$',r'$10^4 M_\odot$',r'$10^5 M_\odot$')
s1 = plt.scatter([], [], c='k', alpha=0.7, s=np.sqrt(1e3)/norm,label=starlabels[0])
s2 = plt.scatter([], [], c='k', alpha=0.7, s=np.sqrt(1e4)/norm,label=starlabels[1])
s3 = plt.scatter([], [], c='k', alpha=0.7, s=np.sqrt(1e5)/norm,label=starlabels[2])
h, l =plt.gca().get_legend_handles_labels()
legend1=plt.legend(h[0:5],l[0:5],loc='upper left',fontsize=14,framealpha=0.3)
plt.legend(h[5:],l[5:],loc='upper right',fontsize=13)
plt.gca().add_artist(legend1)
plt.tight_layout()
#plt.savefig('D:/yeongu/plots/paperplot/new/SFR_height_size_new2.png', dpi=300)
#plt.savefig('D:/yeongu/plots/paperplot/new/SFR_height_size_new2.eps', format='eps',dpi=300,rasterized=True)
plt.show()


'''
        plt.hist(sp['v3'],bins=np.arange(-200,200,5),log=1) #,normed=True
        plt.ylim(0.8,5*1e2)
        plt.xlim(-250,250)
        plt.title('%s_%s' % (labell[j],tidx))
        #plt.savefig('D:/yeongu/plots/sf_v_dis/star_vz_distribution_%s_%s' % (labell[j],tidx))
        #plt.show()

        plt.scatter(sp['x1']/1000.,sp['x3']/1000.,s=5,c=sp['v3'],cmap='jet')
        plt.xlim(-0.512,0.512)
        plt.ylim(-3.584,3.584)
        plt.colorbar()
        plt.clim(-200,200)
        plt.title('%s' % tidx)
        #plt.show()
        plt.close()
    '''

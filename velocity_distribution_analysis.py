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
import seaborn as sea

for tidx in range(250,500):

    data = pd.read_pickle('F:/yeongu/RPS_8pc_ICM1/vel_d/%s/vel_distribution_%s_573' % (tidx,tidx)) # 448, 573
    Temp = data['Temperature']
    temp = Temp[0]
    d = data['density']
    d = d[0]
    vel = data['velocity_z']
    vel = vel[0]
    vel_c = vel[temp < 20000]
    d_c = d[temp < 20000]
    vel_w = vel[temp > 20000]
    d_w = d[temp > 20000]
    #print len(vel_c), len(vel_w), len(vel)
    #print len(vel_w)/len(vel)**2., len(vel_c)/len(vel)**2.
    #vel = vel[Temp>10000]
    #sea.kdeplot(np.ravel(vel),kernel='gau',bw='silverman')
    #plt.show()

    n1, bins1, patches1 = plt.hist(np.ravel(vel),bins=np.arange(vel.min(),vel.max()+5,5),log=True,label='total',facecolor='k',alpha=0.5)#,weights=np.ravel(d))#,normed=1)#,weights=np.ravel(d))
    n3, bins3, patches3 = plt.hist(vel_w,bins=bins1,log=True,label='Hot',facecolor='crimson',alpha=0.7)#,weights=d_w)#,density=len(vel_w)/len(vel)**2.)#,normed=1)#,weights=np.ravel(d))
    n2, bins2, patches2 = plt.hist(vel_c,bins=bins1,log=True,label='Cold/Warm',facecolor='blue',alpha=0.7)#,weights=d_c)#,density=len(vel_c)/len(vel)**2.)#,normed=1)#,weights=np.ravel(d))
    nn = n2.max() + n3.max()
    #print n1.max(), nn
    plt.legend(loc=0)
    plt.xlabel(r'Velocity_z [km/s]')
    plt.title('ICM1_Volume-Weight [z=1kpc] / T = %s ' % tidx)
    #plt.xlim(-1000,1000) #0kpc
    plt.xlim(-250,1000) #1kpc
    plt.ylim(0.5,10**3.8) # Volume
    #plt.ylim(1e-5,200) # Mass
    plt.savefig('D:/yeongu/plots/vel_d/vel_d_Vol_1kpc_%s.png' % tidx ,dpi=100)#,transparent=True)
    #plt.show()
    plt.close()
    print tidx

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

labell = (r'No ICM','ICM00','ICM0','ICM1','ICM2','ICM3','ICM4') # 'NonICM',r'No ICM',
labelll = ('Total','Cold&Warm','Hot')
k = 5
cold = cPickle.load(open('surf_cold_%s.pkl' % labell[k],'rb'))
hot = cPickle.load(open('surf_hot_%s.pkl' % labell[k],'rb'))
label = labell[k]
print cold.shape # tiem / ([hist0,hist0_c,hist0_h,histh,histh_c,histh_h,hist1,hist1_c,hist1_h,hist2,hist2_c,hist2_h,hist3,hist3_c,hist3_h]) # -1000,1000,5

ct = np.mean(cold, axis=0)
c1 = np.mean(cold[0:70], axis=0)
c2 = np.mean(cold[70:170], axis=0)
c3 = np.mean(cold[170:-1], axis=0)
ht = np.mean(hot, axis=0)
h1 = np.mean(hot[0:70], axis=0)
h2 = np.mean(hot[70:170], axis=0)
h3 = np.mean(hot[170:-1], axis=0)

xx = np.arange(-997.5, 997.5, 5)
# print len(xx)

plt.figure(figsize=(15, 7))

C = ['k', 'navy', 'blue', 'deepskyblue', 'skyblue']
alpha = [1, 0.8, 0.7, 1, 1]
labell = ['0 kpc', '0.5 kpc', '1 kpc', '2 kpc', '3 kpc']

j = 0
for i in (448, 510, 573, 698, 823):  # 0,0.5,1,2,3 kpc
    if k == 0:
        coldt = np.sum(ct[i - 3:i + 3] + ct[896 - i - 3:896 - i + 3], axis=0)
        cold1 = np.mean(c1[i - 3:i + 3] + c1[896 - i - 3:896 - i + 3], axis=0)
        cold2 = np.mean(c2[i - 3:i + 3] + c2[896 - i - 3:896 - i + 3], axis=0)
        cold3 = np.mean(c3[i - 3:i + 3] + c3[896 - i - 3:896 - i + 3], axis=0)
        hott = np.sum(ht[i - 3:i + 3] + ht[896 - i - 3:896 - i + 3], axis=0)
        hot1 = np.mean(h1[i - 3:i + 3] + h1[896 - i - 3:896 - i + 3], axis=0)
        hot2 = np.mean(h2[i - 3:i + 3] + h2[896 - i - 3:896 - i + 3], axis=0)
        hot3 = np.mean(h3[i - 3:i + 3] + h3[896 - i - 3:896 - i + 3], axis=0)

    else:
        coldt = np.sum(ct[i - 3:i + 3], axis=0)
        cold1 = np.mean(c1[i - 3:i + 3], axis=0)
        cold2 = np.mean(c2[i - 3:i + 3], axis=0)
        cold3 = np.mean(c3[i - 3:i + 3], axis=0)
        hott = np.sum(ht[i - 3:i + 3], axis=0)
        hot1 = np.mean(h1[i - 3:i + 3], axis=0)
        hot2 = np.mean(h2[i - 3:i + 3], axis=0)
        hot3 = np.mean(h3[i - 3:i + 3], axis=0)

    plt.subplot(2, 4, 1)
    plt.semilogy(xx, coldt, c=C[j], label=labell[j], alpha=alpha[j])
    plt.title('Total : Cold/Warm')
    if k == 0 or k == 1 or k == 2:
        plt.ylim(1e-5, 8 * 1e2)
    else:
        plt.ylim(1e-5, 2 * 1e2)
    plt.legend(loc=0)
    plt.tick_params(direction='in')
    plt.xlim(-250,250)

    plt.subplot(2, 4, 2)
    plt.semilogy(xx, cold1, c=C[j], label=labell[j], alpha=alpha[j])
    plt.title('250~320')
    if k == 0 or k == 1 or k == 2:
        plt.ylim(1e-5, 8 * 1e2)
    else:
        plt.ylim(1e-5, 2 * 1e2)
    plt.tick_params(direction='in')
    plt.xlim(-250,250)

    plt.subplot(2, 4, 3)
    plt.semilogy(xx, cold2, c=C[j], label=labell[j], alpha=alpha[j])
    plt.title('320~420')
    if k == 0 or k == 1 or k == 2:
        plt.ylim(1e-5, 8 * 1e2)
    else:
        plt.ylim(1e-5, 2 * 1e2)
    plt.tick_params(direction='in')
    plt.xlim(-250,250)

    plt.subplot(2, 4, 4)
    plt.semilogy(xx, cold3, c=C[j], label=labell[j], alpha=alpha[j])
    if k != 5:
        plt.title('420~500')
    else:
        plt.title('420-471')

    if k == 0 or k == 1 or k == 2:
        plt.ylim(1e-5, 8 * 1e2)
    else:
        plt.ylim(1e-5, 2 * 1e2)
    plt.tick_params(direction='in')
    plt.xlim(-250,250)

    plt.subplot(2, 4, 5)
    plt.semilogy(xx, hott, c=C[j], label=labell[j], alpha=alpha[j])
    plt.title('Total : Hot')
    plt.ylim(1e-7, 1e1)
    plt.legend(loc=0)
    plt.tick_params(direction='in')

    plt.subplot(2, 4, 6)
    plt.semilogy(xx, hot1, c=C[j], label=labell[j], alpha=alpha[j])
    plt.title('250~320')
    plt.ylim(1e-7, 1e1)
    plt.tick_params(direction='in')

    plt.subplot(2, 4, 7)
    plt.semilogy(xx, hot2, c=C[j], label=labell[j], alpha=alpha[j])
    plt.title('320~420')
    plt.ylim(1e-7, 1e1)
    plt.tick_params(direction='in')

    plt.subplot(2, 4, 8)
    plt.semilogy(xx, hot3, c=C[j], label=labell[j], alpha=alpha[j])
    if k != 5:
        plt.title('420~500')
    else:
        plt.title('420-471')
    plt.ylim(1e-7, 1e1)
    plt.tick_params(direction='in')

    j = j + 1
# plt.ylim(1e-5,2*1e2)
# plt.suptitle('ICM1')
plt.tight_layout()

# label = 'ICM1'
plt.savefig('velocity_%s_modi.png' % label, dpi=500)
#plt.savefig()
plt.show()

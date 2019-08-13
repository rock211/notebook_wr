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


from mpl_toolkits import axes_grid1


from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle
crit_00 = np.log10(0.2036)
crit_0 = np.log10(0.8146)
crit_1 = np.log10(3.2576) # K/cm3
crit_2 = np.log10(6.5152)
crit_3 = np.log10(13.034)
crit_4 = np.log10(32.585)
model = 'ICM00'
#data = pd.read_pickle('F:/yeongu/RPS_8pc_%s/surf_hist/surf_distribution_RPS_8pc_%s_2' %(model,model))
data = pd.read_pickle('F:/yeongu/RPS_8pc_%s/surf_hist/surf_distribution_RPS_8pc_%s_3' %(model,model))
#print data
if model == 'ICM1':
    crit = crit_1
elif model == 'ICM2':
    crit = crit_2
elif model == 'ICM0':
    crit = crit_0
elif model == 'ICM3':
    crit = crit_3
elif model == 'ICM4':
    crit = crit_4
else :
    crit = crit_00

frac_low = []
frac_high = []
stop = 500
for z in range(250,stop):
    data_z = data[z-250:z+1-250]

    under_c = data.loc[str(z):str(z),'-4.995':str(round(crit,2)+0.005)]
    x_under = np.array(map(float,under_c.columns.values))

    over_c = data.loc[str(z):str(z),str(round(crit,2)+0.015):]
    x_over = np.array(map(float,over_c.columns.values))

    N_tot = sum(data_z.values[0])
    N_l = sum(under_c.values[0])
    N_h = sum(over_c.values[0])

    frac_l = N_l/float(N_tot)
    frac_h = N_h/float(N_tot)



    ############## surface density histogram #####################
    plt.bar(x_under,under_c.values[0],width=0.01,facecolor='b',alpha=0.7,label='Low')
    plt.bar(x_over,over_c.values[0],width=0.01,facecolor='r',alpha=0.7,label='High')
    if model == 'ICM1':
        plt.axvline(round(crit,2),ls='--',c='k',label=r'$\Sigma_{crit}$') # ICM1
        plt.text(-2.8,400*0.95,r'Critical density = %s (~3.25M$_{\odot} pc^{-2}$)' % round(crit,3)) #data_z.values[0].max() / ICM1
    elif model =='ICM2' :
        plt.axvline(round(crit, 3), ls='--', c='k', label=r'$\Sigma_{crit}$')  # ICM2
        plt.text(-2.8, 400 * 0.95, r'Critical density = %s (~6.51M$_{\odot} pc^{-2}$)' % round(crit, 3)) # ICM2
    elif model == 'ICM0' :
        plt.axvline(round(crit, 3), ls='--', c='k', label=r'$\Sigma_{crit}$')  # ICM0
        plt.text(-2.8, 400 * 0.95, r'Critical density = %s (~0.814M$_{\odot} pc^{-2}$)' % round(crit, 3))  # ICM0
    elif model == 'ICM3':
        plt.axvline(round(crit, 3), ls='--', c='k', label=r'$\Sigma_{crit}$')  # ICM3
        plt.text(-2.8, 400 * 0.95, r'Critical density = %s (~13M$_{\odot} pc^{-2}$)' % round(crit, 3))  # ICM3
    elif model == 'ICM4':
        plt.axvline(round(crit, 3), ls='--', c='k', label=r'$\Sigma_{crit}$')  # ICM4
        plt.text(-2.8, 400 * 0.95, r'Critical density = %s (~32.58M$_{\odot} pc^{-2}$)' % round(crit, 3))  # ICM4
    else:
        plt.axvline(round(crit, 3), ls='--', c='k', label=r'$\Sigma_{crit}$')  # ICM00
        plt.text(-2.8, 400 * 0.95, r'Critical density = %s (~0.203M$_{\odot} pc^{-2}$)' % round(crit, 3))  # ICM00

    plt.text(-2.8,400*0.9,'Fraction_low = %s' % round(frac_l,3))
    plt.xlim(-3,3)
    plt.ylim(0,400)
    plt.xlabel(r'$log$ surface density [M$_{\odot} pc^{-2}$]')
    plt.ylabel('Number')
    plt.title('%s_%s' % (model,z))
    plt.legend(loc=0)
    plt.savefig('D:/yeongu/plots/surfd_hist_new/ICM00/Hist_%s_%s.png' % (model,z), dpi=100)
    #plt.show()
    plt.close()
    print z
    ###############################################################


######################## fraction history #################
    frac_high.append(frac_h)
    frac_low.append(frac_l)
    print z
xx = range(250,stop)

plt.plot(xx,frac_high,c='r',alpha=0.7,label='High')
plt.plot(xx,frac_low,c='b',alpha=0.7,label='Low')

plt.xlabel('time')
plt.ylabel('fraction')
plt.ylim(0,1)
plt.legend(loc=1)
plt.title('%s' % model)

plt.show()
############################################################

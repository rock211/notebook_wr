import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy.constants as c
import astropy.units as u
import sys
import pyathena as pa
from mpl_toolkits import axes_grid1
import copy
from matplotlib.colors import LogNorm
from six.moves import cPickle as pickle

unit = pa.set_units(muH=1.4271)
# print(unit)
# print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
kb = 1.3806504 * 1e-16  # boltzmann constant erg/K

simid_t = (
'RPS_8pc_noICM_newacc', 'RPS_8pc_ICM0_newacc', 'RPS_8pc_ICM1_newacc', 'RPS_4pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc',
'RPS_4pc_ICM2_newacc', 'RPS_8pc_ICM3_newacc')  # 'MHD_8pc_new' ,
# labell = ('No ICM','Very Weak', 'Weak', 'Strong','Very Strong' ,'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
labell = ('No ICM', 'P1', 'P3', 'P3h', 'P7', 'P7h', 'P14', 'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
C = ('dimgrey', 'coral', 'royalblue')
S = ('-.', '--', '-')

# overplot Starformat
# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])
s_cut = 0.9
for j in (0, 10):

    basedir = 'G:/yeongu/'

    simid = simid_t[j]
    Mj = []
    M_star_c = 0.0
    Mass_tot = []
    Mc = []
    SFE = []

    if j == 6 or j == 5:
        stop = 474
    else:
        stop = 499

    for tidx in (250,275,300,325,350,375,400,425,450):#range(250, stop):  # time step

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # name of original data fields we stored from the simulation
        # print(ds.field_list)

        # It also has predefined data fields can be calculated from the original data.
        # print(ds.derived_field_list)

        ############################## general properties ####################################
        #T1 = ds.read_all_data('T1');
        #coolftn = pa.coolftn()
        #temp = coolftn.get_temp(T1)
        #tempc = copy.copy(temp)

        if j != 0:
            scalar = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1

        #######################################################################################

        #comp = ds.read_all_data('specific_scalar3')  # ism = 0 / icm = 1
        #compc = copy.copy(comp)
        # d = ds.read_all_data('density')*unit['density'].value # density
        # dc = copy.copy(d)

        # pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
        # prec = copy.copy(pre)
        nd = ds.read_all_data('number_density')  # number density
        ndc = copy.copy(nd)
        # scalar cut data
        # d = (d[comp < cut])
        # P = (pre[comp < cut]/kb) # pressure(P/kb)
        # nd = (nd[comp < cut]) # number density
        # T = (tem[comp < cut]) # temperature
        # MagF = np.log10(magF[comp < cut]) # magnetic field
        # MagP = np.log10(magP[comp < cut]) # magnetic pressure

        ######### ISM ##########
        # P = pre[scalar < s_cut]
        # d = d[scalar < s_cut]
        if j!=0:
            nd = nd[scalar < s_cut]
            ndc = ndc[scalar < s_cut]
        #temp = temp[scalar < s_cut]

        # P = np.ravel(P)
        # d = np.ravel(d)
        nd = np.log10(np.ravel(nd))
        ndc = np.ravel(ndc)
        #T = np.ravel(temp)
        print nd.min(),nd.max()
        hist, edge = np.histogram(nd,bins=np.arange(-6,2.1,0.05),weights=ndc)
        hist = hist/float(np.sum(hist))

        plt.plot(edge[0:-1],hist,label='%s' % tidx)
        #plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-7,1)
        plt.xlim(-5,2)
        plt.xlabel('Log[n ($cm^{-3}$)]')
        plt.ylabel('Probability')
        plt.title('%s_%s' % (simid,tidx))
        #plt.savefig('D:/yeongu/plots/1d_den/%s_%s.png' % (simid, tidx))

    #plt.close()
    plt.legend()

    plt.show()

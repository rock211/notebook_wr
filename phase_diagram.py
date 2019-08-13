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



# ## Unit system
#
# The unit system we choose for this simulation is
# * [length] = pc
# * [velocity] = km/s
# * [density] = 1.4271*m_h/cm^3

unit = pa.set_units(muH=1.4271)
#print(unit)
#print (unit['density'].cgs / 1.4271 / c.m_p.cgs, unit['velocity'], unit['length'])
kb = 1.3806504 * 1e-16 #boltzmann constant erg/K

simid_t = ('RPS_8pc_noICM_newacc', 'RPS_8pc_ICM0_newacc', 'RPS_8pc_ICM1_newacc','RPS_4pc_ICM1_newacc', 'RPS_8pc_ICM2_newacc','RPS_4pc_ICM2_newacc', 'RPS_8pc_ICM3_newacc')  # 'MHD_8pc_new' ,
#labell = ('No ICM','Very Weak', 'Weak', 'Strong','Very Strong' ,'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
labell = ('No ICM','P1', 'P3','P3h', 'P7','P7h','P14' ,'ICM1', 'ICM2', 'ICM3', 'ICM4')  # r'No ICM',
C = ('dimgrey','coral','royalblue')
S = ('-.','--','-')

# overplot Starformat
# other units can be easily obtained
# print (unit['mass'], unit['time'], unit['magnetic_field'], unit['temperature'])
for j in (2,4):

    basedir = 'G:/yeongu/'

    simid = simid_t[j]
    Mj = []
    M_star_c = 0.0
    Mass_tot = []
    Mc = []
    SFE = []

    if j == 6 or j==5:
        stop = 474
    else:
        stop = 499

    for tidx in range(251, stop):  # time step

        vtkfname = '%s%s/id0/%s.%04d.vtk' % (basedir, simid, simid, tidx)
        # read in domain information
        ds = pa.AthenaDataSet(vtkfname)

        # name of original data fields we stored from the simulation
        #print(ds.field_list)

        # It also has predefined data fields can be calculated from the original data.
        #print(ds.derived_field_list)

        ############################## general properties ####################################
        T1 = ds.read_all_data('T1'); coolftn = pa.coolftn()
        temp = coolftn.get_temp(T1)
        tempc = copy.copy(temp)

        if j != 0:
            scalar = ds.read_all_data('specific_scalar3') # ism = 0 / icm = 1
            s_cut = 0.5
        #######################################################################################


        #comp = ds.read_all_data('specific_scalar3') # ism = 0 / icm = 1
        #compc = copy.copy(comp)
        #d = ds.read_all_data('density')*unit['density'].value # density
        #dc = copy.copy(d)

        #pre = ds.read_all_data('pressure')*unit['pressure'].value/kb # thermal pressure
        #prec = copy.copy(pre)
        nd = ds.read_all_data('number_density') # number density
        #ndc = copy.copy(nd)
        # scalar cut data
        #d = (d[comp < cut])
        #P = (pre[comp < cut]/kb) # pressure(P/kb)
        #nd = (nd[comp < cut]) # number density
        #T = (tem[comp < cut]) # temperature
        #MagF = np.log10(magF[comp < cut]) # magnetic field
        #MagP = np.log10(magP[comp < cut]) # magnetic pressure
        '''
        ######### ISM ##########
        P = pre[scalar < s_cut]
        d = d[scalar < s_cut]
        nd = nd[scalar < s_cut]
        temp = temp[scalar < s_cut]

        P = np.ravel(P)
        d = np.ravel(d)
        nd = np.ravel(nd)
        T = np.ravel(temp)
        
        ######### ICM #########
        Pc = prec[scalar > s_cut]
        dc = dc[scalar > s_cut]
        ndc = ndc[scalar > s_cut]
        tempc = tempc[scalar > s_cut]
        scalar3 = scalar[scalar > s_cut]

        Pc = np.ravel(Pc)
        dc = np.ravel(dc)
        ndc = np.ravel(ndc)
        Tc = np.ravel(tempc)
        sc3 = np.ravel(scalar3)
        '''


        ############# Mixing gas #################
        nd = ds.read_all_data('number_density')
        #print nd.max(), nd.min()
        #plt.hist(np.ravel(scalar),bins=np.arange(0.5,1.2,0.005),log=True)#,weights=np.ravel(nd))
        #plt.show()
        velz = ds.read_all_data('velocity')[:,:,:,2]

        scalar_range = (scalar <= 0.99) & (scalar >=0.01) # mixing gas 0.99~0.01

        scalar_mix = scalar[scalar_range]

        nd_mix = nd[scalar_range] ; nd_mix = np.ravel(nd_mix)
        velz_mix = velz[scalar_range] ; velz_mix = np.ravel(velz_mix)
        temp_mix = temp[scalar_range] ; temp_mix = np.ravel(temp_mix)
        #print nd_mix.max(), nd_mix.min()
        #nd_mix = np.ravel(nd) ; velz_mix = np.ravel(velz) ; temp_mix = np.ravel(temp)
        ######################################################################################
        #print velz_mix.max(),velz_mix.min()
        ############## Mixing figure (T-Vz, D-Vz) ###############
        plt.figure(figsize=(6,8))
        xl = -5.5 ; xh = 3 ; xbin = xh - xl
        yl = -1 ; yh = 4 ; ybin = yh - yl
        xbins = np.linspace(xl,xh,xbin*30+1)
        ybins = np.linspace(yl,yh,ybin*30+1)

        #bg_color = 'black'; fg_color = 'white'
        #fig = plt.figure(facecolor=bg_color, edgecolor=fg_color)
        #axes = plt.axes(facecolor=bg_color)
        #axes.xaxis.set_tick_params(color=fg_color,direction='in'); axes.yaxis.set_tick_params(color=fg_color,direction='in')

        H, xe, ye = np.histogram2d(np.log10(nd_mix), np.log10(velz_mix), weights=nd_mix*scalar_mix, bins=(xbins, ybins))
        H = H/np.sum(H) # Normalize
        H = H.T # Transpose

        plt.subplot(2,1,1)
        plt.imshow(H,origin='low', extent=[xe[0], xe[-1], ye[0], ye[-1]],norm=LogNorm(),cmap='cool_r',aspect='auto')
        plt.colorbar()
        plt.clim(1e-6, 0.01)
        plt.title(r'V$_z$-D / %s_%s' % (labell[j],tidx))
        plt.xlabel(r'N density (log n $[cm^{-3}])$')
        plt.ylabel(r'log V$_z$ [km/s]')
        plt.ylim(-1, 4)
        plt.tight_layout()


        xl = 0 ; xh = 8.5 ; xbin = xh - xl
        yl = -1 ; yh = 4 ; ybin = yh - yl
        xbins = np.linspace(xl,xh,xbin*30+1)
        ybins = np.linspace(yl,yh,ybin*30+1)

        H2, xe2, ye2 = np.histogram2d(np.log10(temp_mix), np.log10(velz_mix), weights=nd_mix*scalar_mix, bins=(xbins, ybins))
        H2 = H2/np.sum(H2) # Normalize
        H2 = H2.T # Transpose

        plt.subplot(2,1,2)
        plt.imshow(H2, origin='low', extent=[xe2[0], xe2[-1], ye2[0], ye2[-1]], norm=LogNorm(), cmap='cool_r',aspect='auto')
        plt.colorbar()
        plt.clim(1e-6, 0.01)
        plt.title(r'V$_z$-T / %s_%s' % (labell[j],tidx))
        plt.xlabel(r'log T [K]')
        plt.ylabel(r'log V$_z$ [km/s]')
        plt.ylim(-1,4)
        #Normalize problem!
        #plt.xlim(xl,xh)
        #plt.ylim(yl,yh)
        #cbar = plt.colorbar()
        #plt.clim(1e-6, 0.01)
        plt.tight_layout()
        #plt.show()
        print tidx
        plt.savefig('D:/yeongu/plots/mix/mixing_%s_%s_0.01-0.99.png' % (simid,tidx))
        plt.close()


        '''
        ############## P-D figure ###############

        xl = -6 ; xh = 3 ; xbin = xh - xl
        yl = -3 ; yh = 6.5 ; ybin = yh - yl
        xbins = np.linspace(xl,xh,xbin*50+1)
        ybins = np.linspace(yl,yh,ybin*50+1)

        bg_color = 'black'; fg_color = 'white'
        hst = pa.hst_reader('%s/hst/%s.hst' % (ds.dir, ds.id))
        hratio = np.interp(ds.domain['time'], hst.time, hst.heat_ratio)
        coolftn = pa.coolftn()
        nden = coolftn.heat * hratio / coolftn.cool
        peq = coolftn.temp * nden

        fig = plt.figure(facecolor=bg_color, edgecolor=fg_color)
        axes = plt.axes(facecolor=bg_color)
        axes.xaxis.set_tick_params(color=fg_color,direction='in'); axes.yaxis.set_tick_params(color=fg_color,direction='in')

        H, xe, ye = np.histogram2d(np.log10(nd), np.log10(P), weights=d, bins=(xbins, ybins))
        H = H/np.sum(H) # Normalize
        H = H.T # Transpose

        H2, xe2, ye2 = np.histogram2d(np.log10(ndc), np.log10(Pc), weights=dc*sc3, bins=(xbins, ybins))
        H2 = H2/np.sum(H2) # Normalize
        H2 = H2.T # Transpose

        #plt.plot(np.log10(nden), np.log10(peq),'w--',lw=1)
        #plt.imshow(H,origin='low', extent=[xe[0], xe[-1], ye[0], ye[-1]],norm=LogNorm(),cmap='RdYlBu')
        #cbar = plt.colorbar()
        #plt.clim(1e-6, 0.01)

        plt.imshow(H2, origin='low', extent=[xe[0], xe[-1], ye[0], ye[-1]], norm=LogNorm(), cmap='BuPu')
        #Normalize problem!
        plt.xlim(xl,xh)
        plt.ylim(yl,yh)
        cbar = plt.colorbar()
        plt.clim(1e-6, 0.01)
        plt.title(r'P-D_%s_%s' % (labell[j],tidx))
        #plt.text(0,-1.5,r'scalar cut = %s' % cut) #,color='white'
        #plt.text(0, -1.5, '%s' % labell[j]) # for nonICM
        plt.text(0,-2,r'Mass Weighted') #,color='white'
        plt.xlabel(r'N density (log n $[cm^{-3}])$')
        plt.ylabel(r'Pressure (log $P/k_b$ $[Kcm^{-3}]$)')
        plt.tight_layout()
        #plt.show()
        print tidx
        plt.savefig('D:/yeongu/plots/pd_new/scalarweight_ICM_P-D_only_%s_%s.png' % (labell[j],tidx))
        plt.close()



        ############## T-D figure ###############
        bg_color = 'black'
        fg_color = 'white'

        fig = plt.figure(facecolor=bg_color, edgecolor=fg_color)
        axes = plt.axes(facecolor=bg_color)
        axes.xaxis.set_tick_params(color=fg_color,direction='in')
        axes.yaxis.set_tick_params(color=fg_color,direction='in')

        xl = -6 ; xh = 3 ; xbin = xh - xl
        yl = 0 ; yh = 8.5 ; ybin = yh - yl
        xbins = np.linspace(xl,xh,xbin*50+1)
        ybins = np.linspace(yl,yh,ybin*50+1)
        #print np.array(nd), tem.shape, d.shape
        #H, xe, ye, img = plt.hist2d(np.log10(nd),np.log10(T),weights=d,norm=LogNorm(),bins=(xbins,ybins),cmap='BuPu')


        H, xe, ye = np.histogram2d(np.log10(nd),np.log10(T),weights=d,bins=(xbins,ybins))
        H = H/np.sum(H) # Normalize
        H = H.T # Transpose
        #print H.max(), np.min(np.nonzero(H))
        #print np.max(H)
        #plt.imshow(H,origin='low', extent=[xe[0], xe[-1], ye[0], ye[-1]],norm=LogNorm(),cmap='RdYlBu')
        #plt.xlim(xl,xh)
        #plt.ylim(2,yh)
        #cbar = plt.colorbar()
        #plt.clim(1e-6,0.01)

        H2, xe2, ye2 = np.histogram2d(np.log10(ndc), np.log10(Tc), weights=dc*sc3, bins=(xbins, ybins))
        H2 = H2/np.sum(H2) # Normalize
        H2 = H2.T # Transpose

        plt.imshow(H2, origin='low', extent=[xe2[0], xe2[-1], ye2[0], ye2[-1]], norm=LogNorm(), cmap='BuPu')
        cbar = plt.colorbar()
        plt.clim(1e-6,0.01)
        plt.title(r'T-D phase diagram_%s_%s' % (labell[j],tidx))
        #plt.text(0,7,r'scalar cut = %s' % cut)
        #plt.text(0, 7, 'nonICM') # for non ICM
        #plt.text(0,7.5, r'Mass Weighted')
        #plt.tight_layout()
        plt.xlabel(r'log n $[cm^{-3}]$')
        plt.ylabel(r'log T [K]')
        plt.tight_layout()
        #plt.show()

        plt.savefig('D:/yeongu/plots/pd_new/scalarweight_ICM_T-D_only_%s_%s.png' % (labell[j],tidx))
        plt.close()
        print tidx
        '''
'''
    ############## MagP-D figure ###############
        bg_color = 'black'
        fg_color = 'white'

        fig = plt.figure(facecolor=bg_color, edgecolor=fg_color)
        axes = plt.axes(facecolor=bg_color)
        xl = -6 ; xh = 3 ; xbin = xh - xl
        yl = -7 ; yh = 2 ; ybin = yh - yl

        xbins = np.linspace(xl,xh,xbin*100+1)
        ybins = np.linspace(yl,yh,ybin*100+1)
        #print np.array(nd), tem.shape, d.shape
        #H, xe, ye, img = plt.hist2d(np.log10(nd),np.log10(T),weights=d,norm=LogNorm(),bins=(xbins,ybins),cmap='BuPu')
        H, xe, ye = np.histogram2d(np.log10(nd),np.log10(magF),weights=d,bins=(xbins,ybins))
        H = H/np.sum(H) # Normalize
        H = H.T # Transpose

        #print np.max(H)
        plt.imshow(H,origin='low', extent=[xe[0], xe[-1], ye[0], ye[-1]],norm=LogNorm(),cmap='RdYlBu')
        cbar = plt.colorbar()
        plt.clim(0.0000000001,0.005)
        
        plt.imshow(H2, origin='low', extent=[xe[0], xe[-1], ye[0], ye[-1]], norm=LogNorm(), cmap='BuPu')
        #plt.xlim(xl,xh)
        #plt.ylim(yl,yh)
        cbar = plt.colorbar()
        plt.clim(0.0000000001,0.005)
        plt.title(r'B-D phase diagram_%s_%s' % (labell[j],tidx))
        #plt.text(0,7,r'scalar cut = %s' % cut)
        #plt.text(0, 7, 'nonICM') # for non ICM
        #plt.text(0,7.5, r'Mass Weighted')
        #plt.tight_layout()
        plt.xlabel(r'log n $[cm^{-3}]$')
        plt.ylabel(r'B $[\mu Gauss]$')
        #plt.show()

        plt.savefig('D:/yeongu/plots/bd_new/B-D_MassW_%s_%s.png' % (labell[j],tidx))
        plt.close()
        print tidx
'''

    






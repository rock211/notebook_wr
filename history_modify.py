
# # History variables
#
# * Volume average := $$[q]_p \equiv\sum q \Theta(p)\Delta V / \sum \Delta V$$
# * Area average := $$\langle q_k \rangle_p \equiv \frac{\sum_{i,j}q_k\Theta(p)\Delta x \Delta y}{L_x L_y}$$
# * Default variables:
#     * time, dt
#     * mass = $[\rho]$
#     * totalE = $[E]$
#     * x(1/2/3)Mom = $[\rho v_{(1/2/3)}]$
#     * x(1/2/3)KE = $[E_{k,(1/2/3)}=\rho v_{(1/2/3)}^2/2]$
#     * x(1/2/3)ME = $[B_{1/2/3}^2/2]$
#     * gravPE = $[\rho\Phi]$
# * Additional variables:
#     * heat_ratio = $[\Gamma/\Gamma_0]$
#     * x2dke = $[\delta E_{k,2}=\rho (v_2+q\Omega x)^2/2]$
#     * x(1/2/3)KE_2p = $[\delta E_{k,(1/2/3)}]_{2p}$, where $2p=c+u+w$
#     * F3(h2/h1/w/u/c) =
#     $0.5 ( \langle (\rho v_3)_{k=ke+1}\rangle_{(h2/h1/w/u/c)} - \langle (\rho v_3)_{k=ks}\rangle_{(h2/h1/w/u/c)})$
#     * H2 = $[\rho z^2]$
#     * H2(h2/h1/w/u/c) = $[\rho z^2](h2/h1/w/u/c)$
#     * P = $[P]$
#     * Pth(_2p) = $ \langle P_{k = kmid-1} + P_{k =kmid}\rangle_{(2p)} /2$, where $kmid= ks+ Nz/2 -1$
#     * Pturb(_2p) = $\langle (\rho v_3^2)_{k = kmid-1} + (\rho v_3^2)_{k =kmid}\rangle_{(2p)}/2$
#     * nmid(_2p) = $\langle \rho_{k = kmid-1} + \rho_{k =kmid}\rangle_{(2p)}/2$
#     * Vmid_2p = $\langle \Theta(2p) \rangle/2$
#     * V(h2/h1/w/u/c) = $[\Theta(h2/h1/w/u/c)]$
#     * M(h2/h1/w/u/c) = $[\rho\Theta(h2/h1/w/u/c)]$
#     * B(1/2/3) = $[B_{i,(1/2/3)}]$
#     * sfr(10/40/100) = $\Sigma_{\rm SFR} (\Delta t =  (10/40/100){\rm Myr})$
#     * msp = $[M_{\rm sp}/\Delta V]$ if a grid zone (i,j,k) has star particles
#     * mghost = $[\rho]$ if a grid zone (i,j,k) is in control volume of star particles

# # packages to be used
#  * matplotlib http://matplotlib.org/
#  * numpy
#  * pandas http://pandas.pydata.org/
#  * astropy http://www.astropy.org/

# In[ ]:

#get_ipython().magic(u'matplotlib inline')

import matplotlib
matplotlib.__version__

import numpy
numpy.__version__

import pandas
pandas.__version__

import astropy
astropy.__version__

import sys
sys.path.insert(0,'../')

from pyathena import ath_hst

ath_hst.__file__


# # ath_hst.py
#
# * Athena history dump can be found under the directory "id0/" with an extension ".hst"
# * To read the history dump, "ath_hst.py" can be used.
# * There are two functions you can use:
#     * hst = ath_hst.read(hstfilename)
#         * uses standard ascii io to read in the history dump
#         * returns a distionary with keys corresponding to the history variable
#         * using each key, you can retrive sequence of data as a numpy array
#     * hstp = ath_hst.read_w_pandas(hstfilename,write=True)
#         * uses pandas packages to read in the history dump
#         * returns pandas DataFrame contains all the information
#         * if write=True, it automatically write the DataFrame to "pickle"
#         * if there is a pickle file with filename = hstfilename + '.p', and the pickle file is newer than the original history file, it automatically read the "pickle" file, which is much faster than original ascii file
#         * SN history dump can also be accessed with this function


from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import astropy

# Call data
simid_t = ('R8_8pc_metal' , 'RPS_8pc_ICM1' , 'RPS_8pc_ICM2')
labell = (r'No ICM','n1','n2')
C = ('dimgrey','coral','royalblue')

for j in range(3) :
    basedir='F:/yeongu/'
    id=simid_t[j] # 'MHD_8pc_new' , 'RPS_8pc_n1e-4_v1414' , 'RPS_8pc_n2e-4_v1414'
    hstfilename= basedir+id+'/hst/'+id+'.hst'

    hstp=ath_hst.read_w_pandas(hstfilename)

#print hstp.mass
    hstp = hstp[(hstp.time > 250)]
    hstp = hstp[(hstp.time < 500)]
    #print hstp
#print(hstp.mass)
#print(hstp.B1)

#print(np.where(hstp.mass>0.042))
    #hstp.plot(x='time',y='mass')
    #hstp.plot(x='time', y=['Mc', 'Mu', 'Mw', 'Mh1', 'Mh2'])
    '''
    plt.plot(hstp.time,hstp.mass,label='Total')
    plt.plot(hstp.time,hstp.Mc,label='Cold')
    plt.plot(hstp.time, hstp.Mu,label='Unstable')
    plt.plot(hstp.time, hstp.Mw,label='Warm')
    plt.plot(hstp.time, hstp.Mh1,label='H1')
    plt.plot(hstp.time, hstp.Mh2,label='H2')
    plt.title('Mass_Frac_%s' % labell[j])
    plt.legend(loc=1)
    # plt.savefig('Mass_Frac.png')
    #plt.show()
    #plt.close()
    '''
    # Mass history

    plt.plot(hstp.time,hstp.mass,label='%s_Gas' % labell[j],color=C[j],alpha=0.8)
    plt.plot(hstp.time,hstp.msp,label='%s_star' % labell[j],color=C[j],linestyle='--')

plt.xlabel(r'Time [Myr]')
plt.ylabel('mass')
plt.title('Total mass variation')
plt.legend(loc=0)
#plt.savefig('mass_time.png')
plt.show()

# SF history
'''
hstp.plot(x='time',y=['sfr10'])
plt.title('SF10')
plt.savefig('SF10.png')
plt.show()

hstp.plot(x='time',y=['sfr40'])
plt.title('SF40')
plt.savefig('SF40.png')
plt.show()

hstp.plot(x='time',y=['sfr100'])
plt.title('SF100')
#plt.savefig('SF100.png')
plt.show()

# Mass Fractions from the history dump

hstp.plot(x='time',y=['Mc','Mu','Mw','Mh1','Mh2'])
#plt.savefig('Mass_Frac.png')
plt.show()

# to convert the mean density of each phase to fraction
phase=['c','u','w','h1','h2']
for p in phase:
    hstp['f'+p]=hstp['M'+p]/hstp['mass']
    hstp.plot(x='time',y='f'+p,ax=plt.gca())
plt.savefig('Mean_density_to_frac.png')
plt.title('Mass Fraction')
plt.show()

# # Scale heights # Do not under stand yet

hstp['H']=np.sqrt(hstp['H2']/hstp['mass'])
hstp.plot(x='time',y='H')
#plt.savefig('H2_scale_height.png')
plt.show()

# In[20]:

phase=['c','u','w','h1','h2']
for p in phase:
    hstp['H'+p]=np.sqrt(hstp['H2'+p]/hstp['M'+p])
    hstp.plot(x='time',y='H'+p,ax=plt.gca())
plt.yscale('log')
plt.show()
'''
# # SN data
# * id = id of the host star particle
# * time = time of explosion
# * age = age of the host star particle; for runaway, this is just indicator of clock. SN exploded as age exceeds zero
# * mage = mass-weighted age of the host star particle; for runaway, this is time since runaway creation
# * mass = mass of the host star particle. can be used to distinguish runaways
# * (x1/x2/x3) = position of the host star particle
# * (x1sn/x2sn/x3sn) = position of SN explosion; now it is identical to the star particle position, but possibly we can add a distribution from the center
# * (n/v1/v2/v3/e)avg = mean gas properties within SNR
# * vol = volume within SNR
# * radius = radius of SNR
# * SFUV = $\Sigma_{\rm FUV}$ for the host star particle
# * SNRate = SN rate calculated based on the host star particle properties #what is unit?
# * SNprob = randomly generated probability; this should be smaller than SNRate for explosion
# * runaway = 1 if SN from runaways
# * parent = parent star cluster id of the runaway
# * mode = feedback type
# * active = active flag of the host star particle
# * fm = mass within the SNR/M_sf

sn=ath_hst.read_w_pandas(hstfilename.replace('.hst','.sn'))
print sn['time']
#print sn['SNRate'][]
plt.scatter(sn['time'][(sn['SNRate'] < 3000) & (sn['SNRate'] > 0)],sn['SNRate'][(sn['SNRate'] < 3000) & (sn['SNRate'] > 0)],s=1)
plt.show()
# Position history

plt.plot(sn['time'],sn['x3'],'.')
#plt.savefig('posi_x3.png')
plt.show()

# 3D distribution of Star

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(sn['x1'], sn['x2'], sn['x3'])

plt.title('posi_all')
#plt.savefig('posi_all.png')
plt.show()

# In[29]:

runaway=sn['mass'] == 0.0
rsn=sn[runaway]
csn=sn[~runaway]


# In[30]:

plt.plot(rsn['time'],rsn['x3'],'.')
plt.plot(csn['time'],csn['x3'],'.')
#plt.savefig('7.png')
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM

sn = np.genfromtxt('SCPUnion2.1_mu_vs_z.txt')

#print sn

tt = np.genfromtxt('chi_100_re.txt')


z = sn[:,1]
DM = sn[:,2]
DM_err = sn[:,3]

ho = 70
zz = np.arange(z.min(),1.5,0.001)
print zz
DMm = []
for k in range(len(zz)):
    d_L = LambdaCDM(H0=ho, Om0=0.278, Ode0=0.723)
    d_L = d_L.luminosity_distance(zz[k]).value
    # print k, z[k],d_L
    DM_m = 5 * np.log10(d_L) + 25
    DMm.append(DM_m)
    #print k
plt.plot(zz,DMm,c='r',linewidth=2)
plt.scatter(z,DM,facecolor='white',edgecolors='b',linewidths=1)
plt.errorbar(z,DM,yerr=DM_err,color='b',alpha=0.5,capsize=2,fmt='None')
plt.text(1.2,42,r'$\Omega_m$=0.278',size=14)
plt.text(1.2,41,r'$\Omega_{\Lambda}$=0.723',size=14)
plt.xlabel(r'redshift [z]')
plt.ylabel(r'Distance Modulus')
plt.title(r'SN data')
plt.show()

size = 100
chi = np.zeros((8,8))
#print chi
for i in np.arange(0.26,0.34,1./size): # Omega Matter
    for j in np.arange(0.26,0.34,1./size): # Omega DarkEnergy
        #chisq = 0
        chis = []
        for k in range(len(z)) :
            #print chisq
            d_L = LambdaCDM(H0=ho,Om0=i,Ode0=j)
            d_L = d_L.luminosity_distance(z[k]).value
            #print k, z[k],d_L
            DM_m = 5*np.log10(d_L)+25

            chisq = ((DM[k]-DM_m)/DM_err[k])**2
            #plt.scatter(z[k],DM_m)
            chis.append(chisq)
            #print i,j,chis
        #plt.scatter(z,DM)
        #plt.show()
        #print chis
        print i,j, np.sum(chis)
        x=int((i-0.26)*size)
        y=int((j-0.26)*size)
        print x,y
        chi[x,y]=np.sum(chis)


print chi, len(chi)


plt.imshow(chi,origin='lower')

plt.show()
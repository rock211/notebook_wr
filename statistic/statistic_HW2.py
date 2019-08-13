import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM

sn = np.genfromtxt('SCPUnion2.1_mu_vs_z.txt')

#print sn


z = sn[:,1]
DM = sn[:,2]
DM_err = sn[:,3]

plt.scatter(z,DM,facecolor='white',edgecolors='b',linewidths=1)
plt.errorbar(z,DM,yerr=DM_err,color='b',alpha=0.5,capsize=2,fmt='None')
plt.xlabel(r'redshift [z]')
plt.ylabel(r'Distance Modulus')
plt.title(r'SN data')
plt.show()
#z = z[DM_err > 0]
#DM = DM[DM_err > 0]
#DM_err = DM_err[DM_err > 0]
#print len(z)
ho = 70

size = 100
chi = np.zeros((size,size))
#print chi
for i in np.arange(0,1,1./size): # Omega Matter
    for j in np.arange(0,1,1./size): # Omega DarkEnergy
        #chisq = 0
        chis = []
        for k in range(len(z)) : # Put this part at initial of while loop
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
        chi[int((i)*size),int((j)*size)]=np.sum(chis)
    np.savetxt('chi_100_re_%s.txt' % i, chi)
    print chi, np.amin(chi)

print chi, len(chi)

plt.imshow(chi,origin='lower')
np.savetxt('chi_100_re.txt',chi)
plt.show()
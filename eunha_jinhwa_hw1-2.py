import numpy as np
#import astropy.constants as c
#import astropy.units as u
from scipy.integrate import quad
import matplotlib.pyplot as plt

c = 300000. # km/s
H0 = 70. # km/s/Mpc

dh = c/H0
#om = 0.3
#ol = 0.7
#ok = 1-om-ol

#Ez = np.sqrt(om*(1+z)**3+ol+ok*(1+z)**2)

#dc = dh*


def Ez(z, i, j, k):  # Equation 1/E(z)
    return 1. / np.sqrt(i * (1. + z) ** 3 + j + k * (1. + z) ** 2)



result2 = []
z = np.arange(0,10,0.01)
print z
################ Prob 3. Find OmegaM & OmegaL ####################
om = (0.0, 0.3, 0.7, 1.0, 1, 0.3, 0)
ol = (1, 0.7, 0.3, 0, 1, 0, 0)
lb = ('Flat','Flat','Flat','Flat, EdS','Closed','Open','Milne')
ls = ('-','-','-','-','--','--','-.')
m=0
for i in om:  # Make 0.01 step OmegaM between 0~1
    #for j in (0.3,0.7):  # Make 0.01 step OmegaL between 0~1

        j= ol[m]
        k = 1 - i - j

        l = 0
        result = []
        while l < len(z):

            ez, err = quad(Ez, 0., z[l], args=(i, j, k))  # Calculate integral(1/E(z))

            if k > 0:

                dm = dh / np.sqrt(k) * np.sinh(np.sqrt(k) * ez * dh / dh)

            elif k == 0:

                dm = ez * dh

            elif k < 0:

                dm = dh / np.sqrt(np.abs(k)) * np.sin(np.sqrt(np.abs(k)) * ez * dh / dh)
            result.append(dm/(1 + z[l])/dh)
            l = l + 1


            print l

        plt.plot(z,result,ls=ls[m],label='$\Omega_M$=%s, $\Omega_{\Lambda}=%s$ (%s) ' % (i,j,lb[m]))
        m = m + 1
#print result2
#plt.plot(result)
#plt.text(3,0.05,'$\Omega_M+\Omega_{\Lambda}$=1')
plt.legend(loc=0)
plt.xlabel('redshift z')
plt.ylabel('Angular size ') # [c/H$_0$]$^{-1}$
#plt.ylim(0,10)
#plt.xscale('log')
#plt.yscale('log')
plt.show()
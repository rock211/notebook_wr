from scipy.optimize import leastsq
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('cluster_member.dat')
z = data[:,3]
z_err = data[:,4]
n, bins, patches = plt.hist(z,bins=30, alpha=0.4, histtype='stepfilled')
n_bins = np.zeros(len(bins)-1)

for i in range(len(bins)-1):
	n_bins[i] = (bins[i+1] + bins[i])/2.

def gaussian(x, coeffs):
	return coeffs[0] + coeffs[1] * np.exp( - ((x-coeffs[2])/coeffs[3])**2 )
x0 = np.array([1, 1, 1, 1])

def residual(coeffs, y, x):
	return y - gaussian(x, coeffs)

xx, flag = leastsq(residual, x0, args = (n, n_bins))
print xx
x = x = np.linspace(min(z), max(z), 150)
y = gaussian(x, xx)
print np.mean(z), np.std(z)

plt.plot(x,y,'r-', linewidth=2)
plt.text(0.172, 18, '# $of$ $bins$ $=$ $%s$ \n$mean$ $=$ $%s$ \n$std$ $=$ $%s$' %(len(n), round(xx[2],4), round(xx[3]*(-1.),4)), fontsize=16)
plt.title('$Cluster$ $Redshift$', fontsize=16)
plt.ylabel('$Number$', fontsize=15)
plt.xlabel('$redshift$ $z$', fontsize=15)
plt.show()


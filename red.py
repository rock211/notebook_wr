import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.mlab as mlab
from scipy.stats import norm

###### zshift ######
data = np.loadtxt('cluster_member.dat')
z = data[:,3]
z_err = data[:,4]
n, bins, patches = plt.hist(z,bins=20, alpha=0.5, histtype='stepfilled')
mu,sigma = norm.fit(z)
x = np.linspace(min(z), max(z), 150)
y = mlab.normpdf(x,mu,sigma)
print n, bins, patches
'''
###### velocity #####
C0 = 299792 #km/s
vel = C0*(((1+z)**2-1)/((1+z)**2+1)) # km/s
mu_vel, sigma_vel = norm.fit(vel)
#n_v, bins_v, patches_v = plt.hist(vel, bins=20, alpha=0.5, histtype='stepfilled')
x_v = np.linspace(min(vel), max(vel), 150)
y_v = mlab.normpdf(x_v,mu_vel,sigma_vel)
'''
###### plotting ######
plt.plot(x,y*0.5,'r-',linewidth=2)
plt.text(0.172, 75, '$mean$ $=$ $%s$ \n$std$ $=$ $%s$' %(round(mu,3), round(sigma,3)), fontsize=16)
plt.xlabel('$redshift$', fontsize=15)
plt.title('$Cluster$ $Redshift$', fontsize=16)
plt.ylabel('$Number$', fontsize=15)

#####################

'''
plt.plot(x_v, y_v,'r-',linewidth=2)
plt.text(0.172, 75, '$mean$ $=$ $%s$ \n$std$ $=$ $%s$' %(round(mu,3), round(sigma,3)), fontsize=16)
plt.xlabel('$velocity$', fontsize=15)
plt.title('$Cluster$ $Redshift$', fontsize=16)
plt.ylabel('$Number$', fontsize=15)
'''
plt.show()

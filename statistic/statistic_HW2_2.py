import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM

data = np.genfromtxt('xy.dat')

x = data[:,0]
y = data[:,1]
x_err = data[:,2]
y_err = data[:,3]


plt.scatter(x,y,facecolor='white',edgecolors='b',linewidths=1)
plt.errorbar(x,y,x_err,y_err,color='b',alpha=0.5,capsize=2,fmt='None')
plt.xlabel(r'x')
plt.ylabel(r'y')
plt.title(r'Data from xy.dat')
plt.show()

print data.shape

size = 1000
chi_t = np.zeros((size,size))
chi_y = np.zeros((size,size))

for a in np.arange(-1,0,1./size): # y_intercept # if y only -1~0, tot -1.5~-0.5
    for b in np.arange(1,2,1./size): # slope # if y only 1 ~2, tot 1.5~2.5
        chis = []
        chis_y =[]
        for i in range(len(x)) :
            resi = y[i]-b*x[i]-a
            #err_t = np.sqrt(y_err[i]**2+(b**2)*x_err[i]**2)
            chisq_y = resi**2/y_err**2
            #chisq_t = resi**2/err_t**2
            #chis.append(chisq_t)
            chis_y.append(chisq_y)
            #print i,j,chis
        #plt.scatter(z,DM)
        #plt.show()
        #print chis
        #print a,b, np.sum(chis)
        #chi_t[int((a+1.5)*size),int((b-1.5)*size)]=np.sum(chis)
        chi_y[int((a + 1) * size), int((b - 1) * size)] = np.sum(chis_y)
    print a
print chi_t, len(chi_t)

#plt.imshow(chi_t,origin='lower')
#plt.colorbar()
#np.savetxt('chi2_1000_t.txt',chi_t)
np.savetxt('chi2_1000_y.txt',chi_y)
#plt.show()
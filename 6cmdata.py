import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from astropy.io import fits
base = 'F:/fits_bum/'
hdul = fits.open(base+'NGC4522_I.fits')
#print hdul.info()
zz = hdul[0].data[0][0][730:780,730:780]
xx = np.array(range(50))
yy = np.array(range(50))

xii, yii = np.linspace(xx.min(), xx.max(), 100), np.linspace(yy.min(), yy.max(), 100)
xii, yii = np.meshgrid(xii, yii)
rbf = scipy.interpolate.Rbf(xx, yy, zz, function='linear')
zii = rbf(xii, yii)

#print zz.shape
plt.imshow(zii, vmin=zz.min(), vmax=zz.max(), origin='lower',cmap='RdBu')
#plt.imshow(totI,origin='lower',extent=[0,100,0,100])
#plt.show()
data = np.genfromtxt('6cmdata.txt')
#print data[:,2].shape
x = data[:,0]
y = data[:,1]
z = data[:,2]/1000000.
xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)
print zi
plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',alpha=0.5)
           #extent=[x.min(), x.max(), y.min(), y.max()])
#print zi.shape
#plt.contour([xi,yi],zi)
plt.scatter(x, y, c=z,cmap='BuPu',s=1.5)
plt.colorbar()
plt.show()
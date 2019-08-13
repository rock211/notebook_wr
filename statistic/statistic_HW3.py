import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.mplot3d import axes3d

chi = np.genfromtxt('chi_1000.txt')
chi = (chi.T)/2
h = 1 # correspond to 0.001
omM = 279 # correspond to 0.279
omL = 724 # corrspond to 0.724
#print chi.shape

print chi.shape
print chi[724,279]
chi_f=[]
for i in range(len(chi)):
    chi_f.append(chi[999-i,i])

#deriv = (chi_f[omM+h]-2*chi_f[omM]+chi_f[omM-h])/(h**2)
#print np.sqrt(1/(deriv*1000*1000))
#print 1/(deriv*1000*1000)

'''
pp1 = np.exp(-np.array(chi_f))/np.sum(np.exp(-np.array(chi_f)))
#print np.sort(pp)
pp = np.sort(pp1)
ppp =[]
pppp=0
for i in range(len(pp)):
    pppp += np.sum(pp[len(pp)-i-1])
    if 0.7 > pppp >= 0.679 :
        sig1 = pp[len(pp)-i-1]
        print '1sigma',pp[len(pp)-i-1]
    if 0.953 >pppp >= 0.949:
        sig2 = pp[len(pp)-i-1]
        print '2sigma',pp[len(pp) - i - 1]
    if 0.991 > pppp >=0.99:
        sig3 = pp[len(pp)-i-1]
        print '3sigma',pp[len(pp) - i - 1]
    ppp.append(pppp)
xx = np.linspace(0,0.999,1000)
plt.plot(xx,pp1)
plt.axvline(np.argmax(pp1)/1000.,c='k',lw=1)
plt.axvline((2*np.argmax(pp1)-np.argwhere(pp1==sig1))/1000.,c='r',linewidth=0.5,ls='--')
plt.axvline(np.argwhere(pp1==sig1)/1000.,c='r',linewidth=0.5,ls='--')
plt.axvline((2*np.argmax(pp1)-np.argwhere(pp1==sig2))/1000.,c='g',linewidth=0.5,ls='--')
plt.axvline(np.argwhere(pp1==sig2)/1000.,c='g',linewidth=0.5,ls='--')
plt.axvline((2*np.argmax(pp1)-np.argwhere(pp1==sig3))/1000.,c='b',linewidth=0.5,ls='--')
plt.axvline(np.argwhere(pp1==sig3)/1000.,c='b',linewidth=0.5,ls='--')
s1=(np.argmax(pp1)-np.argwhere(pp1==sig1))/1000.
s2=(np.argmax(pp1)-np.argwhere(pp1==sig2))/1000.
s3=(np.argwhere(pp1==sig3)-np.argmax(pp1))/1000.

plt.text(0.35,0.03,r'red : $1\sigma$=%s'%float(s1))
plt.text(0.35,0.0275,r'greed : $2\sigma$=%s'%float(s2))
plt.text(0.35,0.025,r'blue : $3\sigma$=%s'%float(s3))
plt.xlim(0.2,0.4)
plt.xlabel(r'$\Omega_M$',fontsize=15)
plt.ylabel(r'Probability')
plt.show()
'''
'''
#################### Multivariate Gaussian #######################
plt.figure(figsize=(8,8))
mibun2M = (chi[omL,omM+h]-2*chi[omL,omM]+chi[omL,omM-h])/(h**2)
mibun2L = (chi[omL+h,omM]-2*chi[omL,omM]+chi[omL-h,omM])/(h**2)
mibun22 = ((chi[omL+h,omM+h]-chi[omL-h,omM+h])-(chi[omL+h,omM-h]-chi[omL-h,omM-h]))/(4*h**2)
Hess = np.matrix([[mibun2M,mibun22],[mibun22,mibun2L]])*1000*1000
sig = np.linalg.inv(Hess)
print 'Hessian =', Hess
print 'sigma = ',sig
x=np.linspace(0,0.999,1000)-omM/1000.
y=np.linspace(0,0.999,1000)-omL/1000.
O_mm,O_ll = np.meshgrid(x,y)
f_xy = 1/(2*np.pi*np.sqrt(np.linalg.det(sig)))*np.exp(-0.5*(O_mm**2*Hess[0,0] + O_mm*O_ll*Hess[1,0] + O_mm*O_ll*Hess[0,1] + O_ll**2*Hess[1,1]))
#plt.imshow(f_xy,cmap='RdYlBu_r',origin='lower',extent=[0,1,0,1])
fv_f = np.ravel(f_xy); fv_fs = np.sort(fv_f)[::-1]
fv_sum = np.sum(fv_f); fv_new_m = f_xy/fv_sum
fv_news = fv_fs/fv_sum
su = np.zeros(len(fv_f))

for i in range(len(fv_f)):
	J = fv_fs[0:i+1]
	su[i] = np.sum(J)
	if su[i] < fv_sum*0.68 :
		num_sig1 = i
	elif su[i] < fv_sum*0.95 :
		num_sig2 = i
	elif su[i] < fv_sum*0.99 :
		num_sig3 = i

cvel_m = (fv_news[num_sig3], fv_news[num_sig2], fv_news[num_sig1])

plt.contour(fv_new_m, levels=cvel_m, colors ='royalblue',extent=[0,1,0,1])
plt.axvline(omM/1000., 0,1, linestyle='--',linewidth=1.5, color='brown',alpha=0.7)
plt.axhline(omL/1000., 0,1, linestyle='--',linewidth=1.5,color='brown',alpha=0.7)
plt.title('Multivariate Gaussain')
plt.xlabel(r'$\Omega_M$',fontsize=15)
plt.ylabel(r'$\Omega_{\Lambda}$',fontsize=15)
plt.arrow(0,0.999,0.999,-0.999,lw=0.5,ls='-')
plt.arrow(0,0,0.999,0.5,lw=0.5,ls='-')
plt.text(0.83,0.2,'Flat')
plt.text(0.25,0.21,'Accelerating',rotation =30)
plt.text(0.25,0.16,'Decelerating',rotation =30)
plt.savefig('multi_G_2.png',dpi=500)
#plt.show()
###################################################################
'''
################### sampling method ####################
chi_sm = chi - np.amin(chi)

chi_sml = np.exp(-chi_sm)

fv_f = np.ravel(chi_sml); fv_fs = np.sort(fv_f)[::-1]
fv_sum = np.sum(fv_f); fv_new_s = chi_sml/fv_sum
fv_news = fv_fs/fv_sum ;
print 'value',np.argmin(fv_new_s), np.argmax(fv_new_s)
su = np.zeros(len(fv_f))

fig = plt.figure()
ax = fig.gca(projection='3d'); xx=np.linspace(0,0.999,1000); yy=np.linspace(0,0.999,1000)
X,Y = np.meshgrid(xx,yy)
ax.plot_surface(X,Y,fv_new_s,cmap='RdYlBu')
ax.set_xlabel(r'$\Omega_M$')
ax.set_ylabel(r'$\Omega_\Lambda$')
plt.show()

for i in range(len(fv_f)):
        J = fv_fs[0:i+1]
        su[i] = np.sum(J)
        if su[i] < fv_sum*0.68 :
                num_sig1 = i
        elif su[i] < fv_sum*0.95 :
                num_sig2 = i
        elif su[i] < fv_sum*0.99 :
                num_sig3 = i
cvel_s = (fv_news[num_sig3], fv_news[num_sig2], fv_news[num_sig1])
print cvel_s
plt.figure()
plt.figure(figsize=(8,8))
plt.contour(fv_new_s, levels=cvel_s, colors ='crimson',extent=[0,1,0,1])
plt.axvline(omM/1000., 0,1, linestyle='--',linewidth=1.5, color='brown',alpha=0.7)
plt.axhline(omL/1000., 0,1, linestyle='--',linewidth=1.5,color='brown',alpha=0.7)
plt.xlabel(r'$\Omega_M$',fontsize=15)
plt.ylabel(r'$\Omega_{\Lambda}$',fontsize=15)
plt.title('Sampling Method')
plt.arrow(0,0.999,0.999,-0.999,lw=0.5,ls='-')
plt.arrow(0,0,0.999,0.5,lw=0.5,ls='-')
plt.text(0.83,0.2,'Flat')
plt.text(0.25,0.21,'Accelerating',rotation =30)
plt.text(0.25,0.16,'Decelerating',rotation =30)
#c1=plt.Circle((np.argmax(pp1)/1000.,1-np.argmax(pp1)/1000.),(np.argmax(pp1)-np.argwhere(pp1==sig1))/1000.,fc=None,ec='r',linewidth=0.5)
#c2=plt.Circle((np.argmax(pp1)/1000.,1-np.argmax(pp1)/1000.),(np.argmax(pp1)-np.argwhere(pp1==sig2))/1000.,fc=None,ec='g',linewidth=0.5)
#c3=plt.Circle((np.argmax(pp1)/1000.,1-np.argmax(pp1)/1000.),(np.argmax(pp1)-np.argwhere(pp1==sig3))/1000.,fc=None,ec='b',linewidth=0.5)
#ax.add_artist(c1); ax.add_artist(c2); ax.add_artist(c3)
#plt.savefig('sampling.png',dpi=400)
#plt.show()

########################################################################
plt.figure(figsize=(8,8))
plt.contour(fv_new_m, levels=cvel_m, colors ='royalblue',extent=[0,1,0,1])
plt.contour(fv_new_s, levels=cvel_s, colors ='crimson',extent=[0,1,0,1])
plt.axvline(omM/1000., 0,1, linestyle='--',linewidth=1.5, color='brown',alpha=0.7)
plt.axhline(omL/1000., 0,1, linestyle='--',linewidth=1.5,color='brown',alpha=0.7)
plt.xlabel(r'$\Omega_M$',fontsize=15)
plt.ylabel(r'$\Omega_{\Lambda}$',fontsize=15)
plt.title(r'Gaussian vs. Sampling Method')
plt.arrow(0,0.999,0.999,-0.999,lw=0.5,ls='-')
plt.arrow(0,0,0.999,0.5,lw=0.5,ls='-')
plt.text(0.83,0.2,'Flat')
plt.text(0.25,0.21,'Accelerating',rotation =30)
plt.text(0.25,0.16,'Decelerating',rotation =30)
#plt.savefig('vs2.png',dpi=400)
#plt.show()

'''
angle = 180/np.pi*np.arctan((2*abs(sig[0,1]))/(sig[0,0]-sig[1,1]))/2
print angle
plt.figure()
ax = plt.gca()
ell = Ellipse(xy=(omM,omL),width=2*np.sqrt(sig[0,0])*1.52*1000,height=2*np.sqrt(sig[1,1])*1.52*1000,angle=angle)
plt.contour(chi)
plt.imshow(chi,origin='lower',cmap='RdYlBu_r')
ax.add_patch(ell)
plt.xticks([0,200,400,600,800,1000],[0,0.2,0.4,0.6,0.8,1])
plt.yticks([0,200,400,600,800,1000],[0,0.2,0.4,0.6,0.8,1])
plt.xlabel(r'$\Omega_M$',fontsize=15)
plt.ylabel(r'$\Omega_{\Lambda}$',fontsize=15)
plt.arrow(0,999,999,-999,lw=0.5)
plt.show()
'''
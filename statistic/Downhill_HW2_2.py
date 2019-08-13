import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM

data = np.genfromtxt('xy.dat')
#chi =np.genfromtxt('chi2_1000_y.txt') # chi2_1000_t, chi_sn_100.txt
sn = np.genfromtxt('SCPUnion2.1_mu_vs_z.txt')
z = sn[:,1]
DM = sn[:,2]
DM_err = sn[:,3]

def Chi(i,j) :
    d_L = LambdaCDM(H0=70, Om0=i, Ode0=j)
    d_L = d_L.luminosity_distance(z).value
    DM_m = 5 * np.log10(d_L) + 25
    chisq = ((DM - DM_m) / DM_err) ** 2
    return  np.sum(chisq)


A = [0.3,0.9]
chiA =Chi(A[1],A[0]) ; A.append(chiA)
B = [0.6,0.50]
chiB =Chi(B[1],B[0]) ; B.append(chiB)
C = [0.8,0.7]
chiC =Chi(C[1],C[0]) ; C.append(chiC)
#print A,B,C
alpha = 1; gamma = 2; rho = 0.5; sigma = 0.5
x = np.reshape([A,B,C],(3,3))
x = x[np.argsort(x[:,2])]
#print Chi(0.385,0.105)
#print x

tol = 1e-8
m=0
while np.sqrt((x[0,0]-x[2,0])**2+(x[0,1]-x[2,1])**2) > tol :
    m += 1
    #chi[int((i) * size), int((j) * size)] = np.sum(chis)
    plt.scatter(x[0, 0], x[0, 1], s=5, c='r');
    plt.scatter(x[1, 0], x[1, 1], s=5, c='g');
    plt.scatter(x[2, 0], x[2, 1], s=5, c='b')
    # print x[:,0:2]
    t1 = plt.Polygon(x[:, 0:2])
    plt.gca().add_patch(t1)
    plt.xlim(0,1)
    plt.ylim(0,1)
    #plt.savefig('cosmo_trig_%s,png' % m,dpi=300)
    #plt.close()
    #plt.show()

    if np.sqrt((x[0,0]-x[2,0])**2+(x[0,1]-x[2,1])**2) < tol:
        break

    x_o = ((x[0,0])+(x[1,0]))/2 ; y_o = ((x[0,1])+(x[1,1]))/2 # find center of two small point
    x_r = x_o + alpha*(x_o-(x[2,0])) ; y_r = y_o + alpha*(y_o-(x[2,1])) # find reflection point
    Chi_r = Chi(x_r,y_r)
    #print 'reflection' , x_r, y_r, Chi_r
    #print 'new center' , x_o, y_o
    #chi_r = chi[y_r,x_r] # chi value of reflection point

    if x[0,2] <= Chi_r < x[1,2]: #reflection, if chi_new is intermediate
        print 'reflection', x_r, y_r
        x[2,0]=x_r
        x[2,1]=y_r
        x[2,2]=Chi_r

    elif Chi_r < x[0,2] : # if chi_new is best
        x_ex = x_o + gamma* ((x[2,0])-x_o); y_ex = y_o + gamma * ((x[2,1])-y_o)

        chi_ex = Chi(x_ex,y_ex)

        if chi_ex < Chi_r : # expansion confirmed
            print 'expansion', x_ex, y_ex
            x[2,0]=x_ex
            x[2,1]=y_ex
            x[2,2]=chi_ex
        else :              # expansion rejected, take reflection point
            print 'expansion reject', x_r, y_r
            x[2,0] = x_r
            x[2,1] = y_r
            x[2,2] = Chi_r

    else : # if reflection is worst, then contract

        x_c = (x_o + rho*((x[2,0])-x_o)); y_c = (y_o + rho*((x[2,1])-y_o))
        chi_c = Chi(x_c,y_c)

        if chi_c < x[2,2] : # if contraction is better than worst, confirmed
            print 'contraction confirmed', x_c, y_c, chi_c
            x[2,0]=x_c
            x[2,1]=y_c
            x[2,2]=chi_c

        else: # shrink
            print 'shrink'
            x[1,0] = (x[0,0])+sigma*((x[1,0])-(x[0,0]));
            x[1,1]=(x[0,1])+sigma*((x[1,1])-(x[0,1])) ;
            x[1,2]=Chi(x[1,0],x[1,1])

            x[2,0] = (x[0,0])+sigma*((x[2,0])-(x[0,0]));
            x[2,1]=(x[0,1])+sigma*((x[2,1])-(x[0,1])) ;
            x[2,2]=Chi(x[2,0],x[2,1])

    x = x[np.argsort(x[:, 2])]
    #plt.contour(chi)

    #plt.xlim(0,len(chi));plt.ylim(0,len(chi))
    #plt.savefig('%s.png' % m)
    #plt.show()
    #plt.close()
    #print x
#chival = [chi[A],chi[B],chi[C]]
x = x[np.argsort(x[:,2])]
print x[0,0]+x[0,1]
print x

'''
plt.imshow(chi,origin='lower',cmap='RdYlBu_r')
plt.colorbar()
plt.xticks([0,49,99,149,199],['0.6','0.65','0.7','0.75','0.8'])
plt.ylabel(r'Matter Density')
plt.yticks([0, 49, 99, 149, 199], ['0.2', '0.25', '0.3', '0.35', '0.4'])
plt.xlabel(r'Dark Energy')
plt.clim(560,600)

#plt.contour(chi,colors='white',levels=np.arange(1000,1150,3), alpha=0.5) #, levels=np.arange(560,800,3)
plt.scatter(x[0,0],x[0,1],s=5,c='r');plt.scatter(x[1,0],x[1,1],s=5,c='g');plt.scatter(x[2,0],x[2,1],s=5,c='b')
#print x[:,0:2]
t1=plt.Polygon(x[:,0:2])
plt.gca().add_patch(t1)
'''
'''
############# for HW2_ Problem 2 ###############
print 'if total, b =',x[0,0]/len(chi)+1.5,'a =',x[0,1]/len(chi)-1.5
print 'if y only, b =', x[0,0]/len(chi)+1,'a =',x[0,1]/len(chi)-1
b_t=x[0,0]/len(chi)+1.5
b_y=1.598 #x[0,0]/len(chi)+1
a_t=x[0,1]/len(chi)-1.5
a_y=-0.664 #x[0,1]/len(chi)-1
print data[:,0].min()
xx = np.linspace(data[:,0].min()-2,data[:,0].max()+2,1000)

plt.scatter(data[:,0],data[:,1],facecolor='white',edgecolors='b',linewidths=1.5)
plt.errorbar(data[:,0],data[:,1],data[:,2],data[:,3],fmt='None',color='b',alpha=0.5,capsize=2)
plt.plot(xx,b_t*xx+a_t,label='Total Error : y=%sx%s' % (b_t,a_t),c='r',alpha=0.8)
plt.plot(xx,b_y*xx+a_y,label='Y Error only: y=%sx%s' % (b_y,a_y),c='g',alpha=0.8)
plt.legend()
plt.xlim(-3,3); plt.ylim(-8,6)
#plt.show()
'''
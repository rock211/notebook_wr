import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('xy.dat')
chi =np.genfromtxt('chi2_1000_y.txt') # chi2_1000_t, chi_sn_100.txt
sn = np.genfromtxt('SCPUnion2.1_mu_vs_z.txt')
z = sn[:,1]
DM = sn[:,2]
DM_err = sn[:,3]

A = [70*len(chi)/100,37*len(chi)/100]
chiA =chi[A[1],A[0]] ; A.append(chiA)

B = [36*len(chi)/100,50*len(chi)/100]
chiB =chi[B[1],B[0]] ; B.append(chiB)

C = [50*len(chi)/100,22*len(chi)/100]
chiC =chi[C[1],C[0]] ; C.append(chiC)

alpha = 1; gamma = 2; rho = 0.5; sigma = 0.5
x = np.reshape([A,B,C],(3,3))

print x

tol = 0.5
m = 1
while np.sqrt((x[0,0]-x[2,0])**2+(x[0,1]-x[2,1])**2) > tol :
    m += 1
    chis = []
    for k in range(len(z)):  # calculate chisuqare
        # print chisq
        d_L = LambdaCDM(H0=ho, Om0=i, Ode0=j)
        d_L = d_L.luminosity_distance(z[k]).value
        # print k, z[k],d_L
        DM_m = 5 * np.log10(d_L) + 25

        chisq = ((DM[k] - DM_m) / DM_err[k]) ** 2
        # plt.scatter(z[k],DM_m)
        chis.append(chisq)
        # print i,j,chis
    # plt.scatter(z,DM)
    # plt.show()
    # print chis
    print i, j, np.sum(chis)
    chi[int((i) * size), int((j) * size)] = np.sum(chis)

    if np.sqrt((x[0,0]-x[2,0])**2+(x[0,1]-x[2,1])**2) < tol:
        break

    x = x[np.argsort(x[:,2])]
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
    print x
    x_o = (int(x[0,0])+int(x[1,0]))/2 ; y_o = (int(x[0,1])+int(x[1,1]))/2 # find center of two small point
    x_r = x_o + alpha*(x_o-int(x[2,0])) ; y_r = y_o + alpha*(y_o-int(x[2,1])) # find reflection point
    print 'new center' , x_o, y_o
    #chi_r = chi[y_r,x_r] # chi value of reflection point

    if x_r > len(chi) or x_r < 0 or y_r > len(chi) or y_r < 0 :

        x[1, 0] = x[0, 0] + sigma * (x[1, 0] - x[0, 0]);
        x[1, 1] = x[0, 1] + sigma * (x[1, 1] - x[0, 1]);
        x[1, 2] = x[0, 2] + sigma * (x[1, 2] - x[0, 2])

        x[2, 0] = x[0, 0] + sigma * (x[2, 0] - x[0, 0]);
        x[2, 1] = x[0, 1] + sigma * (x[2, 1] - x[0, 1]);
        x[2, 2] = x[0, 2] + sigma * (x[2, 2] - x[0, 2])
    else:
        chi_r = chi[y_r, x_r]

    if x[0,2] <= chi_r < x[1,2]: #reflection, if chi_new is intermediate
        print 'reflection', x_r, y_r
        x[2,0]=x_r
        x[2,1]=y_r
        x[2,2]=chi_r

    elif chi_r < x[0,2] : # if chi_new is best
        x_ex = x_o + gamma* (int(x[2,0])-x_o); y_ex = y_o + gamma * (int(x[2,1])-y_o)

        chi_ex = chi[y_ex,x_ex]

        if chi_ex < x[2,2] : # expansion confirmed
            print 'expansion', x_ex, y_ex
            x[2,0]=x_ex
            x[2,1]=y_ex
            x[2,2]=chi_ex
        else :              # expansion rejected, take reflection point
            print 'expansion reject', x_r, y_r
            x[2,0] = x_r
            x[2,1] = y_r
            x[2,2] = chi_r

    elif x[1,2] < chi_r: # if reflection is worst, then contract

        x_c = int(x_o + rho*(int(x[2,0])-x_o)); y_c = int(y_o + rho*(int(x[2,1])-y_o))
        chi_c = chi[y_c,x_c]

        if chi_c < x[2,2] : # if contraction is better than worst, confirmed
            print 'contraction confirmed', x_c, y_c
            x[2,0]=x_c
            x[2,1]=y_c
            x[2,2]=chi_c

        else: # shrink
            print 'shrink'
            x[1,0] = int(x[0,0])+sigma*(int(x[1,0])-int(x[0,0]));
            x[1,1]=int(x[0,1])+sigma*(int(x[1,1])-int(x[0,1])) ;
            x[1,2]=int(x[0,2])+sigma*(int(x[1,2])-int(x[0,2]))

            x[2,0] = int(x[0,0])+sigma*(int(x[2,0])-int(x[0,0]));
            x[2,1]=int(x[0,1])+sigma*(int(x[2,1])-int(x[0,1])) ;
            x[2,2]=int(x[0,2])+sigma*(int(x[2,2])-int(x[0,2]))

    #plt.contour(chi)

    plt.xlim(0,len(chi));plt.ylim(0,len(chi))
    #plt.savefig('%s.png' % m)
    #plt.show()
    plt.close()
    #print x
#chival = [chi[A],chi[B],chi[C]]
x = x[np.argsort(x[:,2])]
print x

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

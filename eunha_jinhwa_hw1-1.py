import numpy as np
import matplotlib.pyplot as plt

def rk4(f):
    return lambda t, y, dt: (
            lambda dy1: (
            lambda dy2: (
            lambda dy3: (
            lambda dy4: (dy1 + 2*dy2 + 2*dy3 + dy4)/6
            )( dt * f( t + dt  , y + dy3   ) )
        )( dt * f( t + dt/2, y + dy2/2 ) )
        )( dt * f( t + dt/2, y + dy1/2 ) )
        )( dt * f( t       , y         ) )

ome = np.array([[0.001,1],[0.3,0.7], [0.5,0.5], [1.0,0.0], [6.0,0.0]]) #Matter & D.E array
om = [0.3,0.5,1.0,5.0]
daaa = np.genfromtxt('data4.txt')
print daaa[0:843,1]
#print(ome[0,1]) #check the array

initial =[]
age = []
lb = ('','$\Lambda$CDM','Flat','EdS','Closed')
for i in range(5):

     da = rk4(lambda t, a: ((ome[i,0]/a+ome[i,1]*(a**2)+(1-ome[i,0]-ome[i,1])))**0.5) #Equation with D.E
     results = []
    #if da >= 0.:

     t, a, dt = 1., 1., -0.001
     while a >= 0.:

          results.append([t,a,-da(t, a, dt)])
          t, a = t+dt, a+ da(t, a , dt)

     t, a, dt = 1., 1., 0.001
     while t <= 4.:

          results.append([t,a,da(t, a, dt)])
          t, a = t+dt, a+ da(t, a , dt)


     results.sort()
     data = np.array(results)
     print(data)
     if i ==4:
         data[844:1687,1]=daaa[0:843,1][::-1]

     plt.plot(data[:,0],data[:,1],'-',label=r'$\Omega_M$ = %s, $\Omega_{\Lambda}$ = %s (%s)' % (ome[i,0], ome[i,1],lb[i]))  #plot all in one plot
     #np.savetxt('data%s.txt' % i, data)

     initial.append(data[0,0])  #save initial point for Prob4

     h = data[data > 1. ]
    #print(h)

#print initial
plt.axvline(1,ls='--',c='k')
plt.axhline(1,ls='--',c='k')
plt.xlabel(r"$t/t_H$",fontsize=10)
plt.ylabel("Scale Factor a",fontsize=10)
plt.xlim(0.1,3.)
plt.ylim(0,4)
plt.text(2.2,0.25, '$H_0=70km/s/Mpc$', fontsize = 10)
plt.xticks(np.arange(0.,3.000001,0.3))
plt.legend(loc='best',fontsize=10)
#plt.savefig("hw2_3")
plt.show()

'''


################ Prob 3_Initial Age calculate #################

for init in initial:
 Age = (1-init)*139.8
 print Age

############ Prob 4 Plot Omega_M & Omega_L ############

age = np.loadtxt('data0.txt')  # load OmegaM=0.3, OmegaL=0.7 data
print(age) # check the data

ds = 0.0001

Omega = []

for i in range(len(age)): #Equation without Hubble constant. Because i exclude Hubble constant at 'da' also

    Ome_M = 0.3/(((age[i,2]/ds)**2)*age[i,1])  
    Ome_L = 0.7/(((age[i,2]/ds)**2)/age[i,1]**2)
    Omega.append([Ome_M,Ome_L])

data3 = np.array(Omega)
np.savetxt('data_omega.txt',data3)

################## Prob 4_Future of Universe ####################

for i in range(len(age)):  #Find 'a' in after another Hubble time
    if age[i,0] >= 1.9999:
        n = age[i,1]  # n is 'a' at after another Hubble time
        print(n)
        break

z = -1.+1./n #redshift of another Hubble time
print(z)

for i in range(len(age)): #Find a value of Density parameter after another Hubble time
    if -1.+1./age[i,1] <= z:
        M = 0.3/(((age[i,2]/ds)**2)*age[i,1])
        L = 0.7/(((age[i,2]/ds)**2)/age[i,1]**2)
        print(M)
        print(L)
        break

########plot Prob 4_1########

plt.plot(-1+1/age[:,1],data3[:,0],label=r'$\Omega_M(z)$')
plt.plot(-1+1/age[:,1],data3[:,1],label=r'$\Omega_{\Lambda}(z)$')
plt.plot((-2,0), (0.3,0.3), ls = '--', color = 'k')
plt.plot((-2,0), (0.7,0.7), ls = '--', color = 'k')
plt.arrow(2., 0.3, -1.75, 0., fc="k", ec="k",head_width=0.03, head_length=0.1)
plt.arrow(2., 0.7, -1.75, 0., fc="k", ec="k",head_width=0.03, head_length=0.1)
plt.text(2.15,0.3, '$\Omega_M(z)$=0.3', fontsize = 13)
plt.text(2.15,0.7, '$\Omega_{\Lambda}(z)$ = 0.7', fontsize = 13)
plt.legend(loc='best',fontsize=17)
plt.xlim(-2,10)
plt.xticks(np.arange(-2, 10.1, 1))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.axvline(x=0, ymin = 0.0, ymax = 1.0, ls = '--', color = 'k')
#plt.minorticks_on()
plt.xlabel("z",fontsize=15)
plt.ylabel("Density Parameter",fontsize=15)
plt.savefig("hw2_4_1")
plt.show()

########plot Prob 4_2########

plt.plot(-1+1/age[:,1],data3[:,0],label=r'$\Omega_M(z)$')
plt.plot(-1+1/age[:,1],data3[:,1],label=r'$\Omega_{\Lambda} (z)$')
plt.plot((-2,z), (M,M), ls = '--', color = 'k')
plt.plot((-2,z), (L,L), ls = '--', color = 'k')

plt.arrow(2., 0.2, z-1.8, -0.185, fc="k", ec="k",head_width=0.03, head_length=0.1, shape='left' )
plt.arrow(2., 0.3, z-1.8, -0.265, fc="k", ec="k",head_width=0.03, head_length=0.1, shape='left' )
plt.arrow(2., 0.8, z-1.8, 0.155, fc="k", ec="k",head_width=0.03, head_length=0.1, shape='left' )
plt.text(2.1,0.22, 'z=-0.59095(after another Hubble Time)', fontsize=13)
plt.text(2.1,0.32, '$\Omega_M(z)$=0.02849', fontsize = 13)
plt.text(2.1,0.78, '$\Omega_{\Lambda}(z)$ = 0.97142', fontsize = 13)

plt.legend(loc='best',fontsize=17)
plt.xlim(-2,10)
plt.xticks(np.arange(-2, 10.1, 1))
plt.yticks(np.arange(0, 1.1, 0.1))
plt.axvline(x=z, ymin = 0.0, ymax = 1.0, ls = '--', color = 'k')
#plt.minorticks_on()
plt.xlabel("z",fontsize=15)
plt.ylabel("Density Parameter",fontsize=15)
plt.savefig("hw2_4_2")
plt.show()


################ Prob 5_Age at diff z #################

for i in range(len(age)):
    if age[i,1] >= 0.24999 and age[i,1] <= 0.2501 :
        print((age[i,0]-0.0359)*13.98)

    elif age[i,1] >= 0.33329 and age[i,1] <= 0.3334:
        print((age[i,0]-0.0359)*13.98)

    elif age[i,1] >= 0.49998:
        print((age[i,0]-0.0359)*13.98)
        break
'''

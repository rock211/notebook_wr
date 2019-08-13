from urllib import urlretrieve
import time
import json
import os

'''
bodypath = 'G:/yeongu/'
path = bodypath + 'RPS_4pc_ICM2_newacc/slice'
for i in range(250,477) :
    for j in range(2,4) :
        file_name = 'RPS_4pc_ICM2_newacc.0'+ str(i) + '.scal' + str(j) + '.p'
        #file_name = 'RPS_4pc_ICM2_newacc.0' + str(i) + '.surf' + '.p'
        url = 'http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_4pc_ICM2_newacc/slice/' + file_name
        fullfilename = os.path.join(path, file_name)
        urlretrieve(url, fullfilename)
        #print j,i
        print i
'''
'''
bodypath = 'D:/yeongu/'
path = bodypath + 'RPS_8pc_ICM2_newacc/surf'
for i in range(250,499) :
    #for j in ('.whole','.phase1','.phase2','.phase3','.phase4','.phase5') :

        file_name = 'RPS_8pc_ICM2_newacc.0'+ str(i) + '.surf.p'
        url = 'http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_8pc_ICM2_newacc/surf/' + file_name
        fullfilename = os.path.join(path, file_name)
        urlretrieve(url, fullfilename)
        print i
'''
'''
bodypath = 'G:/yeongu/'
path = bodypath + 'RPS_4pc_ICM2_newacc/starpar'
for i in range(250,501) :
    #for j in ('.whole','.phase1','.phase2','.phase3','.phase4','.phase5') :

        file_name = 'RPS_4pc_ICM2_newacc.0'+ str(i) + '.starpar.vtk'
        url = 'http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_4pc_ICM2_newacc/starpar/' + file_name
        fullfilename = os.path.join(path, file_name)
        urlretrieve(url, fullfilename)
        print i

    #http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_8pc_ICM1/zprof_merged/RPS_8pc_ICM1.phase1.zprof.nc
    #http: // tigress - web.princeton.edu / ~changgoo / RPS_8pc / R8_8pc_metal / surf / R8_8pc_metal.0250.surf.p
'''

model = ('noICM','ICM0','ICM1','ICM2','ICM3')
bodypath = 'G:/yeongu/'
for k in (1,4):
    for i in ('','-icm') :
        for j in ('.phase1','.phase2','.phase3','.phase4','.phase5') :

            file_name = 'RPS_8pc_'+model[k]+'_newacc'+ j+i +'.zprof.nc'
            path = bodypath + 'RPS_8pc_'+model[k]+'_newacc/zprof_merged'
            url = 'http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_8pc_'+model[k]+'_newacc/zprof_merged' + '/' + file_name
            fullfilename = os.path.join(path, file_name)
            urlretrieve(url, fullfilename)
        print model[k],i,j
#http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_8pc_ICM1_newacc/zprof_merged/RPS_8pc_ICM1_newacc.phase1-icm.zprof.nc
#print i


#http://tigress-web.princeton.edu/~changgoo/RPS_8pc/MHD_8pc_new/id0/MHD_8pc_new.0250.vtk
#http://tigress-web.princeton.edu/~changgoo/RPS_8pc/MHD_8pc_new/id1/MHD_8pc_new-id1.0250.vtk

#http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_8pc_ICM1/id0/RPS_8pc_ICM1.0250.vtk
#http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_8pc_ICM1/zprof/RPS_8pc_ICM1.0250.phase1.zprof
#http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_8pc_ICM1/surf/RPS_8pc_ICM1.0250.scal0.p
#http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_8pc_ICM1/starpar/RPS_8pc_ICM1.0250.starpar.vtk
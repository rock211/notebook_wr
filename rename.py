from urllib import urlretrieve

import json
import os
import sys

bodypath = 'F:/yeongu/'
idx = range(1,224) # 1 61 61 122 122 183 183 224

for j in map(str,idx) : #'4','8','12','13'

    #os.mkdir(bodypath + 'RPS_8pc_ICM4/id%s' % j)
    print 'ICM3',j
    for i in range(250,501) :
        os.rename('F:/yeongu/RPS_8pc_ICM0/id%s/RPS_8pc_ICM0.0%s.vtk' % (j,i),'F:/yeongu/RPS_8pc_ICM0/id%s/RPS_8pc_ICM0-id%s.0%s.vtk' % (j,j,i))

        print(i)

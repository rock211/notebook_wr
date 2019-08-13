from urllib import urlretrieve
from time import sleep
import datetime
import json
import os
import sys
#http = urllib3.PoolManager()
#print(sys.version_info[0])
#'0','1','2','3','4','5','6','7','8','9','10','11','12','13'

#sleep(5000)

time = datetime.datetime.now()
print time

bodypath = 'G:/yeongu/'
idx = (5,29) # 5,9,14
#idx = (23,26)
for j in map(str,idx) : #'4','8','12','13'

    #os.mkdir(bodypath + 'RPS_4pc_ICM2_newacc/id%s' % j) # RPS_4pc_ICM1_newacc
    print j
    for i in range(474,477) :
        url1 = 'http://tigress-web.princeton.edu/~changgoo/RPS_8pc/RPS_4pc_ICM2_newacc/id'
        url2 = j
        url3 = '/'
        url4 = 'RPS_4pc_ICM2_newacc'
        url5 = '-id%s' % j
        url6 = '.0'
        url7 = str(i)
        url8 = '.vtk'

        if j=='0':
            url = url1 + url2 + url3 + url4 + url6 + url7 + url8  # id0
            file_name = url4 + url6 + url7 + url8 # id0
        else:
            url = url1 + url2 + url3 + url4+ url5+ url6 + url7 + url8 # + url6 + url7 + url8
            file_name = url4+ url5 + url6 + url7 + url8

        print file_name
        path = bodypath + url4+ '/id%s' % j # + url2
        #http: // tigress - web.princeton.edu / ~changgoo / RPS_8pc / RPS_8pc_ICM2 / id1 / RPS_8pc_ICM2 - id1.0250.vtk
        #http: // tigress - web.princeton.edu / ~changgoo / RPS_8pc / RPS_8pc_ICM2 / id0 / RPS_8pc_ICM2.0250.vtk
        fullfilename = os.path.join(path,file_name)
        urlretrieve(url, fullfilename)
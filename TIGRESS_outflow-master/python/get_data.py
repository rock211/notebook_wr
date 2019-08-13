import sys,glob,os
import urllib.request

# check directory
if not os.path.isdir('../data'): os.mkdir('../data')

file_list=['../data/MHD_4pc_new_joined.0450.vtk',
           '../data/MHD_4pc_new.merged.nc',
           '../data/MHD_4pc_new.hst',
           '../data/MHD_4pc_new.sn',
           '../data/mpdf_hot.nc',
           '../data/mpdf_warm.nc',
           '../data/mpdf_int.nc']
# check file 
data_url = 'http://tigress-web.princeton.edu/~changgoo/TIGRESS_example_data/outflow' 
for file_name in file_list:
    print("checking file: {}".format(file_name),end=' -- ')
    if os.path.isfile(file_name):
        print("found!")
    else:
        print("downloading",end='... ')
        url = file_name.replace('../data',data_url)
        urllib.request.urlretrieve(url, file_name)
        print("complete!")

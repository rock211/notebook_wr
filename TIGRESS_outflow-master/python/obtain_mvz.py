# %load obtain_mvz.py
import yt
from yt import derived_field
from yt.units.unit_object import Unit
from yt.units import kpc, pc, kboltz, mh
import numpy as np
import sys,os
import matplotlib.pyplot as plt
from math import *

import cooling
cf=cooling.coolftn()

#def _temp(field, data):
#    return  (mh*data['gas',"pressure"]/data['gas',"density"]/kboltz)

def _T1(field, data):
        return data["gas","pressure"]/data["gas","density"]*mh/kboltz

def _mu(field, data):
        T1=data["gas","T1"].d
        temp=cf.get_temp(T1)
        return temp/T1

def _temperature(field,data):
        return data["gas","T1"]*data["gas","mu"]
    
yt.add_field(("gas", "T1"), function=_T1, units='K')
yt.add_field(("gas", "mu"), function=_mu, units='')
yt.add_field(("gas", "temperature"), function=_temperature, units='K')


unit_base={"length_unit": (1.0,"pc"),
           "time_unit": (1.0,"s*pc/km"),
           "mass_unit": (2.38858753789e-24,"g/cm**3*pc**3"),
           "velocity_unit": (1.0,"km/s"),
           "magnetic_unit": (5.4786746797e-07,"gauss")}

tfin = 545
tini = 199
time = tini

base_dir='/tigress/changgoo/'
pid='MHD_4pc_new'
data_dir='{}{}/id0/'.format(base_dir,pid)
out_dir='{}{}/MVZ/'.format(base_dir,pid)
while time<=tfin:
    num = int(time)
    filename = '{}{}.{:04d}.vtk'.format(data_dir,pid,num)
    outfile_warm = '{}{}{}_mvz.{:04d}.h5'.format(out_dir,'warm/',pid,num)
    outfile_hot  = '{}{}{}_mvz.{:04d}.h5'.format(out_dir,'hot/',pid,num)
    outfile_int  = '{}{}{}_mvz.{:04d}.h5'.format(out_dir,'int/',pid,num)
    print(outfile_warm)
    if not os.path.isfile(outfile_warm):
        ds = yt.load(filename, units_override=unit_base)
        data = ds.all_data()
 
        cut_warm = data.cut_region(['(obj["temperature"]<=2.e4) & (obj["temperature"]>5050)']) 
        cut_int  = data.cut_region(['(obj["temperature"]>2.e4) & (obj["temperature"]<0.5e6)']) 
        cut_hot  = data.cut_region(['(obj["temperature"]>=0.5e6)']) 
 
        pdf_warm = yt.create_profile(cut_warm,['velocity_z','z'],fields='cell_mass',
                           extrema={'velocity_z':(-200,200),'z':(-3584,3584)},
                           n_bins=(400,1792),logs={'velocity_z':False,'z':False},
                           weight_field=None,fractional=False, units={'velocity_z':"km/s", 'z':'pc'})
 
        pdf_int = yt.create_profile(cut_int,['velocity_z','z'],fields='cell_mass',
                           extrema={'velocity_z':(-500,500),'z':(-3584,3584)},
                           n_bins=(400,1792),logs={'velocity_z':False,'z':False},
                           weight_field=None,fractional=False, units={'velocity_z':"km/s", 'z':'pc'})
 
 
        pdf_hot = yt.create_profile(cut_hot,['velocity_z','z'],fields='cell_mass',
                           extrema={'velocity_z':(-500,500),'z':(-3584,3584)},
                           n_bins=(400,1792),logs={'velocity_z':False,'z':False},
                           weight_field=None,fractional=False, units={'velocity_z':"km/s", 'z':'pc'})
 
        pdf_warm.save_as_dataset(outfile_warm)
        pdf_int.save_as_dataset(outfile_int)
        pdf_hot.save_as_dataset(outfile_hot)
    time = time + 1.

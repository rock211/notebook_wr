import glob
import os
import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm,SymLogNorm,NoNorm,Normalize
import numpy as np

import pyathena.yt_analysis.ytathena as ya
from pyathena import read_starvtk,texteffect,set_units
from pyathena.yt_analysis.scatter_sp import scatter_sp
from pyathena.yt_analysis.plot_projection import *
import pyathena.yt_analysis.plot_projection as pproj

problem_id=('RPS_8pc_noICM_newacc','RPS_4pc_ICM1_newacc','RPS_8pc_ICM2_newacc')##'RPS_8pc_ICM4' #
i=1
base='G:yeongu/'

surfnames=[]
bin = 12
for itime in range(250,472,bin): #250-472 #250-401
    surfnames.append('{}{}/surf/{}.{:04d}.surf.p'.format(base,problem_id[i],problem_id[i],itime))

fig=pproj.plot_projection_all(surfnames,axis='y',runaway=False,iscal=3,norm_factor=0.5,old='no')#,pngfname='D:/yeongu/plots/paperplot/surf_%s_%s.eps' % (problem_id[i], bin))
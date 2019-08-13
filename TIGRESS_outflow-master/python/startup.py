import numpy as np
import pandas as pd
import xarray as xr
import glob,os
import astropy.constants as c
import astropy.units as u

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,SymLogNorm,Normalize

from matplotlib.ticker import MultipleLocator
def toggle_xticks(axes,visible=False):
    plt.setp([ax.get_xticklabels() for ax in axes],visible=visible)
def toggle_yticks(axes,visible=False):
    plt.setp([ax.get_yticklabels() for ax in axes],visible=visible)

def texteffect(fontsize=12):
    try:
        from matplotlib.patheffects import withStroke
        myeffect = withStroke(foreground="w", linewidth=3)
        kwargs = dict(path_effects=[myeffect], fontsize=fontsize)
    except ImportError:
        kwargs = dict(fontsize=fontsize)
    return kwargs

from matplotlib import rc_params_from_file

rc=rc_params_from_file('../python/paperstyle.mplstyle')
mpl.rcParams.update(rc)

from ath_hst import read_w_pandas as hst_reader

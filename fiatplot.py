#!/usr/bin/env python2.7

# ...plots a .fiat file. see also fitsplot.
# usage: python fiatplot.py fiatfile 'mag' 'fwhm'

import sys
import warnings
import numpy as np
from matplotlib import *
from pylab import *

warnings.filterwarnings('ignore', '', FutureWarning)

def get_column_names(fiatfile):
    """
    returns list of column names.
    """
    col_names = []
    comments = 0
    with open(fiatfile) as catalog:
        for line in catalog:
            if line.startswith('#'):
                comments += 1
                if 'ttype' in line.lower():
                    col_names.append(line.split('=')[1].split()[0].lower())
            else:
                return col_names

def fiat_plot(fiatfile, x_col, y_col):
    """
    x_col and y_col are names of the columns to be graphed.
    this function doesn't return anything; it just plots.
    """
    x_data = []
    y_data = []
    cols = get_column_names(fiatfile)
    x_loc = cols.index(x_col)
    y_loc = cols.index(y_col)
    with open(fiatfile) as catalog:
        for line in catalog:
            if not line.startswith('#'):
                x_data.append(float(line.split()[x_loc]))
                y_data.append(float(line.split()[y_loc]))

    scatter( x_data, y_data, s=1, c='blue', edgecolors='none', label='Star')
    xlabel(x_col.capitalize())
    ylabel(y_col.capitalize())
    show()

if __name__ == '__main__':
    if len(sys.argv) == 4: # usr included colnames
        fiat_plot(sys.argv[1], sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 2:
        try:
            fiat_plot(sys.argv[1], 'mag', 'fwhm')
        except ValueError:
            fiat_plot(sys.argv[1], 'mag', 'size')
    else:
        print "usage: fiatplot fiatfile ['xlabel' 'ylabel']"

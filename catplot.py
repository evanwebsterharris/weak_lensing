#!/usr/bin/env python2.7

import sys
import warnings
import numpy as np
from matplotlib import *
from pylab import *

warnings.filterwarnings('ignore', '', FutureWarning)
   
def get_column_names(catfile):
    """
    python catplot.py catfile [colnames=['mag_iso', 'fwhm_image']]

    returns list of column names and number of comment lines.
    input = sextracted catfile
    """
    col_names = []
    with open(catfile) as catalog:
        for line in catalog:
            if line.startswith('#'):
                line = line.split()
                try:
                    test = int(line[1]) # see if it names a column
                    col_names.append(line[2].lower())
                except ValueError:
                    pass
            else:
                return col_names

def cat_plot(catfile, x_col, y_col):
    """
    x_col and y_col are names of the columns to be graphed.
    this function doesn't return anything; it just plots.
    """
    x_data = []
    y_data = []
    cols = get_column_names(catfile)
    x_loc = cols.index(x_col)
    y_loc = cols.index(y_col)
    with open(catfile) as catalog:
       for line in catalog:
           if not line.startswith('#'):
               x_data.append(float(line.split()[x_loc]))
               y_data.append(float(line.split()[y_loc]))

    scatter( x_data, y_data, s=1, c='blue', edgecolors='none', label='Star')
    xlabel('Magnitude')
    ylabel('FWHM')
    show()

if __name__ == '__main__':
    if len(sys.argv) == 1:
        help = ("usage: python catplot.py catfile "
                "[colnames=['mag_iso', 'fwhm_image']]")
        print help
    elif len(sys.argv) == 3: # usr included colnames
        cat_plot(sys.argv[1], *sys.argv[2])
    else:
        cat_plot(sys.argv[1], 'mag_iso', 'fwhm_image')

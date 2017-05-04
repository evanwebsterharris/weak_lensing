# program to remove stars from an image by plotting size of object vs.
# magnitude and removing the characteristic band of objects whose size
# does not change with magnitude.
#
# this version uses binning. 
#
# could be improved by favoring lower FWHM values, because theoretically,
# stars should constitute the lowest value. also could check for consistent
# distribution of stars along entire FWHM line. also could check for highest
# density of object distribution along the magnitude axis---this is the
# "mean."
#
# could pick a FWHM that has fewest empty bins, and then check the vertical
# spread of this range---check if that spread is even.
#
# need also to remove FWHM outliers, and possible magnitude outliers too.
# should pick out just relevant elements of object arrays---we don't need
# e.g. the object lens magnitude.
#
# should weight importance with how close mean is to step.
#
# also won't generalize too well---assumes all headers are the same for
# example

import sys
import glob
import subprocess
import fileinput
import warnings
import numpy as np
from scipy.stats import linregress
from matplotlib import *
from pylab import *

warnings.filterwarnings('ignore', '', FutureWarning)

def get_stars(input_file, accuracy=0, FWHM_stepsize=0.1,
              create_file=False, create_galaxy_file=False,
              plot_objects=False, overwrite=False):
    """ 
    With default parameters, returns a .cat file containing the stars in an
    image. Uses SEXtractor to pick out objects, unless a .cat file is
    passed as an argument in the first place.

    Parameters
    ----------
    input_file: can be either .fits or .cat file
    accuracy: if we are sextracting, accuracy gives the number of objects
        to skip, for efficiency. higher = less accurate.
    FWHM_stepsize: how small a stepsize to use when picking out stars.
        The default is probably good.
    number_of_bins: how many bins to divide the image into.
    create_file: whether to create a new .cat file of the stars.
        (default = True)
    create_galaxy_file: whether to create a new .cat file of the galaxies.
        (default = False)
    plot_objects: shows a nice graph of the objects.

    Returns
    -------
    star_list = python list of stars as sextractor objects.
    if create_file==True, the function also creates a new file with
    the stars, in the form of baseimagename_stars.cat. If you want to make
    this file compatible with dlscirc, run 'sex2fiat baseimagename_stars.cat
    > baseimagename.star'
    """

    base_name = ('.').join(input_file.split('.')[:-1])
    warnings.filterwarnings('ignore', '', RuntimeWarning) # for scipy

    if input_file.split('.')[-1] == 'fits': # i.e. input is fitsfile
        if not base_name + '.cat' in glob.glob('*'):
            process = subprocess.Popen('sex.py {0}'.format(input_file))
            process.wait()
        else:
            print "{0} file already exists; skipping sextraction...".format(
                base_name + '.cat')
            cat_file = base_name + '.cat'
        
    else:
        print ".cat file already exists; skipping sextraction..."
        cat_file = input_file

    with open(cat_file) as catalog:
        i = 1
        objects = []
        for line in catalog:
            i += 1
            if line.startswith('#'):
                line_list = line.split()
                if line_list[2] == 'MAG_ISO':
                    m = int(line_list[1]) - 1 # find location of mag
                if line_list[2] == 'MAG':
                    m = int(line_list[1]) - 1 # find location of mag
                elif line_list[2] == 'FWHM_IMAGE':
                    f = int(line_list[1]) - 1 # find location of fwhm
                elif line_list[2] == 'SIZE':
                    f = int(line_list[1]) - 1 # find location of fwhm
                elif line_list[2] == 'NUMBER':
                    n = int(line_list[1]) - 1 # find location of obj. number
                elif line_list[2] == 'CLASS_STAR':
                    cs = int(line_list[1]) - 1 # etc.
                elif line_list[2] == 'X_IMAGE':
                    x_image = int(line_list[1]) - 1
                elif line_list[2] == 'X':
                    x_image = int(line_list[1]) - 1
                elif line_list[2] == 'Y_IMAGE':
                    y_image = int(line_list[1]) - 1
                elif line_list[2] == 'Y':
                    y_image = int(line_list[1]) - 1
                else:
                    pass
            elif i % 2**accuracy == 0:
                obj = line.split()
                objects.append([float(o) for o in obj])
            else:
                continue
    
    mag_min, mag_max = np.percentile([o[m] for o in objects], [0,99])
    FWHM_min, FWHM_max = np.percentile([o[f] for o in objects], [0,99])
    objects = [o for o in objects if (o[m] > mag_min and o[m] < mag_max)\
               and (o[f] > FWHM_min and o[f] < FWHM_max)]
    
    sys.stdout.write("We have %i objects.\n" % len(objects))
    
    FWHM_range = numpy.arange(
        *np.percentile([i[f] for i in objects], [0,80]),
        step=FWHM_stepsize)
    print "Range is", FWHM_range[0], "to", FWHM_range[-1], "pixels."
    # ^we assume the star band is in the lower 80th percentile
    
    all_mags = [o[m] for o in objects]
    all_FWHM = [obj[f] for obj in objects]
    mean_mag = np.mean(all_mags)
    mean_std = np.std([o[m] for o in objects])
    test_mag = np.percentile(all_mags, 50)
    possible_stars = [o for o in objects if o[m] < (test_mag - 0.5*mean_std)
                     ]

    print "{0} possible stars.".format(len(possible_stars))
    # ^hmm

    number_of_bins = 2*int(sqrt((float(len(objects)))))
    # (more objects => more bins)
    print "{0} bins.".format(number_of_bins)

    tol = round(float(
            np.diff(np.percentile([i[f] for i in objects], [25,75])))/14, 2)
     #^FWHM tolerance for counting an object as "close1
    
    print "Tolerance is", tol, "pixels."

    #possible_stars = objects
    bins = list(numpy.linspace(
        min([i[m] for i in possible_stars]),\
        max([i[m] for i in possible_stars]),
        number_of_bins)) # number_of_bins is up to you. using lists for now.
        # ^we are looking in the lower 50th percentile for stars.
    
    number_of_empty_bins = number_of_bins # to start, assume all bins empty
    
    for step in FWHM_range:
        close_objects = [
            o for o in possible_stars if np.abs(o[f] - step) < tol]

        if close_objects: # i.e. if there are any objects within the tol
            star_list = [
                o for o in close_objects if o[m] < (mean_mag - 0.5*mean_std)
                ]

            #if len(star_list) > 1:
            #    slope, ipt, res, pp, slope_std = linregress(
            #        [o[m] for o in star_list],
            #        [o[23] for o in star_list])
            #else:
            #    slope = 0
        
            # ^this whole block weighs how flat the star band is.

            if len(star_list) > 1:
                magnitudes = [s[m] for s in star_list]
                #sizes = [o[f] for o in close_objects]
                indices = np.digitize(magnitudes, bins)
                all_bins = range(1,len(bins))
                empty_bins = [b for b in all_bins if b not in indices]

                star_range = (star_list[0][m] + star_list[-1][m])
                middle = int((star_range)/2)
                mean_starmag = [
                    np.mean([s[m] for s in star_list[:middle]]),
                    np.mean([s[m] for s in star_list[middle:]])]
                mean_starsize = [
                    np.mean([s[f] for s in star_list[:int(step)]]),
                    np.mean([s[f] for s in star_list[int(step):]])]
                true_mean_starsize = np.mean([s[f] for s in star_list])

                horizontal_spread = abs(abs(mean_starmag[0] - middle) -\
                                        abs(mean_starmag[1] - middle))/\
                                        abs(all_mags[-1] - all_mags[0])
                vertical_spread = abs(abs(mean_starsize[0] - step) -\
                                      abs(mean_starsize[1] - step))/\
                                      abs(FWHM_range[-1] - FWHM_range[0])
                close_to_center = abs(true_mean_starsize - step)/\
                                      abs(FWHM_range[-1] - FWHM_range[0])
                #print float(horizontal_spread), float(vertical_spread),\
                #     float(len(empty_bins))/number_of_bins

            #if len(empty_bins)*(1 + np.abs(slope))\
            #< number_of_empty_bins: # want a straight star line
                metric = (float(len(empty_bins))/number_of_bins +\
                0.25*horizontal_spread + vertical_spread +\
                close_to_center + step/abs(FWHM_range[-1]))

                if metric < number_of_empty_bins: 
                    number_of_empty_bins = metric # change later
                    target_list = close_objects
                    final_FWHM = step
                    sys.stdout.flush()
                    sys.stdout.write("\r" + ' '*80)
                    sys.stdout.write(
                        "\r" + "Cutoff FWHM is " + str(final_FWHM) +\
                        " pixels.")
                    sys.stdout.flush()
        else: 
            continue
    
    print "\nDone picking stars."
    target_mags = [obj[m] for obj in target_list]
    target_FWHM = [obj[f] for obj in target_list]
    
    average_mag = sum(target_mags)/float(len(target_mags))
    average_FWHM = sum(target_FWHM)/float(len(target_FWHM))
    mag_std = np.std(target_mags)
    FWHM_std = np.std(target_FWHM)
    
    star_list = [o for o in target_list if o[m] < (test_mag - 0.5*mean_std)]
    galaxy_list = [o for o in objects if o[f] > final_FWHM + tol]
    
    star_mags = [s[m] for s in star_list]
    star_FWHM = [s[f] for s in star_list]
    star_locs = zip([str(s[x_image]) for s in star_list], 
                    [str(s[y_image]) for s in star_list])
    # ^ maybe izip
    # ^^there is definitely a better way to do this
    
    galaxy_mags = [g[m] for g in galaxy_list]
    galaxy_FWHM = [g[f] for g in galaxy_list]
    galaxy_locs = zip([str(s[x_image]) for s in galaxy_list], 
                    [str(s[y_image]) for s in galaxy_list])
    #galaxy_indices = [str(int(g[n])) for g in galaxy_list]

    print len(star_list), "stars found."
    print len(empty_bins), "empty bins."
    #target_slope, ipt, res, pp, slope_std = linregress(star_mags,star_FWHM)
    #print "Slope of stars is", "{:.4f}".format(target_slope),\
    #      "pixels, with standard deviation of",\
    #      "{:.4f}".format(slope_std),"\b."
    
    if overwrite == True:
        print "Writing stars to {0}...".format(base_name + '.cat')
        if not cs:
            print "CLASS_STAR column not found."
            sys.exit()
        with open(cat_file, 'r+') as catalog:
             lines = catalog.readlines()
             catalog.seek(0)
             catalog.truncate()
             for line in lines:
                if line.startswith('#'):
                    catalog.write(line)
                    continue
                #elif line.split()[n] in star_indices\
                elif (line.split()[x_image], line.split()[y_image])\
                in star_locs and not line.startswith('#'): 
                    newline = line.split()
                    newline[cs] = '1'
                else:
                    newline = line.split()
                    newline[cs] = '0'

                catalog.write(' '.join(newline) + '\n')

    if create_file == True:
        print "Writing to {0}...".format(base_name + '_stars.cat')
        with open(base_name + '_stars.cat', 'a') as target_file:
            with open(cat_file, 'r') as catalog:
                for line in catalog: # go thru original catalog.
                    if line.startswith('#'):
                        target_file.write(line)
                    #elif line.split()[n] in star_indices:
                    elif (line.split()[x_image], line.split()[y_image])\
                    in star_locs:
                        target_file.write(line)
                    else:
                        continue

    if create_galaxy_file == True:
        print "Writing to {0}...".format(base_name + '_galaxy.cat')
        with open(base_name + '_galaxy.cat', 'a') as target_file:
            with open(cat_file, 'r') as catalog:
                for line in catalog: # go thru original catalog.
                    if line.startswith('#'):
                        target_file.write(line)
                    elif (line.split()[x_image], line.split()[y_image])\
                    in galaxy_locs:
                        target_file.write(line)
                    else:
                        continue

    if plot_objects == True:
        print "Plotting objects..."

        scatter(
           all_mags, all_FWHM, s=1, c='green',
           edgecolors='none', label='Throwaways')
        scatter(galaxy_mags, galaxy_FWHM, c='black', s=1, edgecolors='none',
                label='Galaxies')
        scatter(star_mags, star_FWHM, c='blue', s=1, edgecolors='none',
                label='Stars')
        axhline(y=final_FWHM + tol, c='red', label='Cutoff Line')
        axhline(y=final_FWHM - tol, c='red')
        axvline(x=test_mag, c='black', label='Average Magnitude')
        legend(loc=0)
        xlabel('Isophotal Magnitude')
        ylabel('FWHM')
        show()

    else:
        return [star_list, galaxy_list]
        # i.e. the list of the data of each star

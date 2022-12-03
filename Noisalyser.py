#!/usr/bin/python3
import datetime
import os
import re
import shutil
import sys
import time
from os import listdir
from os.path import isfile, join
import csv

import matplotlib
from scipy.optimize import curve_fit

#backend which prevent figures showing
matplotlib.use('Agg')

if os.name != 'nt':
    print('Not using Windows...')
    print('OS:', os.name)
    matplotlib.use('Agg')
# matplotlib.use('TkAgg') # makes plots the active window
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy import stats
from matplotlib import animation
from FlowCytometryTools import FCMeasurement, PolyGate
# Documentation: https://eyurtsev.github.io/FlowCytometryTools/tutorial.html
from pylab import *
import gc

csv_header = ['Filename','Directory index','Number of data points','Mean','Median','Mode','STD','CV']
csv_content = [csv_header]

print('Greetings!')
# Plotting parameters
xs = 18
legend_args = [1, 20, True]  # location, fontsize, frameon
plt.rcParams['axes.linewidth'] = 2.0


def n_plot(xlab, ylab, xs, ys):
    plt.minorticks_on()
    plt.tick_params(axis='both', which='major', labelsize=ys, direction='in', length=6, width=2)
    plt.tick_params(axis='both', which='minor', labelsize=ys, direction='in', length=4, width=2)
    plt.tick_params(axis='both', which='both', top=True, right=True)
    plt.xlabel(xlab, fontsize=xs)
    plt.ylabel(ylab, fontsize=ys)
    plt.tight_layout()
    return None


# List only the files in a directory
def file_list(mypath):
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    return onlyfiles


# List everything in a folder
def folder_list(mypath):
    A = os.listdir(mypath)
    return A


# Flags found here
flag_path = 1
flag_ftime_g1 = 0
flag_ftime_g2 = 0
flag_savefig = 0
flag_showfig = 0
flag_apply_g1 = 1
flag_apply_g2 = 0

# main_chan = 'GFP-A'
# main_chan = 'BL1-A'
# main_chan = 'BL2-A'

main_chan = input('Please type name of channel you wish to analyse:   ')

if flag_path == 1:
    # Obtain the path of the current directory
    dir_path = r'C:\Users\mc01535\OneDrive - University of Surrey\11NOV2021 MCW BIG GROUP'
    # dir_path = r'C:\Users\mc01535\OneDrive - University of Surrey\11AUG2021 MCW Drive backup\19MAR2021'

    print('Current path: ', dir_path + '\n')
elif flag_path == 2:
    try:
        path = str(input('Please input a path   '))
    except Exception as e:
        print(e)
        sys.exit('Something went wrong, invalid path')
else:
    sys.exit('Invalid path')

files = file_list(dir_path)
files = [x for x in files if '.fcs' in x]
print('ID, name')
for i, j in enumerate(files):
    print(i, '  ', j)
try:
    # option = int(input('Please enter ID of file to load    '))
    option = 0
except Exception as e:
    print(e)
    sys.exit('Something went wrong, invalid input!')

print('Number of files =')
print(i)

print('Name of last file =')
print(j)

print('File selected =')
print(files[option])

print(option)

while option <= i:
    print('Looping...')

    print('Number of files =')
    print(i)

    print('Name of last file =')
    print(j)

    print('File selected =')
    print(files[option])

    print(option)

    datafile = join(dir_path, files[option])
    name = files[option][:-4].replace(' ', '_')
    sample = FCMeasurement(ID='Test Sample', datafile=datafile)
    print('\nChannel names:   ')
    print(sample.channel_names)
    print('\nChannel information:   ')
    print(sample.channels)
    print('\n\n')

    # Transforming the sample
    print('Transforming the sample...')
    # Documentation found here: https://eyurtsev.github.io/FlowCytometryTools/API/FlowCytometryTools.FCMeasurement.transform.html#FlowCytometryTools.FCMeasurement.transform
    tsample = sample.transform('hlog', channels=[main_chan, 'FSC-A', 'SSC-A', 'FSC-A'], b=500.0)

    # plot FSC vs SSC
    print('Plotting FSC vs SSC')
    tsample.plot(['FSC-A', 'SSC-A'], gates=None, bins=100)
    plt.title('Ungated transformed', size=xs)
    n_plot('FSC-A', 'SSC-A', xs, xs)
    if flag_savefig == 1:
        plt.savefig(name + '_' + 'Ungated_transformed.pdf')
    if flag_showfig == 1:
        show()
    close()

    # Perform Gate 1
    if flag_ftime_g1 == 1:
        print('Please use interface to create the correct gate!')
        print('Select x: FSC-A vs y: SSC-A')
        try:
            tsample.view_interactively(backend='wx')
        except Exception as e:
            print(e)
            sys.exit('Something went wrong, wxpython must be installed!')
    gate1_old = PolyGate([(8.083e+03, 4.555e+03), (8.014e+03, 6.210e+03), (7.642e+03, 5.730e+03), (7.382e+03, 4.915e+03),
                      (7.270e+03, 3.548e+03), (7.371e+03, 2.253e+03), (7.749e+03, 2.829e+03), (7.951e+03, 3.716e+03),
                      (7.966e+03, 4.124e+03)], ('FSC-A', 'SSC-A'), region='in', name='gate1')

    gate1 = PolyGate([(7.126e+03, 4.597e+01), (7.741e+03, 1.091e+02), (8.610e+03, 1.387e+02), (7.933e+03, 4.992e+01),
                      (7.880e+03, 5.190e+01)], ('FSC-A', 'FSC-W'), region='in', name='gate1')



    # Applying gate 1
    if flag_apply_g1 == 1:
        g1tsample = tsample.gate(gate1)
    else:
        g1tsample = tsample

    # plot FSC vs SSC
    print('Plotting gated FSC vs SSC')
    tsample.plot(gate1.channels, gates=[gate1], bins=100)
    plt.title('Ungated transformed', size=xs)
    n_plot(gate1.channels[0], gate1.channels[1], xs, xs)
    if flag_savefig == 1:
        plt.savefig(name + '_' + 'Gated_1_transformed.pdf')
    if flag_showfig == 1:
        show()
    close()


    # Perform Gate 2
    if flag_ftime_g2 == 1:
        print('Please use interface to create the correct gate!')
        print('Select x: GFP-A vs y: SSC-A')
        try:
            g1tsample.view_interactively(backend='wx')
        except Exception as e:
            print(e)
            sys.exit('Something went wrong, wxpython must be installed!')
    gate2 = PolyGate([(6.720e+03, 6.804e+03), (2.054e+03, 9.088e+02), (1.807e+04, -4.516e+02), (1.913e+04, 1.088e+04),
                      (5.950e+03, 1.098e+04), (5.950e+03, 1.098e+04)], (main_chan, 'SSC-A'), region='in', name='gate3')

    # Applying gate 2
    if flag_apply_g2 == 1:
        print("applying gate 2")
        g2g1tsample = g1tsample.gate(gate2)
    else:
        print("not applying gate 2")
        g2g1tsample = g1tsample

    # plot GFP-A SSC-A
    print('Plotting gated GFP-A SSC-A')
    tsample.plot(gate2.channels, gates=[gate2], bins=100)
    plt.title('Ungated transformed', size=xs)
    n_plot(gate2.channels[0], gate2.channels[1], xs, xs)
    if flag_savefig == 1:
        plt.savefig(name + '_' + 'Gated_1_transformed.pdf')
    if flag_showfig == 1:
        show()
    close()

    # plot Gate 1 counts of GFPA
    print('Plotting Gate 1 counts of GFPA')
    g1tsample.plot([main_chan], gates=None, bins=100)
    plt.title('Gated 1 transformed', size=xs)
    n_plot(main_chan, 'Counts', xs, xs)
    if flag_savefig == 1:
        plt.savefig(name + '_' + 'CountGated_1_transformed.pdf')
    if flag_showfig == 1:
        show()
    close()

    # plot Gate 1 and 2 counts of GFPA
    print('Plotting Gate 1 and 2 counts of GFPA')
    g2g1tsample.plot([main_chan], gates=None, bins=100)
    plt.title('Gated 1,2 transformed', size=xs)
    n_plot(main_chan, 'Counts', xs, xs)
    if flag_savefig == 1:
        plt.savefig(name + '_' + 'CountGated_1_2_transformed.pdf')
    if flag_showfig == 1:
        show()
    close()

    # plot Time vs FSC
    print('Plotting x: Time vs y: FSC')
    tsample.plot(['Time', 'FSC-A'], kind='scatter', gates=None)
    plt.title('Time series scatter', size=xs)
    n_plot('Time', 'FSC-A', xs, xs)
    if flag_savefig == 1:
        plt.savefig(name + '_' + 'Time_series_scatter.pdf')
    if flag_showfig == 1:
        show()
    close()

    # plot Time vs FSC
    print('Plotting x: Time vs y: SSC')
    tsample.plot(['Time', 'SSC-A'], kind='scatter', gates=None)
    plt.title('Time series scatter', size=xs)
    n_plot('Time', 'SSC-A', xs, xs)
    if flag_savefig == 1:
        plt.savefig(name + '_' + 'Time_series_scatter.pdf')
    if flag_showfig == 1:
        show()
    close()

    print('\n\nNow calculating the statistics')
    A = g2g1tsample.data[[main_chan, 'FSC-A']][:].values
    data_shape = A.shape
    # Number of data points
    data_num = len(g2g1tsample.data[[main_chan]][:].values)
    # Mean Florescence
    data_mean = np.mean(g2g1tsample.data[[main_chan]][:].values)
    # Median Florescence
    data_median = np.median(g2g1tsample.data[[main_chan]][:].values)
    # Mode Florescence
    data_mode = stats.mode(sample.data[[main_chan]][:].values)[0][0]
    # Standard deviation of a population
    data_pop_STD = np.std(g2g1tsample.data[[main_chan]][:].values)
    # %CV (data_pop_STD/MEAN)
    data_CV_pop = np.divide(data_pop_STD, data_mean) * 100.0  # to get percentage
    # Standard deviation of a population
    data_sample_STD = np.std(g2g1tsample.data[[main_chan]][:].values, ddof=1)
    # %CV (data_pop_STD/MEAN)
    data_CV_sample = np.divide(data_sample_STD, data_mean) * 100.0  # to get percentage
    print('Number of data points: ', data_num)
    print('Mean: %s' % float('%.4g' % data_mean))
    print('Median: %s' % float('%.4g' % data_median))
    print('Mode: %s' % float('%.4g' % data_mode[0]))
    print('STD_pop: %s' % float('%.4g' % data_pop_STD))
    print('CV_pop %s' % float('%.4g' % data_CV_pop))
    print('STD_sample: %s' % float('%.4g' % data_sample_STD))
    print('CV_sample: %s' % float('%.4g' % data_CV_sample))

    # Make an array to save to
    Arr = []
    Arr.append('Number of data points: ' + str(data_num) + '\n')
    Arr.append('Mean: ' + str(data_mean) + '\n')
    Arr.append('Median: ' + str(data_median) + '\n')
    Arr.append('Mode: ' + str(data_mode[0]) + '\n')
    Arr.append('STD_pop: ' + str(data_pop_STD) + '\n')
    Arr.append('CV_pop: ' + str(data_CV_pop) + '\n')
    Arr.append('STD_sample: ' + str(data_sample_STD) + '\n')
    Arr.append('CV_sample: ' + str(data_CV_sample) + '\n')
    # Arr.extend(': '+str()+'\n')

    # print(Arr)
    # print()
    # print(str(data_num))
    # print(str(data_mean))
    # print(str(data_median))
    # print(str(data_mode[0]))
    # print(str(data_pop_STD))
    # print(str(data_CV_pop))
    # print(str(data_sample_STD))
    # print(str(data_CV_sample))
    # np.savetxt(files[option][:-4] + 'stats.txt', Arr, fmt='%s')

    Arr2 = []
    Arr2.append(files[option])
    Arr2.append(str(option))
    Arr2.append(str(data_num))
    Arr2.append(str(data_mean))
    Arr2.append(str(data_median))
    Arr2.append(str(data_mode[0]))
    Arr2.append(str(data_pop_STD))
    Arr2.append(str(data_CV_pop))
    Arr2.append(str(data_sample_STD))
    Arr2.append(str(data_CV_sample))

    csv_content.append(Arr2)

    print('Printing Arr2')
    print(Arr2)

    option = (option + 1)
    print(option)
    plt.close('all')
    gc.collect()

#csv_content = [Arr2,Arr2]
csv_file = open("21FEB2022 11NOV22 red YL3-A group", 'w')

csv_writer = csv.writer(csv_file, delimiter=",")

for row in csv_content:
    csv_writer.writerow(row)

csv_file.close()
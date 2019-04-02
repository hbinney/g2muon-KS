from collections import defaultdict
import csv
import sys 
import ROOT as r
import numpy as np
import sys 

from KSutils import *


def main(): 

    if(len(sys.argv)) < 3:
        print('KS.py <tsv file> <num of intervals>')
        sys.exit(1)

    try:
        fileName = str(sys.argv[1])
        intervals = int(sys.argv[2])
    except ValueError as val_error:
        print(val_error)
        sys.exit(1)

    print("building dictionaries.....")
    totalDict, runNumMin, runNumMax = fill_dictionary(fileName)

    print("getting info from root file.....")
    energybins, energydict, energytot, timebins, timedict, timetot = \
        get_arrs_from_root('T0_cluster_hists.root')

    print("making plots.....")
    plot_hist_intrange(energybins, energydict, energytot, 'cluster energy histogram, different integral ranges', 'cluster energy [MeV]', 'energyhist')
    plot_hist_intrange(timebins, timedict, timetot, 'cluster time histogram, different integral ranges', 'cluster time [c.t.]', 'timehist')
    
    runRange = runNumMax - runNumMin
    diff = int(runRange/intervals)
    minmax = []
    for i in range(intervals):
        minmax.append((runNumMin + i*diff, runNumMin + ((i + 1)*diff) - 1))

    for title in totalDict.keys():
        plot_hist_cdf_runrange(totalDict[title], title, minmax)

    #plot_ks_scan(integrals[0], diff, 200, 'T0 integral sum ks stats', 'integral')

if __name__ == '__main__':
    sys.exit(main())

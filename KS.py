from collections import defaultdict
import csv
import sys 
import ROOT as r
import numpy as np
import sys 

from KSutils import *


def main(): 

    if(len(sys.argv)) < 4:
        print('KS.py <tsv file> <root cluster file> <num of intervals> <KS scan run overlap>')
        sys.exit(1)

    try:
        tsvFileName = str(sys.argv[1])
        rootFileName = str(sys.argv[2])
        intervals = int(sys.argv[3])
        overlap = int(sys.argv[4])
    except ValueError as val_error:
        print(val_error)
        sys.exit(1)

    print("building dictionaries.....")
    totalDict, runNumMin, runNumMax = fill_dictionary(tsvFileName)

    print("getting info from root file.....")
    energybins, energydict, energytot, timebins, timedict, timetot = \
        get_arrs_from_root(rootFileName)

    print("making plots.....")
    plot_hist_intrange(energybins, energydict, energytot, 'cluster energy histogram, different integral ranges', 'cluster energy [MeV]', 'energyhist')
    plot_hist_intrange(timebins, timedict, timetot, 'cluster time histogram, different integral ranges', 'cluster time [c.t.]', 'timehist')
    
    runRange = runNumMax - runNumMin
    diff = int(runRange/intervals)
    minmax = []
    for i in range(intervals):
        minmax.append((runNumMin + i*diff, runNumMin + ((i + 1)*diff) - 1))

    for title in totalDict.keys():
        plot_hist_cdf_runrange(totalDict[title], title, runNumMin, runNumMax, intervals)
        plot_ks_scan(totalDict[title], title, runNumMin, runNumMax, intervals, overlap)

if __name__ == '__main__':
    sys.exit(main())

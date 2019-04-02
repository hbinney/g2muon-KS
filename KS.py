from collections import defaultdict
import csv
import sys 
import ROOT as r
import numpy as np
import sys 

from KSutils import *


def main(): 


    #integrals = defaultdict(list)
    #times = defaultdict(list)
    #rmss = defaultdict(list)
    #ctags = defaultdict(list)
    
    print("building dictionaries.....")

    if(len(sys.argv)) < 3:
        print('KS.py <tsv file> <num of intervals>')
        sys.exit(1)
        
    try:
        fileName = str(sys.argv[1])
        intervals = int(sys.argv[2])
    except ValueError as val_error:
        print(val_error)
            
            
    runNumMin = 1000000
    runNumMax = 0

    try:
        with open(fileName, 'r') as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')
            totalDict = defaultdict(lambda: defaultdict(list))
            titles = next(tsvin)
            #number of physics entries (after run num and pulse index)
            #numPhysicsEntries = len(titles)-2
            next(tsvin)
            for row in tsvin:
                runNum = int(row[0])
                pulseNum = int(row[1])
                runNumMin = min(runNumMin, runNum)
                runNumMax = max(runNumMax, runNum)
                #loop over physics entries (not run num or pulse num)
                for physNum in range(2, len(titles)):
                    totalDict[str(titles[physNum])][pulseNum].append((runNum, float(row[physNum])))
                #integrals[pulsenum].append((runnum, float(row[2])))
                #times[pulsenum].append((runnum, float(row[3])))
                #rmss[pulsenum].append((runnum, float(row[4])))
                #ctags[pulsenum].append((runnum, int(row[5])))
    except FileNotFoundError as fnf_error:
        print(fnf_error)

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
        plot_hist_cdf_runrange(totalDict[title], str(title).strip(), minmax)

    #plot_hist_cdf_runrange(integrals, 'T0 integral sum histogram and cdf', 'integral sum [ADC ct]', 'integral', minmax)
    #plot_hist_cdf_runrange(times, 'T0 beam time histogram and cdf', 'time [ct]', 'time', minmax)
    #plot_hist_cdf_runrange(rmss, 'T0 RMS histogram and cdf', 'RMS [ct]', 'RMS', minmax)
    #plot_hist_cdf_runrange(ctags, 'ctag histogram and cdf', 'ctag', 'ctag', minmax)
    
    #plot_ks_scan(integrals[0], diff, 200, 'T0 integral sum ks stats', 'integral')

if __name__ == '__main__':
    sys.exit(main())

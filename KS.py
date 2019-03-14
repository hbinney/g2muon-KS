from collections import defaultdict
import csv
import sys 
import ROOT as r
import numpy as np

from KSutils import *


def main(): 
    integrals = defaultdict(list)
    times = defaultdict(list)
    rmss = defaultdict(list)
    ctags = defaultdict(list)
    
    print("building dictionaries.....")
    with open('endgame_nodqc_longer.tsv', 'r') as tsvin:
        next(tsvin)
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            runnum = int(row[0])
            pulsenum = int(row[1])
            integrals[pulsenum].append((runnum, float(row[2])))
            times[pulsenum].append((runnum, float(row[3])))
            rmss[pulsenum].append((runnum, float(row[4])))
            ctags[pulsenum].append((runnum, int(row[5])))
    
    print("getting info from root file.....")
    energybins, energydict, energytot, timebins, timedict, timetot = \
        get_arrs_from_root('T0_cluster_hists.root')

    print("making plots.....")
    plot_hist_intrange(energybins, energydict, energytot, 'cluster energy histogram, different integral ranges', 'cluster energy [MeV]', 'energyhist')
    plot_hist_intrange(timebins, timedict, timetot, 'cluster time histogram, different integral ranges', 'cluster time [c.t.]', 'timehist')
    
    #minmax = ((16908, 17063), (17064, 17218), (17219, 17373), (17374, 17528))
    #minmax = ((16908, 17130), (17131, 17300), (17301, 17405), (17406, 17528))
    intervals = 8
    diff = int(620/intervals)
    minmax = []
    for i in range(intervals):
        minmax.append((16908 + i*diff, 16908 + ((i + 1)*diff) - 1))

    plot_hist_cdf_runrange(integrals, 'T0 integral sum histogram and cdf', 'integral sum [ADC ct]', 'integral', minmax)
    plot_hist_cdf_runrange(times, 'T0 beam time histogram and cdf', 'time [ct]', 'time', minmax)
    plot_hist_cdf_runrange(rmss, 'T0 RMS histogram and cdf', 'RMS [ct]', 'RMS', minmax)
    plot_hist_cdf_runrange(ctags, 'ctag histogram and cdf', 'ctag', 'ctag', minmax)
    

    plot_ks_scan(integrals[0], 124/2, 200, 'T0 integral sum ks stats', 'integral')

if __name__ == '__main__':
    sys.exit(main())

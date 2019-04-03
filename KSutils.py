import ROOT as r
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import os.path
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import csv
from scipy import interpolate
from matplotlib.lines import Line2D
import sys

plt.rcParams.update({'figure.max_open_warning': 0})

'''Takes default dict of default dict of list and fills it 
with contents of fileName .tsv file, indexing by physics 
parameter and then bunch number.  Assumes .tsv file structure
created by make_textfile.'''  
def fill_dictionary(fileName):

    totalDict = defaultdict(lambda: defaultdict(list))

    runNumMin = 1000000
    runNumMax = 0

    try:
        with open(fileName, 'r') as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')
            titles = next(tsvin)
            for row in tsvin:
                runNum = int(row[0])
                pulseNum = int(row[1])
                runNumMin = min(runNumMin, runNum)
                runNumMax = max(runNumMax, runNum)
                #loop over physics entries (after run num, sub run num, event num, pulse num)
                for physNum in range(4, len(titles)):
                    totalDict[str(titles[physNum])][pulseNum].append((runNum, float(row[physNum])))
            return totalDict, runNumMin, runNumMax
    except FileNotFoundError as fnf_error:
        print(fnf_error)
        sys.exit(1)


'''Function so that histograms can be plotted with no fill
(outlines only).'''
def histOutline(dataIn, thebins):
    (histIn, binsIn) = np.histogram(dataIn, bins=thebins)

    stepSize = binsIn[1] - binsIn[0]

    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0

    return (bins, data)

'''Takes an array of numbers and returns cdf function
and sorted array for plotting'''
def make_cdf(arr):
    x = np.sort(arr)
    f = np.array(range(len(arr)))/float(len(arr))
    return x, f

'''Finds the (x,y) location of the maximum vertical 
difference between two sets of (x,y) coordinates.  An 
approximation of the KS statistic with cdf functions,
used for plotting.  TODO: doesn't always work properly
when the x ranges of the input functions are different.'''
def find_max_diff(x1, f1, x2, f2):
    interp1 = interpolate.interp1d(x1, f1, fill_value="extrapolate")
    interp2 = interpolate.interp1d(x2, f2, fill_value="extrapolate")
    xs = np.linspace(max(min(x1), min(x2)), min(max(x1), max(x2)), 10000)
    f1n = interp1(xs)
    f2n = interp2(xs)
    diffs = [f1-f2 for f1, f2 in zip(f1n, f2n)]
    themax = max(diffs)
    loc = np.argmax(diffs)
    return xs[loc], f1n[loc], f2n[loc]

'''Takes list of tuples (run, val) and run limits and
returns the KS statistic between the full distribution
and the distribution of the value in the range of runs.'''
def find_ks(arr, minrun, maxrun):
    arrtot = [i[1] for i in arr]
    arrlim = [i[1] for i in arr if i[0] >= minrun and i[0] <= maxrun]
    ks, _ = stats.ks_2samp(arrtot, arrlim)
    return ks

'''Takes a dict (indexed by pulse num) of lists of tuples (run, val) 
and plots the ks statistic for a set of run ranges.  Intervals is the 
number of ranges and overlap is the run number overlap between adjacent points.'''
def plot_ks_scan(arr, titlestr, runNumMin, runNumMax, intervals, overlap):
    runRange = runNumMax - runNumMin
    diff = int(runRange/intervals)
    if(titlestr=="T0_int"):
        title = 'T0 integral sum ks stats'
        outstr = 'integral'
    elif(titlestr=="T0_time"):
        title='T0 beam time ks stats'
        outstr = 'time'
    elif(titlestr =="T0_RMS"):
        title = 'T0 RMS ks stats'
        outstr = 'RMS'
    elif(titlestr == "ctag"):
        title = 'ctag ks stats'
        outstr = 'ctag'
    else:
        print('Unknown title string: ' + titlestr)
        sys.exit(1)
    pdf = PdfPages('ksplots/ksscan_'+outstr+'.pdf')
    for pulse in range(8):
        arrpulse = arr[pulse]
        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(1,1,1)
        mins = np.arange(runNumMin, runNumMax, diff)
        xs = range(len(mins))
        kslist = []
        ticks = []
        for min in mins:
            max = min+overlap
            if max > runNumMax:
                max = runNumMax
            ks = find_ks(arrpulse, min, max)
            kslist.append(ks)
            ticks.append(str(int(min))+'-'+str(int(max)))
        plt.xlabel('run range')
        plt.xticks(xs, ticks, rotation=30)
        plt.ylabel('ks stat')
        plt.title('pulse '+str(pulse)+' '+title)
        plt.plot(xs, kslist, 'ko')
        plt.tight_layout()
        pdf.savefig()
        ax.clear()
    pdf.close()

'''Version of plot_hist_cdf for dividing into run ranges.
Takes dictionary of tuples (run, val) arranged by pulse
number, min and max run number, and difference between divisions
and plots duo-plot of histogram and cdfs with ks values for 
a default set of run ranges for each pulse.'''
def plot_hist_cdf_runrange(arr, titlestr, runNumMin, runNumMax, intervals):
    runRange = runNumMax - runNumMin
    diff = int(runRange/intervals)
    minmax = []
    for i in range(intervals):
        minmax.append((runNumMin + i*diff, runNumMin + ((i + 1)*diff) - 1))
    if(titlestr=="T0_int"):
        title = 'T0 integral sum histogram and cdf'
        xstr = 'integral sum [ADC ct]'
        outstr = 'integral'
    elif(titlestr=="T0_time"):
        title='T0 beam time histogram and cdf'
        xstr = 'time [ct]'
        outstr = 'time'
    elif(titlestr =="T0_RMS"):
        title = 'T0 RMS histogram and cdf',
        xstr = 'RMS [ct]'
        outstr = 'RMS'
    elif(titlestr == "ctag"):
        title = 'ctag histgram and cdf'
        xstr = 'ctag'
        outstr = 'ctag'
    else:
        print('Unknown title string: ' + titlestr)
        sys.exit(1)
    fig = plt.figure(figsize=(8,4))
    ax1 = fig.add_subplot(1,1,1)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(1,1,1)
    colors = ['C6', 'C3', 'C1', 'C2', 'C0', 'C4', 'C5', 'C7', 'C8', 'C9']
    legend_elements = list()
    pdf = PdfPages('hists/'+outstr+'_nodqc_all.pdf') 
    ksavg = np.zeros(len(minmax))
    for pulse in range(0, 8):
        figarr, axarr = plt.subplots(1,2, figsize=(8,4))
        arrtot = [i[1] for i in arr[pulse]]
        x, f = make_cdf(arrtot)
        n, thebins, patches = axarr[0].hist(arrtot, bins=100, label='all runs, '+str(len(arrtot))+' events', facecolor = 'silver', density=True)
        axarr[1].plot(x, f, label='all runs, '+str(len(arrtot))+' events', color='silver')
        axarr[1].fill_between(x, 0, f, color='silver')
        axarr[0].set_xlim(min(x), max(x))
        axarr[1].set_xlim(min(x), max(x))
        axarr[1].set_ylim(0, 1)
        for i in range(len(minmax)):
            minrun = minmax[i][0]
            maxrun = minmax[i][1]
            arrlim = [i[1] for i in arr[pulse] \
                          if i[0] >= minrun and i[0] <= maxrun]
            axarr[0].set_yscale('log')
            (bins, n) = histOutline(arrlim, thebins)
            n = n/(sum(n)*(bins[2]-bins[0]))
            axarr[0].plot(bins, n, color=colors[i], label='run '+str(minrun)+'-'+str(maxrun)+', '+str(len(arrlim))+' events')
            xlim, flim = make_cdf(arrlim)
            axarr[1].plot(xlim, flim, label='run '+str(minrun)+'-'+str(maxrun), color = colors[i])
            ks = find_ks(arr[pulse], minrun, maxrun)
            ksavg[i]+=ks/len(minmax)
            ax1.plot(pulse, ks, color = colors[i], marker = 'o', linestyle='-', label='run '+str(minrun)+'-'+str(maxrun))
            if(pulse==0):
                legend_elements.append(Line2D([0], [0], marker='o', color=colors[i], label='run '+str(minrun)+'-'+str(maxrun)))
            if(outstr=='time' or outstr=='RMS'):
                axarr[1].text(0.6, 0.55-0.05*i, 'ks = %.2f'%ks, transform=axarr[1].transAxes, color=colors[i])
            else:
                axarr[1].text(0.2, 0.75-0.05*i, 'ks = %.2f'%ks, transform=axarr[1].transAxes, color=colors[i])
        figarr.suptitle(titlestr+', pulse '+str(pulse))
        axarr[0].set_xlabel(xstr)
        axarr[0].set_ylabel('normalized counts')
        axarr[1].set_xlabel(xstr)
        axarr[1].set_ylabel('cdf')
        handles, labels = axarr[0].get_legend_handles_labels()
        if(outstr=='time' or outstr=='RMS'):
            anchor = (0.86, 0.98)
        else:
            anchor = (0.14, 0.98)
        figarr.legend(handles, labels, loc='upper center', bbox_to_anchor=anchor, \
                          fontsize='x-small', fancybox=True, shadow=True)
        pdf.savefig()
        axarr[0].clear()
        axarr[1].clear()
    xs = range(len(minmax))
    ax2.plot(xs, ksavg, 'ko')
    ax1.legend(handles = legend_elements, loc='upper center', bbox_to_anchor=(0, 1.15), \
                   fontsize='small', fancybox=True, shadow=True)
    ax1.set_ylim(0, 0.4)
    ax1.set_xlabel('pulse num')
    ax1.set_ylabel('ks stat')
    ax1.set_title(outstr+' distribution, KS stat vs pulse num')
    ax2.set_xlabel('subset')
    ax2.set_ylabel('ks stat avg')
    ax2.set_title(outstr+' distribution, KS stat avg vs data subset')
    fig2.savefig('ksplots/'+outstr+'_ksvals_avg_nodqc.png')
    fig.savefig('ksplots/'+outstr+'_ksvals_nodqc.png')
    plt.close(fig)
    plt.close(fig2)
    plt.close(figarr)
    ax1.clear()
    pdf.close()
    ax2.clear()
    
'''Turn root hists into py dictionaries/arrays for use in other functions.
Returns list with histogram bins for time and energy, dictionaries with the 
time/energy for each T0 integral range, and an array with all ranges summed 
together.'''
def get_arrs_from_root(filename):
    f = r.TFile(filename)
    ranges = np.arange(20000, 420000, 20000)
    energyhist = f.Get('T0ClusterHists/energy_'+str(ranges[0]))
    timehist = f.Get('T0ClusterHists/time_'+str(ranges[0]))
    energytot = np.zeros(energyhist.GetNbinsX())
    timetot = np.zeros(timehist.GetNbinsX())
    energydict = defaultdict(list)
    timedict = defaultdict(list)
    energybins = list()
    timebins = list()
    for bin in range(0, energyhist.GetNbinsX()):
        energybins.append(energyhist.GetBinCenter(bin+1))
    for bin in range(0, timehist.GetNbinsX()):
        timebins.append(timehist.GetBinCenter(bin+1))
    for i in range(0, len(ranges)):
        min = ranges[i]
        energyhist = f.Get('T0ClusterHists/energy_'+str(min))
        timehist = f.Get('T0ClusterHists/time_'+str(min))
        clusterenergies = []
        clustertimes = []
        for bin in range(0,energyhist.GetNbinsX()):
            clusterenergies.append(energyhist.GetBinContent(bin+1))
            energytot[bin]+=energyhist.GetBinContent(bin+1)
        for bin in range(0,timehist.GetNbinsX()):
            clustertimes.append(timehist.GetBinContent(bin+1))
            timetot[bin]+=timehist.GetBinContent(bin+1)
        energydict[min] = clusterenergies
        timedict[min] = clustertimes
    return energybins, energydict, energytot, timebins, timedict, timetot

'''Histogram comparison for cluster time/energy with different
integral sum ranges.  Takes dict with structure min integral: 
(list) and array with sum of all integral ranges for comparison.'''
def plot_hist_intrange(bins, dict, arrtot, titlestr, xstr, outstr):
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(1,1,1)
    arrtot = [i/sum(arrtot) for i in arrtot]
    plt.semilogy(bins, arrtot, label='all T0 integrals', color='k', lw=6)
    arr1 = [i/sum(dict[40000]) for i in dict[40000]]
    arr2 = [i/sum(dict[320000]) for i in dict[320000]]
    plt.semilogy(bins, arr1, 'C0', label='T0 integral 40,000-60,000')
    plt.semilogy(bins, arr2, 'C3', label='T0 integral 320,000-340,0000')
    plt.xlabel(xstr)
    plt.ylabel('events')
    plt.title(titlestr)
    ax.legend(fontsize='small', fancybox=True, shadow=True)
    plt.savefig(outstr+'.png')
    ax.clear()


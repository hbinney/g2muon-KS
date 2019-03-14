import ROOT as r
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import os.path
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import csv

def pass_filter(t):
    if(t.integraltubeA+t.integraltubeB > 20000 and \
           (t.timetubeA+t.timetubeB)/2 > 20000 and \
           (t.timetubeA+t.timetubeB)/2 < 25050 and \
           (t.tRMStubeA+t.tRMStubeB)/2 > 20 and \
           (t.tRMStubeA+t.tRMStubeB)/2 < 50 and \
           t.syncPulseTimetubeA > 5300 and \
           t.syncPulseTimetubeA < 5600 and \
           t.syncPulseIntegraltubeA > 8000 and \
           t.ctag > 100 and \
           t.ctag/(t.integraltubeA+t.integraltubeB) > 0.0005):
        return True
    else:
        return False

with open('endgame.txt', 'r') as myfile:
    runs=myfile.read().split('\n')
runstr = 'endgame'

runs = runs[:200]

with open("endgame_dqc_long.tsv", "w") as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    writer.writerow(["run num", "subrun num", "event num", "pulse index", "int sum", "time", "RMS", "ctag"])
    for run in runs:
        print(run)
        if(os.path.isfile(run)):
            f = r.TFile(run)
            try:
                t = f.Get('T0Analyzer/t0Tree')
                num = t.GetEntries()
                for i in range(0,num):
                    if(i%1000==0):
                        print i
                    t.GetEntry(i)
                    if(pass_filter(t)):
                        writer.writerow([t.runNum, t.subRunNum, t.eventNum, t.pulseIndex, t.integraltubeA+t.integraltubeB, \
                                             (t.timetubeA+t.timetubeB)/2, (t.tRMStubeA+t.tRMStubeB)/2, t.ctag])
            except Exception as ex:
                template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                print message



import ROOT as r
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import os.path
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import csv
import sys

runstr = str(sys.argv[1])

try:
    with open(runstr+".txt", "r") as myfile:
        runs=myfile.read().split('\n')
except Exception as ex:
    template = "An exception of type {0} occurred. Arguments:\n{1!r}\nExpected Format of input file is arg.txt"
    message = template.format(type(ex).__name__, ex.args)
    print(message)
    sys.exit()

with open(runstr+".tsv", "w") as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    writer.writerow(["run_num", "subrun_num", "event_num", "pulse_index", "T0_int_sum", "T0_time", "T0_RMS", "ctag"])
    for run in runs:
        print(run)
        if(os.path.isfile(run)):
            f = r.TFile(run)
            try:
                t = f.Get('T0Analyzer/t0Tree')
                num = t.GetEntries()
                for i in range(0,num):
                    if(i%1000==0):
                        print(i)
                    t.GetEntry(i)
                    writer.writerow([t.runNum, t.subRunNum, t.eventNum, t.pulseIndex, t.integraltubeA+t.integraltubeB, \
                                         (t.timetubeA+t.timetubeB)/2, (t.tRMStubeA+t.tRMStubeB)/2, t.ctag])
            except Exception as ex:
                template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                message = template.format(type(ex).__name__, ex.args)
                print(message)
        else:
            print("File "+run+" does not exist\n")

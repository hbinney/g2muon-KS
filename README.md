# g2muon-KS
KS testing for Muon g-2 experiment

Inputs:

- root file containing a TTree which is the output of the T0PulseProcessorAnalyzer

- root file containing histograms which have cluster time and energy info binned by T0 integral.  This is a custom data format which should probably be changed to be more general in the future

Procedure:

- Run make_txtfile.py on the first root file.  This produces a .tsv file with the run number, subrun number, event number, pulse number, and then the physics parameters for each event in columns.

- Run KS.py to make plots from the .tsv file and the second root file. 
  - KS.py \<tsv file\> \<root cluster file\> \<num of intervals\> \<KS scan run overlap\>
  - Num of intervals is the number of different run divisions for the histogram/cdf plot. 
  - KS scan run overlap is the number of overlapping runs between adjacent points in the KS scan.  Currently the scan uses the same num of intervals as the histogram/cdf plots.

Outputs:

- For each physics input, produces histogram/cdf plot and plot of averaged ks (across pulses) using the same subsets as the histogram/cdf plot.  Stored under kssplots/

- For each physics input, produces KS scan with overlapping subsets.  Currently uses the same number of intervals as the histogram/cdf plots.  Stored under kssplots/

- Creates timehist and energyhist comparing cluster energy and time for different T0 integral ranges.

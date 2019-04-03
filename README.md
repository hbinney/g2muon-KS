# g2muon-KS
KS testing for Muon g-2 experiment

Inputs:

  -root file containing a TTree which is the output of the T0PulseProcessorAnalyzer

  -root file containing a TTree which has cluster and energy information

Procedure:

  -Run make_txtfile.py on the first root file.  This produces a .tsv file with the run number, subrun number, event number, pulse number, and then the physics parameters for each event in columns.

  -Run KS.py to make plots from the .tsv file and the second root file.

#!/usr/bin/python
# was: #!/snap/bin/pyroot
# was: #!/usr/bin/python3
# Út 17. září 2024, 13:55:18 CEST

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt
from scipy.stats import skew
import numpy as np

def makeGraph():
    gr = ROOT.TGraphErrors()
    gr.SetMarkerColor(ROOT.kBlack)
    gr.SetLineColor(gr.GetMarkerColor())
    gr.SetLineWidth(2)
    gr.SetLineStyle(1)
    gr.SetMarkerSize(1.)
    gr.SetMarkerStyle(20)
    return gr

# https://stackoverflow.com/questions/61521371/calculate-weighted-statistical-moments-in-python
# but FIXED!
def get_n_weighted_moment(values, ws, n):

    #assert n>0 & (values.shape == weights.shape)
    w_avg = np.average(values, weights = ws)
    w_var = np.sum(ws * (values - w_avg)**2) / np.sum(ws)

    if n==1:
        return w_avg
    elif n==2:
        return w_var
    else:
        w_std = np.sqrt(w_var)
        print('avg, var, std: ', w_avg, w_var, w_std)
        #return np.sum(ws * ((values - w_avg)/w_std)**n)/np.sum(ws)
        cn = np.average( (values - w_avg)**n, weights=ws)
        print(f'c{n}: {cn}')
        return cn

cans = []
stuff = []

##########################################
# https://www.tutorialspoint.com/python/python_command_line_arguments.htm
def main(argv):
    #if len(sys.argv) > 1:
    #  foo = sys.argv[1]

    ### https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    ### https://pymotw.com/2/getopt/
    ### https://docs.python.org/3.1/library/getopt.html
    gBatch = False
    gTag=''
    if len(argv) < 2:
        print('Usage: ')
        print(f'{argv[0]} conex_rootfile.root')
        print('Example:')
        print(f'{argv[0]} simulated_showers/merged/conex_p_E_19_EPOS_merged.root)
        return

    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, g`))

  

    nbx, x1,x2 = 100, 0, 3000
    nby, y1,y2 = 100, 0, 50000000
    hn = 'h2'
    h2 = ROOT.TH2D(hn, ';', nbx, x1, x2, nby, y1, y2)
    h2.SetStats(0)

    filename = argv[1]
    rfile = ROOT.TFile(filename, 'read')
    treename = 'Shower'
    tree = rfile.Get(treename)
    Nevts = tree.GetEntries()
    
    hXmax = ROOT.TH1D('hXmax', ';', 150, 0, 1500)
    hsigma = ROOT.TH1D('hsigmaXmax', ';#sigma_{Xmax}', 200, 0, 500)
    hskew = ROOT.TH1D('hskewXmax', ';skew of Xmax', 200, -1, 3)
  
    for ievt in range(0, Nevts):
        if ievt % 100 == 0: 
            print(f'Processing event {ievt} / {Nevts}')
        ientry = tree.LoadTree(ievt)
        nb = tree.GetEntry(ievt)
        
        Xs = []
        dEdXs = []
        # dE/dX is a weight to each position x of the shower development
        for x,dedx in zip(tree.X,tree.dEdX):
            h2.Fill(x,dedx)
            Xs.append(x)
            dEdXs.append(dedx)
        
        Xmax = Xs[dEdXs.index(max(dEdXs))]
        Xaver = np.average(Xs, weights = dEdXs)
        variance = np.average( (Xs - Xaver)**2, weights = dEdXs )
        print(Xmax, Xaver, variance)
        sigmaXmax = -1.
        skew = -1.
        if variance > 0: 
            sigmaXmax = sqrt(variance)
            skew = get_n_weighted_moment(Xs, dEdXs, 3) / sigmaXmax**3
            print(sigmaXmax,skew)
            skew2 = np.average(((Xs-Xaver)/np.sqrt(variance))**3, weights=dEdXs)
            print(skew, skew2)
        hXmax.Fill(Xmax)
        hsigma.Fill(sigmaXmax)
        hskew.Fill(skew)

    #hnames = ['histo_h', '']
    # hs = []
    # for hname in hnames:
    # 	  h = rfile.Get(hname)
    # 	  hs.append(h)

    ROOT.gStyle.SetPalette(1)

    canname = 'can'
    can = ROOT.TCanvas(canname, canname, 0, 0, 1000, 800)
    cans.append(can)
    can.Divide(2,2)

    can.cd(1)
    ROOT.gPad.SetLogz(1)
    h2.Draw("colz")
    
    can.cd(2)
    hXmax.Draw('hist')
    
    can.cd(3)
    hsigma.Draw('hist')
    
    can.cd(4)
    hskew.Draw('hist')
    
    # opt = ''
    #for h in hs:
    #   h.Draw(opt + 'hist')
    # 	opt = 'same'

    for can in cans:
    	can.Print(can.GetName() + '.pdf')
    	can.Print(can.GetName() + '.png')
    ROOT.gPad.Update()

    ROOT.gApplication.Run()
    return

###################################
###################################
###################################

if __name__ == "__main__":
    # execute only if run as a script"
    main(sys.argv)
    
###################################
###################################
###################################


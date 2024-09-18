#!/usr/bin/python
# was: #!/snap/bin/pyroot
# was: #!/usr/bin/python3
# Út 17. září 2024, 13:55:18 CEST

import ROOT
from math import sqrt, pow, log, exp
import os, sys, getopt
from scipy.stats import skew
import numpy as np

##########################################
def makeGraph():
    gr = ROOT.TGraphErrors()
    #gr.SetMarkerColor(ROOT.kBlack)
    #gr.SetLineColor(gr.GetMarkerColor())
    gr.SetLineWidth(1)
    gr.SetLineStyle(2)
    gr.SetMarkerSize(1.)
    gr.SetMarkerStyle(20)
    return gr

##########################################
def getMax(hs):
    maxy = -1
    for h in hs:
        val = h.GetMaximum()
        if val > maxy:
            maxy = 1.*val
    return maxy

##########################################
from array import array

# chatgpt:
def custom_palette(NCont = 3):
    NRGBs = 3  # Number of RGB points
    # NCont = n  # Number of colors in the palette

    # Define arrays for the RGB points in the gradient
    stops = array('d', [0.00, 0.50, 1.00])  # Position in the gradient
    red   = array('d', [0.00, 1.00, 1.00])  # Red component (0-1 scale)
    green = array('d', [0.00, 0.00, 1.00])  # Green component (0-1 scale)
    blue  = array('d', [1.00, 0.00, 0.00])  # Blue component (0-1 scale)

    # Create the custom color palette with interpolation
    ROOT.TColor.CreateGradientColorTable(NRGBs-1, stops, green, red, blue, NCont)

    # Set the number of colors for smooth shading
    ROOT.gStyle.SetNumberContours(NCont)

    # The palette is now applied for the session
    print("Custom color palette applied")
    return


##########################################
# https://stackoverflow.com/questions/61521371/calculate-weighted-statistical-moments-in-python
# but FIXED additionally!
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
        #print('avg, var, std: ', w_avg, w_var, w_std)
        #return np.sum(ws * ((values - w_avg)/w_std)**n)/np.sum(ws)
        cn = np.average( (values - w_avg)**n, weights=ws)
        #print(f'c{n}: {cn}')
        return cn

cans = []
stuff = []

##########################################
##########################################
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
        print(f'{argv[0]} simulated_showers/merged/conex_p_E_19_EPOS_merged.root')
        return

    if len(argv) > 2:
        gBatch = int(argv[2])

    print('*** Settings:')
    print('tag={:}, batch={:}'.format(gTag, gBatch))

    if gBatch:
        ROOT.gROOT.SetBatch(1)

    
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

    outdir = 'analyzed/'
    os.system('mkdir -p ' + outdir)

    proto = filename.split('/')[-1].replace('merged', 'analyzed')
    outfname = outdir + proto
    tag = proto.replace('.root','')
    
    print('Output file: ', outfname)
    outfile = ROOT.TFile(outfname, 'recreate')
    outfile.cd()
    
    mx1, mx2 = 0, 1500
    hXmax = ROOT.TH1D('hXmax', ';', 150, mx1, mx2)
    hsigma = ROOT.TH1D('hsigmaXmax', ';#sigma_{Xmax}', 200, 0, 500)
    hskew = ROOT.TH1D('hskewXmax', ';skew of Xmax', 200, -1, 3)

    h1s = []
    grs = []

    # Event loop!
    for ievt in range(0, Nevts):
        #for ievt in range(0, 150):
        if ievt % 100 == 0: 
            print(f'Processing event {ievt} / {Nevts}')
        ientry = tree.LoadTree(ievt)
        nb = tree.GetEntry(ievt)
        
        Xs = []
        dEdXs = []
        h1 = ROOT.TH1D(f'h1_shower_{ievt}', ';x,dE/dx', 100, x1, x2)
        h1.SetLineWidth(1)
        h1.SetLineStyle(2)
        gr = makeGraph()
        gr.SetName(f'gr_shower_{ievt}')
        ip = 0
        # dE/dX is a weight to each position x of the shower development
        for x,dedx,n in zip(tree.X,tree.dEdX,tree.N):
            if n > 0 and dedx > 0:
                h2.Fill(x,dedx)
                h1.Fill(x,dedx)
                Xs.append(x)
                dEdXs.append(dedx)
                err = dedx / sqrt(n)
                gr.SetPoint(ip, x, dedx)
                gr.SetPointError(ip, 0, err)
                ip = ip + 1
        h1s.append(h1)
        grs.append(gr)
        
        Xmax = Xs[dEdXs.index(max(dEdXs))]
        Xaver = np.average(Xs, weights = dEdXs)
        variance = np.average( (Xs - Xaver)**2, weights = dEdXs )
        #print(Xmax, Xaver, variance)
        sigmaXmax = -1.
        skew = -1.
        if variance > 0: 
            sigmaXmax = sqrt(variance)
            skew = get_n_weighted_moment(Xs, dEdXs, 3) / sigmaXmax**3
            #print(sigmaXmax,skew)
            skew2 = np.average(((Xs-Xaver)/np.sqrt(variance))**3, weights=dEdXs)
            #print(skew, skew2)
        hXmax.Fill(Xmax)
        hsigma.Fill(sigmaXmax)
        hskew.Fill(skew)

    #hnames = ['histo_h', '']
    # hs = []
    # for hname in hnames:
    # 	  h = rfile.Get(hname)
    # 	  hs.append(h)

    #ROOT.gStyle.SetPalette(1)
    #ROOT.gStyle.SetNumberContours(1000)
    # Now you can call custom_palette() before drawing histograms
    custom_palette(Nevts)
    
    canname = 'can'
    can = ROOT.TCanvas(canname, canname, 0, 0, 1200, 800)
    cans.append(can)
    can.Divide(2,3)
  
    can.cd(5)
    opt = ''
    ymax = getMax(h1s)
    for h1 in h1s:
        h1.SetMaximum(ymax*1.1)
        h1.Draw('hist plc' + opt)
        opt = 'same'
    
    can.cd(6)
    tmp2 = h2.Clone(h2.GetName() + 'gr_tmp')
    tmp2.Reset()
    stuff.append(tmp2)
    opt = ''
    tmp2.SetStats()
    tmp2.Draw('')

    # https://en.wikipedia.org/wiki/Gaisser%E2%80%93Hillas_function
    form = '(x > [1])*[0]*( (x-[1])/([2]-[1]) )^( ([2]-[1])/([3]) ) * exp(([2] - x)/[3])'
    pns = 'Nmax', 'X0', 'Xmax', 'lambda'

    funs = []
    for gr,h1 in zip(grs,h1s):
        gr.Draw('pmc plc' + opt)
        fname = f'fit{h1s.index(h1)}'
        fun = ROOT.TF1(fname, form, x1, x2)
        fun.SetLineWidth(1)
        fun.SetLineStyle(2)
        pvs = [h1.GetMaximum()/8., x1, h1.GetMean(), 150]
        #h1.Fit(fname, '', '0')
        for pn,pv in zip(pns,pvs):
            ip = pns.index(pn)
            fun.SetParameter(ip, pv)
            fun.SetParName(ip, pn)
        #funcp = fun.DrawCopy('same')
        #funcp.SetLineColor(ROOT.kGreen+2)
        #stuff.append(funcp)
        gr.Fit(fname, '0Q', '')
        fun.Draw('same')
        funs.append(fun)
        opt = ''


    can.cd(2)
    hXmax.Draw('hist')
    
    can.cd(3)
    hsigma.Draw('hist')
    
    can.cd(4)
    hskew.Draw('hist')

    can.cd(1)
    ROOT.gPad.SetLogz(1)
    h2.Draw("colz")

    pngdir, pdfdir = 'png/', 'pdf/'
    os.system('mkdir -p ' + pdfdir + ' ' + pngdir)
    
    for can in cans:
    	can.Print(pdfdir + can.GetName() + '_' + tag + '.pdf')
    	can.Print(pngdir + can.GetName() + '_' + tag + '.png')
    ROOT.gPad.Update()

    for gr in grs:
        gr.Write()
    outfile.Write()
    if not gBatch:
        ROOT.gApplication.Run()
    outfile.Close()
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


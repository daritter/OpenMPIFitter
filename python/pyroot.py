#!/usr/bin/python
#coding: utf8

import sys
import os
import ROOT
import __builtin__
import exceptions
import numpy as np

__builtin__.root = ROOT
__builtin__.gROOT = root.gROOT
__CANSIZE = 1000
__CANVAS = 1
__HISTNR = 1

palette = 0

def default_style():
    global palette
    stops = np.array([0.00, 0.11, 0.12, 0.34, 0.35, 0.38, 0.64, 0.65, 0.66, 0.89, 0.91, 1.00],np.double)
    red   = np.array([0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.098039, 0.932954, 0.970904, 0.996205, 0.999109, 0.909982, 0.500000],np.double)
    green = np.array([0.000000, 0.000000, 0.017647, 0.864706, 0.896078, 1.000000, 1.000000, 0.959332, 0.930283, 0.073348, 0.000726, 0.000000],np.double)
    blue = np.array([0.500000, 0.999109, 1.000000, 0.996205, 0.970904, 0.869703, 0.034788, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],np.double)
    palette = ROOT.TColor.CreateGradientColorTable(len(stops),stops,red,green,blue,256)
    ROOT.gStyle.SetNumberContours(256)
    ROOT.gStyle.SetPadGridX(True)
    ROOT.gStyle.SetPadGridY(True)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

default_style()

def get_canvas(name=None, nx=1, ny=1, title=None, ratio=0.66666, size=__CANSIZE):
    global __CANVAS
    if name is None:
        name = "canvas%d" % __CANVAS
        __CANVAS+=1

    if title is None: title = name

    size = size, int(round(ratio*size))

    canvas = gROOT.FindObject(name)
    if(canvas):
        canvas.SetTitle(title)
        canvas.SetWindowSize(*size)
        canvas.Clear()
    else:
        canvas = root.TCanvas(name,title,*size)

    if nx>1 or ny>1:
        canvas.Clear()
        canvas.Divide(nx,ny,1e-6,1e-6)
        canvas.cd(1)

    return canvas

def get_subpad(canvas, x, y, w, h, name=None):
    if name is None:
        name = canvas.GetName() + "_pad%02d" % canvas.GetListOfPrimitives().GetSize()
    s = root.gStyle
    p = root.TPad(name,name,x,y,x+w,y+h)
    p.SetTopMargin(s.GetPadTopMargin()/h)
    p.SetRightMargin(s.GetPadRightMargin()/w)
    p.SetBottomMargin(s.GetPadBottomMargin()/h)
    p.SetLeftMargin(s.GetPadLeftMargin()/w)
    return p


def draw(tree,expr,range="",cut="",option="", name=None, index=None):
    if isinstance(cut,list)   or isinstance(cut,tuple):   cut   = "&&".join(cut)
    if isinstance(range,list) or isinstance(range,tuple): range = "(%s)" % ",".join(map(str,range))
    if index is not None and index>0: option+=" same"
    if name is None:
        global __HISTNR
        name = "hist%02d" % __HISTNR
        __HISTNR += 1

    tree.Draw("%s >> %s %s" % (expr, name, range),cut,option)
    h =  root.gDirectory.Get(name)
    return h

def close_all():
    gROOT.CloseFiles()

class tree(object):
    cuts = 0
    def __init__(self,filename,treename=None,reusefile=True):
        """
        Create wrapper to TTree object

        The name should be of the form filename:treename. The root file
        will be opened if not already open and searched for the given tree.
        An exception will be thrown if either the file cannot be opened or no
        such tree exists
        """

        if treename is None:
            filename,treename = filename.split(":",1)

        filename = os.path.abspath(filename)
        self.rfile = None
        if reusefile:
            fList = gROOT.GetListOfFiles()
            for i in range(fList.GetSize()):
                if fList[i].GetName() == filename:
                    self.rfile = fList[i]
                    break
        if not self.rfile:
            self.rfile = root.TFile(filename)
        if not self.rfile.IsOpen(): raise exceptions.IOError("Could not open rootfile '%s'" % filename)
        self.rtree = self.rfile.Get(treename)
        if not self.rtree: raise Exception("No tree '%s' in %s" %(treename,filename))

        self.treename = treename

    def draw(self, *args, **argk):
        return draw(self.rtree,*args,**argk)

    def add_cut(self, cut):
        tree.cuts += 1
        name = ">>evts_%s_%02d" % (self.treename,tree.cuts)
        self.rtree.Draw(name,cut)
        self.rtree.SetEventList(root.gDirectory.Get(name[2:]))

    def clear_cuts(self):
        self.rtree.SetEventList(None)

class chain(object):
    def __init__(self, filenames, treename):
        self.chain = root.TChain(treename)
        for f in filenames:
            self.chain.Add(f)

    def draw(self, *args, **argk):
        return draw(self.chain,*args,**argk)

def scale(hist, factor=None, range=None):
    hist.Sumw2()
    if factor is None:
        factor = 1./ hist.GetEffectiveEntries()
    hist.Scale(factor)
    if range is not None:
        hist.GetYaxis().SetRangeUser(*range)

def next_pad():
    """Got to next Adjacent or Child Subpad if possible"""
    n = root.gPad.GetNumber()
    if not n:
        root.gPad.cd(1)
    else:
        root.gPad.GetMother().cd(n+1)
        if(root.gPad.GetNumber() == n):
            root.gPad.GetMother().cd(1)

def histRatio(hist1,hist2):
    """Compute Ratio for 2 histograms"""
    xaxis = hist1.GetXaxis()
    result = root.TH1F(hist1.GetName()+"_"+hist2.GetName(),hist1.GetTitle()+" / "+hist2.GetTitle(),
                       hist1.GetNbinsX(),xaxis.GetXmin(),xaxis.GetXmax())
    hist1.Sumw2()
    hist2.Sumw2()
    result.Divide(hist1,hist2)
    return result

if __name__=="__main__":
    argv=["--colors=LightBG","-noconfirm_exit",'-pprint']
    try:
        __IPYTHON__
    except NameError:
        try:
            from IPython.frontend.terminal.embed import InteractiveShellEmbed
            ipshell = InteractiveShellEmbed(banner1="pyroot (ROOT Version %s)" % root.gROOT.GetVersion())
        except NameError:
            from IPython.Shell import IPShellEmbed
            ipshell = IPShellEmbed(argv=["--colors=LightBG","-noconfirm_exit",'-pprint'],
                                   banner="pyroot (ROOT Version %s)" % root.gROOT.GetVersion())

        ipshell()

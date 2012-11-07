#coding: utf8

"""
Convert ROOT stuff to matplotlib plots

due to the extremly ugly plot possibilities of ROOT, this script was written to
use matplotlib to do the actual plotting of ROOT objects.
"""

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as pl
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.cm import jet
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate

def plot_manager(nx=1,ny=1,*args,**argk):
    """
    Generator to return axes to plot on

    This generator will return a new axes for each call to the next()
    function. It will create new figures as needed. The tiling in
    x and y can be specified with nx and ny, all other parameters
    will be passed to the pyplot.figure() call

    This example will loop over a data collection and plot each collection on a
    single subplot, 4 plots per figure:

    >>> plots = plot_manager(2,2)
    >>> for data in collection:
    >>>     plots.next()
    >>>     pyplot.plot(data)
    """
    i = 0
    fig = None
    while True:
        i+=1
        if i>nx*ny or fig==None:
            fig = pl.figure(*args,**argk)
            i=1
        yield fig.add_subplot(nx,ny,i)

def save_all(basename,pdf=True,png=True,single_pdf=False,close=True):
    """Save all figures"""
    if pdf: pp = PdfPages(basename+".pdf")
    for i in pl.get_fignums():
        fig = pl.figure(i)
        if pdf: pp.savefig(fig)
        if png:
            fig.patch.set_alpha(0.0)
            fig.savefig(basename+"-%02d.png" % i)
        if single_pdf:
            fig.savefig(basename+"-%02d.pdf" % i)

    if pdf: pp.close()
    if close: pl.close("all")

def plot(obj, errors=False, fill=None, **argk):
    """Overloaded plot function working for TH1,TH2 and TF1"""
    if isinstance(obj,root.TH2):
        return plot2D(obj,**argk)
    if isinstance(obj,root.TH1):
        if errors: return plotError(obj,**argk)
        if 'smoothing' in argk: return plotSmooth(obj,**argk)
        return plot1D(obj,fill=fill, **argk)
    if isinstance(obj,root.TF1):
        return plotF1(obj,**argk)


def plot1D(hist, axes = None,  fill=None, **argk):
    """Plot 1D root Histogram"""
    if axes is None: axes = pl.gca()

    nbins = hist.GetNbinsX()
    x = np.zeros(2*nbins+2)
    y = np.zeros(2*nbins+2)

    for i in range(nbins):
        x[2*i+1] = hist.GetBinLowEdge(i+1)
        x[2*i+2] = hist.GetBinLowEdge(i+2)
        y[2*i+1] = hist.GetBinContent(i+1)
        y[2*i+2] = hist.GetBinContent(i+1)

    x[0] = x[1]
    x[-1] = x[-2]
    ymin = min(y)
    y[0] = ymin/2.
    y[-1] = ymin/2.

    if axes.get_autoscalex_on():
        axes.set_xlim(x[0],x[-1])
        axes.set_autoscalex_on(False)

    if fill is not None:
        axes.fill_between(x,0,y,color=fill,alpha=0.05)

    return axes.plot(x,y,**argk)

def plotSmooth(hist, axes=None, smoothing=0, samples=1024, **argk):
    if axes is None: axes = pl.gca()
    nbins = hist.GetNbinsX()
    x = np.zeros(nbins)
    y = np.zeros(nbins)
    w = np.zeros(nbins)

    for i in range(nbins):
        x[i] = hist.GetBinCenter(i+1)
        y[i] = hist.GetBinContent(i+1)
        try:
            w[i] = 1./hist.GetBinError(i+1)
        except:
            w[i] = 1e6

    xmin = hist.GetXaxis().GetXmin()
    xmax = hist.GetXaxis().GetXmax()

    spline = interpolate.splrep(x,y,w,k=3,s=smoothing)
    x = np.linspace(xmin,xmax,samples)
    y = interpolate.splev(x,spline)

    if axes.get_autoscalex_on():
        axes.set_xlim(xmin,xmax)
        axes.set_autoscalex_on(False)
    return axes.plot(x,y,**argk)


def plotError(hist, axes = None, sparse=True, **argk):
    """Plot 1D ROOT Histogram with errorbars"""
    if axes is None: axes = pl.gca()
    #hist.Sumw2()
    nbins = hist.GetNbinsX()
    data = np.zeros((nbins,4))
    points = 0
    for bin in range(1,nbins+1):
        if sparse and not hist.GetBinContent(bin): continue
        data[points] = (
            hist.GetBinCenter(bin),
            hist.GetBinContent(bin),
            hist.GetBinError(bin),
            hist.GetBinWidth(bin)/2.0,
        )
        points+=1

    if points==0: return None

    argk.setdefault("fmt",",")
    argk.setdefault("linewidth",0.5)
    argk.setdefault("capsize",1.0)
    argk.setdefault("marker","None")

    xmin = hist.GetXaxis().GetXmin()
    xmax = hist.GetXaxis().GetXmax()
    if axes.get_autoscalex_on():
        axes.set_xlim(xmin,xmax)
        axes.set_autoscalex_on(False)
    return axes.errorbar(*data[:points].T,**argk)

def getHistData(hist):
    nbinsX = hist.GetNbinsX()
    nbinsY = hist.GetNbinsY()
    data = np.zeros((nbinsX,nbinsY))
    for i in range(nbinsX):
        for j in range(nbinsY):
            data[i,j] = hist.GetBinContent(i+1,j+1)

    xmin = hist.GetXaxis().GetXmin()
    xmax = hist.GetXaxis().GetXmax()
    ymin = hist.GetYaxis().GetXmin()
    ymax = hist.GetYaxis().GetXmax()

    return data,(xmin,xmax,ymin,ymax)


def plotRGBChannel(hists,log=False,axes=None):
    if axes is None: axes = pl.gca()

    if len(hists)>3:
        print "Only three channels avaiable"
        return

    data = []
    extent = None
    shape = None
    for h in hists:
        d,e = getHistData(h)
        if extent is None:
            extent=e
        elif e!=extent:
            "Ranges do not match"
            return
        if shape is None:
            shape = d.T.shape
        elif d.T.shape != shape:
            "Shapes do not match"
            return
        data.append(d)

    col = np.zeros(shape+(3,))+1

    for i,d in enumerate(data):
        d = d.T
        channel = np.zeros(shape+(3,))+1
        idx = np.nonzero(d>0)
        if log:
            cnorm = mpl.colors.LogNorm()
        else:
            cnorm = mpl.colors.Normalize()

        cnorm.autoscale(d[idx])
        for j in range(3):
            factor = i==j and 0.25 or 1.0
            channel[idx + (np.zeros(len(idx[0]),'i')+j,)] = 1-factor*cnorm(d[idx])

        col *= channel

    axes.imshow(col,extent=extent,interpolation="nearest",aspect="auto",origin="lower")

def plot2D(hist, axes=None, sparse=True, contour=None, color=None, log=False, **argk):
    """axes 2D ROOT Histogram as colorfield"""
    if axes is None: axes=pl.gca()

    data,(xmin,xmax,ymin,ymax) = getHistData(hist)

    if isinstance(axes,Axes3D):
        x = np.linspace(xmin,xmax,nbinsX)
        y = np.linspace(ymin,ymax,nbinsY)
        y,x = np.meshgrid(y,x)
        argk.setdefault("cstride",1)
        argk.setdefault("rstride",1)
        argk.setdefault("cmap",jet)
        return axes.plot_surface(x,y,data, **argk)

    if contour is not None:
        return axes.contour(data.T, extent=[xmin,xmax,ymin,ymax], origin='lower', **argk)

    if color is not None:
        cconv = mpl.colors.ColorConverter()
        cnorm = mpl.colors.Normalize()
        cnorm.autoscale(data)
        col = np.array( [[cconv.to_rgba(color,alpha=1.0)]*nbinsX]*nbinsY )
        if sparse:
            col[:,:,3] = cnorm(data.T)*0.8 + 0.2
            col[np.nonzero(data.T==0)] = (1,1,1,0)
        else:
            col[:,:,3] = cnorm(data.T)
        data = col.T

    if log:
        argk["norm"] = mpl.colors.LogNorm()

    argk.setdefault("interpolation","nearest")
    argk.setdefault("aspect","auto")

    if sparse is True:
        data= np.ma.masked_equal(data,0)
    elif sparse is not None:
        data= np.ma.masked_less_equal(data,sparse)

    return axes.imshow(data.T, extent=[xmin,xmax,ymin,ymax], origin='lower', **argk)

def plotF1(func, axes=None, use_range=False, samples=200, scale=1, **argk):
    if axes is None: axes = pl.gca()
    x1,x2 = axes.get_xlim()
    if use_range:
        x1 = max(x1,func.GetXmin())
        x2 = min(x2,func.GetXmax())
    x = np.linspace(x1,x2,samples)
    y = np.array([func(xi)*scale for xi in x])
    axes.set_autoscalex_on(False)
    return axes.plot(x,y,**argk)

if __name__ == '__main__':
    import pyroot as pr
    import sys

    h = root.TH1F("g","",100,-10,10)
    g = root.TF1("g","gaus")
    for x in np.random.normal(0,2,1e6):
        h.Fill(x)
    h.Draw()
    h.Fit(g)

    for x in np.random.normal(5,0.5,1e5):
        h.Fill(x)

    plotError(h)
    plotSmooth(h,smoothing=100)
    plotF1(g)

    pl.show()

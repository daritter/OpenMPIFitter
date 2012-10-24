#!/usr/bin/env python
import sys
import os
sys.path.insert(0,os.path.expanduser("~/belle/DspDsmKs/python"))
import matplotlib
from matplotlib import pyplot as pl
matplotlib.rc("path", simplify = False)


#command handling here
filename = sys.argv[1]
sys.argv = sys.argv[:1]

import pyroot as pr
import r2mpl

rootfile = root.TFile(filename + ".root")

def rescale(data,*fits):
    nData = data.GetNbinsX()
    nFit = fits[0].GetNbinsX()
    if(nData != nFit):
        scale = nFit/nData
        if isinstance(fits[0],root.TH2):
            for f in fits: f.Scale(scale**2)
        else:
            for f in fits: f.Scale(scale)


def calc_sigmas(name, data, fit):
    sigmas = data.Clone(name)
    nData = data.GetNbinsX()
    nFit = fit.GetNbinsX()
    if(nData != nFit):
        scale = nFit/nData
        fit = fit.Clone(fit.GetName()+"_resampled")
        if isinstance(fit,root.TH2):
            fit.Rebin2D(scale,scale)
            fit.Scale(1./scale**2)
        else:
            fit.Rebin(scale)
            fit.Scale(1./scale)

    for i in xrange(sigmas.GetSize()):
        d = data.GetBinContent(i)
        f = fit.GetBinContent(i)
        e = data.GetBinError(i)
        if d!=0:
            sigmas.SetBinContent(i, (d-f)/e)
        sigmas.SetBinError(i,1)
    return sigmas

def plotMbcDe(data, label=None, title=None, **argk):
    f = pl.figure(figsize=(8,8))
    a = f.add_subplot(111)
    p = r2mpl.plot(data, a, **argk)
    if label is None:
        label = "Entries / $%s \\times %s GeV^2$" % (data.GetXaxis().GetBinWidth(1), data.GetYaxis().GetBinWidth(1))
    if title is not None:
        a.set_title(title)
    f.colorbar(p).set_label(label)
    a.grid()
    a.set_xlabel("$M_{BC}$ / GeV")
    a.set_ylabel("$\Delta E$ / GeV")

def plot_dfs(name,data,fits,data_axes,sigma_axes,label,title=None,log=False, unit="GeV"):
    colors = {"signal":"r", "mixed":"g", "charged":"b", "all":"b", "signal, +1":"r", "signal, -1":"b"}
    last = None
    for i,(l,f) in enumerate(fits):
        if last is not None:
            last.Add(f)
        else:
            last = f.Clone("tmp")

        r2mpl.plotSmooth(last, axes=data_axes, color=colors[l], samples=2000, label=l, zorder=i+1)
        #r2mpl.plot(last, axes=data_axes, color=colors[l], label=l)

    r2mpl.plot(data,errors=True,axes=data_axes, color="k", label="data", linewidth=0.5, capsize=1.0, zorder=10)

    sigma = calc_sigmas(name + "_sigma",data, last)
    r2mpl.plot(sigma,errors=True,axes=sigma_axes, color="k", linewidth=0.5, capsize=1.0)
    data_axes.grid()
#    data_axes.legend(loc="best")
    data_axes.set_ylim(ymin=0)
    data_axes.set_xticklabels([])
    sigma_axes.set_ylim(-4,4)
    sigma_axes.axhline(-2, color="r")
    sigma_axes.axhline( 0, color="k")
    sigma_axes.axhline(+2, color="r")
    sigma_axes.set_yticks([-4,-3,-2,-1,0,1,2,3,4])
    sigma_axes.set_yticklabels(["","",-2,"","","",2,"",""])
    sigma_axes.grid()
    sigma_axes.set_xlabel(label)
    sigma_axes.set_ylabel("normalized\nresiduals")
    data_axes.set_ylabel("Entries / %s %s" % (data.GetBinWidth(1), unit))
    #data_axes.get_yaxis().set_label_coords(-0.2,0.5)
    #sigma_axes.get_yaxis().set_label_coords(-0.2,0.5)
    loc = matplotlib.ticker.MaxNLocator(nbins=5)
    sigma_axes.get_xaxis().set_major_locator(loc)
    data_axes.get_xaxis().set_major_locator(loc)
    data_axes.set_title(title)

    if log:
        data_axes.set_yscale("log",nonposy="clip")
        ymin=data.GetMinimum(0.5)/2.0
        for (l,f) in fits:
            ymin = min(ymin,f.GetMaximum()/2.0)
        #data_axes.set_ylim(ymin=data.GetMinimum(0.5)/2.0)
        data_axes.set_ylim(ymin)



def make_mbcde_plots(name, title):
    print "make mBC dE plots for ", name, title
    mbcde_data = rootfile.Get(name + "_data")
    mbcde_fit  = rootfile.Get(name + "_fit")
    #fits = [("all",mbcde_fit)]
    fits = []
    for component in ["charged","mixed","signal"]:
        fit = rootfile.Get(name + "_fit_" + component)
        if(fit):
            fits.append((component,fit))

    if not mbcde_data or not mbcde_fit:
        print "Could not find M_BC:dE histograms"
        sys.exit(1)

    scale = mbcde_fit.GetNbinsX()/mbcde_data.GetNbinsX()
    mbcde_fit2 = mbcde_fit.Clone(name+"_fit2")
    mbcde_fit2.Rebin2D(scale,scale)
    rescale(mbcde_data,mbcde_fit2)
    mbcde_sigma = calc_sigmas(name+"_sigma",mbcde_data, mbcde_fit2)

    vmax = max(mbcde_data.GetMaximum(),mbcde_fit2.GetMaximum())

    plotMbcDe(mbcde_data, vmax=vmax, title=title)
    plotMbcDe(mbcde_fit2, vmax=vmax, title=title, sparse=0.01)
    plotMbcDe(mbcde_sigma, vmin=-5, vmax=5, label="normalized residuals", title=title)
    #plotMbcDe(mbcde_data, vmax=vmax, vmin=0.1, log=True, title=title)
    #plotMbcDe(mbcde_fit2, vmax=vmax, vmin=0.1, log=True, sparse=0.1, title=title)

    mbc_data = mbcde_data.ProjectionX()
    de_data = mbcde_data.ProjectionY()
    mbc_fits = []
    de_fits = []
    for (component,fit) in fits:
        mbc_fit = fit.ProjectionX()
        de_fit = fit.ProjectionY()
        rescale(mbc_data,mbc_fit,de_fit)
        mbc_fits.append((component,mbc_fit))
        de_fits.append((component,de_fit))


    f = pl.figure(figsize=(8,4))
    a_Mbc  = f.add_axes((0.13,0.38,0.35,0.55))
    a_sMbc = f.add_axes((0.13,0.13,0.35,0.25))
    a_dE   = f.add_axes((0.60,0.38,0.35,0.55))
    a_sdE  = f.add_axes((0.60,0.13,0.35,0.25))
    plot_dfs("mbc", mbc_data,mbc_fits,a_Mbc,a_sMbc,"$M_{BC}$ / GeV", title=title)
    plot_dfs("de", de_data,de_fits,a_dE,a_sdE,"$\Delta E$ / GeV", title=title)

    f = pl.figure(figsize=(8,4))
    a_Mbc  = f.add_axes((0.13,0.38,0.35,0.55))
    a_sMbc = f.add_axes((0.13,0.13,0.35,0.25))
    a_dE   = f.add_axes((0.60,0.38,0.35,0.55))
    a_sdE  = f.add_axes((0.60,0.13,0.35,0.25))
    plot_dfs("mbc", mbc_data,mbc_fits,a_Mbc,a_sMbc,"$M_{BC}$ / GeV", title=title, log=True)
    plot_dfs("de", de_data,de_fits,a_dE,a_sdE,"$\Delta E$ / GeV", title=title, log=True)

def make_dt_plots(name, title):
    print "Make dT plot for", name, title
    dt_data_p = rootfile.Get(name + "_data_p")
    dt_data_m = rootfile.Get(name + "_data_m")

    fits_p = []
    fits_m = []
    for component in ["charged","mixed","signal"]:
        fit_p = rootfile.Get(name + "_fit_p_" + component)
        if(fit_p):
            rescale(dt_data_p, fit_p)
            fits_p.append((component + ", +1",fit_p))
        fit_m = rootfile.Get(name + "_fit_m_" + component)
        if(fit_m):
            rescale(dt_data_m, fit_m)
            fits_m.append((component + ", -1",fit_m))

    f = pl.figure(figsize=(8,4))
    a_Mbc  = f.add_axes((0.13,0.38,0.35,0.55))
    a_sMbc = f.add_axes((0.13,0.13,0.35,0.25))
    a_dE   = f.add_axes((0.60,0.38,0.35,0.55))
    a_sdE  = f.add_axes((0.60,0.13,0.35,0.25))

    plot_dfs("dt_p", dt_data_p, fits_p, a_Mbc,a_sMbc,"$\Delta t$ / ps", title=title, unit="ps")
    plot_dfs("dt_m", dt_data_m, fits_m, a_dE,a_sdE,"$\Delta t$ / ps", title=title, unit="ps")

    f = pl.figure(figsize=(8,4))
    a_Mbc  = f.add_axes((0.13,0.38,0.35,0.55))
    a_sMbc = f.add_axes((0.13,0.13,0.35,0.25))
    a_dE   = f.add_axes((0.60,0.38,0.35,0.55))
    a_sdE  = f.add_axes((0.60,0.13,0.35,0.25))

    plot_dfs("dt_p", dt_data_p, fits_p, a_Mbc,a_sMbc,"$\Delta t$ / ps", title=title, unit="ps", log=True)
    plot_dfs("dt_m", dt_data_m, fits_m, a_dE,a_sdE,"$\Delta t$ / ps", title=title, unit="ps", log=True)



make_mbcde_plots("mbcde_svd1", "SVD 1")
make_mbcde_plots("mbcde_svd2", "SVD 2")
make_mbcde_plots("mbcde", "Both")

make_dt_plots("dT_svd1", "SVD1")
make_dt_plots("dT_svd2", "SVD2")

r2mpl.save_all(filename, png=False)
#pl.show()

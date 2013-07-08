#!/usr/bin/env python
import sys
import os
import utils
import matplotlib
from matplotlib import pyplot as pl

#command handling here
filename = sys.argv[1]
sys.argv = sys.argv[:1]

import pyroot as pr
import r2mpl

rootfile = root.TFile(filename + ".root")

def rescale(data,*fits):
    """Rescale histogram according to difference in binsize. Used for oversampling of fit function"""
    scaleX = fits[0].GetNbinsX() / data.GetNbinsX()
    scaleY = 1
    if isinstance(fits[0], root.TH2):
        scaleY = fits[0].GetNbinsY() / data.GetNbinsY()

    for f in fits: f.Scale(scaleX*scaleY)

def accumulate(fits):
    """Accumulate all given histograms into one new histogram"""
    last = None
    for l,f in fits:
        if last is None:
            last = f.Clone("tmp")
        else:
            last.Add(f)
    return last

def calc_sigmas(name, data, fit):
    """Calculate normalized residuals between data and fit histogram"""
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

def plot_mBCdE(data, label=None, title=None, **argk):
    """Make a 2D Plot of mBC and dE"""
    f, a = utils.get_plotaxes()
    #f = pl.figure(figsize=(4,4))
    #a = f.add_axes((0.20,0.13,0.70,0.80))
    p = r2mpl.plot(data, a, **argk)
    if label is None:
        label = r"Entries / $\num{%.3g} \times \num{%.3g}\text{ GeV}^2$" % \
                (data.GetXaxis().GetBinWidth(1), data.GetYaxis().GetBinWidth(1))
    if title is not None:
        a.set_title(title)

    formatter = matplotlib.ticker.ScalarFormatter(False,False)
    formatter.set_powerlimits((-3,3))
    cb = f.colorbar(p, format=formatter, fraction=0.2)
    cb.set_label(label)
    a.grid()
    a.set_xlabel("$M_{BC}$ / GeV")
    a.set_ylabel("$\Delta E$ / GeV")

def plot_dfs(name,data,fits,label,title=None,log=False, unit="GeV"):
    """Plot data, fit and normalized residuals"""
    fig, data_axes, sigma_axes = utils.get_plotaxes_stacked()
    colors = {"signal":"r", "misrecon":"m", "bbar":"g", "continuum":"b", "dummy":"r"}
    last = None
    for i,(l,f) in enumerate(fits):
        if last is not None:
            last.Add(f)
        else:
            last = f.Clone("tmp")

        r2mpl.plotSmooth(last, axes=data_axes, color=colors[l], samples=2000, label=l, zorder=i+1)
        #r2mpl.plot(last, axes=data_axes, color=colors[l], label=l, linewidth=0.2)

        #plot resampled histogram for better comparison
        #scale = last.GetNbinsX()/data.GetNbinsX()
        #if(scale!=1):
            #last2 = last.Clone("tmp_resampled")
            #last2.Rebin(scale)
            #last2.Scale(1./scale)
            #r2mpl.plot(last2, axes=data_axes, color=colors[l], label=l, linewidth=0.2)

    r2mpl.plot(data,errors=True,axes=data_axes, color="k", label="data", linewidth=0.5, capsize=1.0, zorder=10)

    sigma = calc_sigmas(name + "_sigma",data, last)
    r2mpl.plot(sigma,errors=True,axes=sigma_axes, color="k", linewidth=0.5, capsize=1.0)
    data_axes.set_ylim(ymin=0)
    sigma_axes.set_ylim(-4,4)
    sigma_axes.axhline(-2, color="r")
    sigma_axes.axhline( 0, color="k")
    sigma_axes.axhline(+2, color="r")
    sigma_axes.set_yticks([-4,-3,-2,-1,0,1,2,3,4])
    sigma_axes.set_yticklabels(["","","$-2$","","$0$","","$2$","",""])
    sigma_axes.set_xlabel(label)
    sigma_axes.set_ylabel("normalized\nresiduals")
    data_axes.set_ylabel(r"Entries / \num{%.3g} %s" % (data.GetBinWidth(1), unit))
    data_axes.set_title(title)

    if log:
        data_axes.set_yscale("log",nonposy="clip")
        ymin=data.GetMinimum(0.5)/2.0
        for (l,f) in fits:
            ymin = min(ymin,f.GetMaximum()/2.0)
        #data_axes.set_ylim(ymin=data.GetMinimum(0.5)/2.0)
        data_axes.set_ylim(ymin)

def plot_asymmetry(name, data_p, data_m, fits_p, fits_m, label):
    """Plot 2 dT flavours and the asymmetry between them"""
    fig, a1, a2 = utils.get_plotaxes_stacked()

    fp = accumulate(fits_p)
    fm = accumulate(fits_m)
    asym_data = data_p.GetAsymmetry(data_m)
    asym_fit  = fp.GetAsymmetry(fm)

    #r2mpl.plotSmooth(fp, axes=a1, color="r", zorder=1, samples=2000, label="$%s=+1$"%label)
    #r2mpl.plotSmooth(fm, axes=a1, color="b", zorder=1, samples=2000, label="$%s=-1$"%label)
    r2mpl.plot(fp, axes=a1, color="r", zorder=1, label="$%s=+1$"%label, linewidth=0.2)
    r2mpl.plot(fm, axes=a1, color="b", zorder=1, label="$%s=-1$"%label, linewidth=0.2)
    r2mpl.plot(data_p, axes=a1, errors=True, color="r", zorder=2, linewidth=0.5, capsize=1.0)#, marker="^")
    r2mpl.plot(data_m, axes=a1, errors=True, color="b", zorder=2, linewidth=0.5, capsize=1.0)#, marker="o")

    #r2mpl.plotSmooth(asym_fit, axes=a2, zorder=1, color="r", samples=2000)
    r2mpl.plot(asym_fit, axes=a2, zorder=1, color="r", linewidth=0.2)
    r2mpl.plot(asym_data, axes=a2, errors=True, color="k", zorder=2, linewidth=0.5, capsize=1.0)
    a1.legend(loc="upper left", numpoints=2)
    a2.set_ylim(-1.0,1.0)
    a1.set_ylabel(r"Entries / \num{%.3g} %s" % (data_p.GetBinWidth(1), "ps"))
    a2.set_xlabel(r"$\Delta t$ / ps")

def make_mBCdE_plots(name, title):
    """Make all mBC and dE plots for a given data set"""
    print "Make mBC dE plots for ", name
    mbcde_data = rootfile.Get(name + "_data")
    mbcde_fit  = rootfile.Get(name + "_fit")
    fits = []
    for component in ["continuum","bbar","misrecon","signal","dummy"]:
        fit = rootfile.Get(name + "_fit_" + component)
        if(fit):
            fits.append((component,fit))

    if not mbcde_data or not mbcde_fit:
        print "Could not find M_BC:dE histograms for", name
        return

    scale = mbcde_fit.GetNbinsX()/mbcde_data.GetNbinsX()
    mbcde_fit2 = mbcde_fit.Clone(name+"_fit2")
    mbcde_fit2.Rebin2D(scale,scale)
    rescale(mbcde_data,mbcde_fit2)
    mbcde_sigma = calc_sigmas(name+"_sigma",mbcde_data, mbcde_fit2)

    vmax = max(mbcde_data.GetMaximum(),mbcde_fit2.GetMaximum())

    #cmap = matplotlib.cm.get_cmap("binary")
    plot_mBCdE(mbcde_data, vmin=1, vmax=vmax, title=title) #, norm=matplotlib.colors.LogNorm())#, cmap=cmap)
    plot_mBCdE(mbcde_fit2, vmin=1, vmax=vmax, title=title) #, norm=matplotlib.colors.LogNorm())#, cmap=cmap)
    rescmap = matplotlib.cm.get_cmap("bwr")
    plot_mBCdE(mbcde_sigma, vmin=-5, vmax=5, label="normalized residuals", title=title, cmap=rescmap)

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

    plot_dfs("mbc", mbc_data, mbc_fits, "$M_{BC}$ / GeV",   title=title)
    #plot_dfs("mbc", mbc_data, mbc_fits, "$M_{BC}$ / GeV",   title=title, log=True)
    plot_dfs("de",  de_data,  de_fits,  "$\Delta E$ / GeV", title=title)
    #plot_dfs("de",  de_data,  de_fits,  "$\Delta E$ / GeV", title=title, log=True)

def make_dT_plots(name, title):
    """Make all deltaT plots for a given data set"""
    print "Make dT plot for", title
    data = rootfile.Get(name + "_data")

    fits = []
    for component in ["continuum","bbar","misrecon","signal","dummy"]:
        fit = rootfile.Get(name + "_" + component)
        if(fit):
            rescale(data, fit)
            fits.append((component,fit))

    if not (data and fits):
        print "Could not make dT plots for", name
        return 0

    plot_dfs("dt", data, fits, "$\Delta t$ / ps", title=title, unit="ps")#, log=True)
    #plot_dfs("dt_m", dt_data_m, fits_m, "$\Delta t$ / ps", title="%s, $%s=-1$" % (title,label), unit="ps")
    #plot_dfs("dt_p", dt_data_p, fits_p, "$\Delta t$ / ps", title="%s, $%s=+1$" % (title,label), unit="ps", log=True)
    #plot_dfs("dt_m", dt_data_m, fits_m, "$\Delta t$ / ps", title="%s, $%s=-1$" % (title,label), unit="ps", log=True)
    #plot_asymmetry(name, dt_data_p, dt_data_m, fits_p, fits_m, label)
    return 1

try:
    make_mBCdE_plots("mbcde_svd1", "SVD1")
except ValueError:
    pl.close("all")

make_mBCdE_plots("mbcde_svd2", "SVD2")
make_mBCdE_plots("mbcde", "SVD1 + SVD2")
r2mpl.save_all(filename + "-mBCdE", png=False, single_pdf=False)

qname = "pm"
qval = ["-1","+1"]
tname = ["q","qe"]
tval  = ["q", "q\eta"]

for svd in range(2):
    for q in range(2):
        for eta in range(2):
            make_dT_plots("dT_svd{0}_q{1}_e{2}".format(svd+1, qname[q], qname[eta]),
                          "SVD{0}, $q={1}$, $\eta={2}$".format(svd+1, qval[q], qval[eta]))
    for qt in range(2):
        for val in range(2):
            make_dT_plots("dT_svd{0}_{1}{2}".format(svd+1, tname[qt], qname[val]),
                          "SVD{0}, ${1}={2}$".format(svd+1, tval[qt], qval[val]))

    r2mpl.save_all(filename + "-dT-svd{0}".format(svd+1), png=False, single_pdf=False)

# rbin plots
for svd in range(2):
    for rbin in range(7):
        for q in range(2):
            for eta in range(2):
                make_dT_plots("dT_svd{0}_rbin{1}_q{2}_e{3}".format(svd+1, rbin, qname[q], qname[eta]),
                              "SVD{0}, rbin{1}, $q={2}$, $\eta={3}$".format(svd+1, rbin, qval[q], qval[eta]))

    r2mpl.save_all(filename + "-dT-svd{0}-rbins".format(svd+1), png=False, single_pdf=False)

#for t,l in [("q","q"),("qe","q\eta")]:
    #nplots = 0
    #nplots += make_dT_plots("dT_svd1_%s_rbin%d" % (t,7), "SVD1", l)
    #nplots += make_dT_plots("dT_svd2_%s_rbin%d" % (t,7), "SVD2", l)
    #for i in range(7):
        #nplots += make_dT_plots("dT_svd1_%s_rbin%d" % (t,i), "SVD1, rbin %d" % i, l)
    #for i in range(7):
        #nplots += make_dT_plots("dT_svd2_%s_rbin%d" % (t,i), "SVD2, rbin %d" % i, l)

    #if nplots>0:
        #r2mpl.save_all(filename + "-dT-"+t, png=False, single_pdf=False)

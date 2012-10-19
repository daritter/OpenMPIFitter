import sys
import os
import subprocess
sys.path.insert(0,os.path.expanduser("~/belle/DspDsmKs/python"))
import matplotlib
from matplotlib import pyplot as pl


#command handling here

filename = "plots.root"
parameters = "ddk-out.par"
events = "~/belle/DspDsmKs/skim/ddk-signal-correct.root"

subprocess.call(["./ddk-plotter","-i",parameters,"-o",filename,events,"--bins=100","--sampling=4"])#,"--mindE=-0.01","--maxdE=0.01"])

sys.argv = sys.argv[:1]
import pyroot as pr
import r2mpl

rootfile = root.TFile(filename)
mbcde_data = rootfile.Get("mbcde_data")
mbcde_fit  = rootfile.Get("mbcde_fit")
#mbcde_fit.Scale(16)
if not mbcde_data or not mbcde_fit:
    print "Could not find M_BC:dE histograms"
    sys.exit(1)

def calc_sigmas(name, data, fit):
    sigmas = data.Clone(name)
    nData = data.GetNbinsX()
    nFit = fit.GetNbinsX()
    if(nData != nFit):
        scale = nFit/nData
        oldfit = fit
        #fit.Scale(scale)
        fit = fit.Clone(fit.GetName()+"_resampled")
        if isinstance(fit,root.TH2):
            fit.Rebin2D(scale,scale)
            oldfit.Scale(scale**2)
        else:
            oldfit.Scale(scale)
            fit.Rebin(scale)

    for i in xrange(sigmas.GetSize()):
        d = data.GetBinContent(i)
        f = fit.GetBinContent(i)
        e = data.GetBinError(i)
        if d!=0:
            sigmas.SetBinContent(i, (d-f)/e)
        sigmas.SetBinError(i,1)
    return sigmas

scale = mbcde_fit.GetNbinsX()/mbcde_data.GetNbinsX()
mbcde_fit2 = mbcde_fit.Clone("mbcde_fit2")
mbcde_fit2.Rebin2D(scale,scale)
mbcde_sigma = calc_sigmas("mbcde_sigma",mbcde_data, mbcde_fit2)

#pars = matplotlib.figure.SubplotParams(0.1,0.10,0.95,0.95,0.1,0.1)
#plmanager = r2mpl.plot_manager(1,1,figsize=(8,8))#,subplotpars=pars)

def plotMbcDe(data, label=None, **argk):
    f = pl.figure(figsize=(8,8))
    a = f.add_subplot(111)
    p = r2mpl.plot(data, a, **argk)
    if label is None:
        label = "Entries / $%s \\times %s GeV^2$" % (data.GetXaxis().GetBinWidth(1), data.GetYaxis().GetBinWidth(1))
    f.colorbar(p).set_label(label)
    a.grid()
    a.set_xlabel("$M_{BC}$ / GeV")
    a.set_ylabel("$\Delta E$ / GeV")

vmax = max(mbcde_data.GetMaximum(),mbcde_fit2.GetMaximum())

plotMbcDe(mbcde_data, vmax=vmax)
plotMbcDe(mbcde_fit2, vmax=vmax)
plotMbcDe(mbcde_sigma, vmin=-5, vmax=5, label="normalized residuals")

plotMbcDe(mbcde_data, vmax=vmax, vmin=0.1, log=True)
plotMbcDe(mbcde_fit2, vmax=vmax, vmin=0.1, log=True, sparse=0.1)

#f = pl.figure(figsize=(8,8))
#a = f.add_subplot(111)
#p = r2mpl.plot(mbcde_fit2, a, vmin=0, vmax=maxV)
#f.colorbar(p)
#f = pl.figure(figsize=(8,8))
#a = f.add_subplot(111)
#p = r2mpl.plot(mbcde_sigma, a, vmin=-5, vmax=5)
#f.colorbar(p)
#f.tight_layout()

def plot_dfs(name,data,fit,data_axes,sigma_axes,label, log=False):
    sigma = calc_sigmas(name + "_sigma",data,fit)
    r2mpl.plotSmooth(fit, axes=data_axes, color="b", samples=2000)
    #r2mpl.plot(fit, axes=data_axes, color="g")
    r2mpl.plot(data,errors=True,axes=data_axes, color="k")
    r2mpl.plot(sigma,errors=True,axes=sigma_axes, color="b")
    data_axes.grid()
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
    data_axes.set_ylabel("Entries / %s GeV" % data.GetBinWidth(1))
    #data_axes.get_yaxis().set_label_coords(-0.2,0.5)
    #sigma_axes.get_yaxis().set_label_coords(-0.2,0.5)
    loc = matplotlib.ticker.MaxNLocator(nbins=5)
    sigma_axes.get_xaxis().set_major_locator(loc)
    data_axes.get_xaxis().set_major_locator(loc)

    if log:
        data_axes.set_yscale("log",nonposy="clip")
        data_axes.set_ylim(ymin=data.GetMinimum(0.5)/2.0)

f = pl.figure(figsize=(8,4))
a_Mbc  = f.add_axes((0.13,0.33,0.35,0.62))
a_sMbc = f.add_axes((0.13,0.13,0.35,0.20))
a_dE   = f.add_axes((0.60,0.33,0.35,0.62))
a_sdE  = f.add_axes((0.60,0.13,0.35,0.20))

mbc_data = mbcde_data.ProjectionX("mbc_data")
mbc_fit = mbcde_fit.ProjectionX("mbc_fit")
plot_dfs("mbc", mbc_data,mbc_fit,a_Mbc,a_sMbc,"$M_{BC}$ / GeV")

de_data = mbcde_data.ProjectionY("de_data")
de_fit = mbcde_fit.ProjectionY("de_fit")
plot_dfs("de", de_data,de_fit,a_dE,a_sdE,"$\Delta E$ / GeV")


f = pl.figure(figsize=(8,4))
a_Mbc  = f.add_axes((0.13,0.38,0.35,0.57))
a_sMbc = f.add_axes((0.13,0.13,0.35,0.25))
a_dE   = f.add_axes((0.60,0.38,0.35,0.57))
a_sdE  = f.add_axes((0.60,0.13,0.35,0.25))

mbc_data = mbcde_data.ProjectionX("mbc_data")
mbc_fit = mbcde_fit.ProjectionX("mbc_fit")
plot_dfs("mbc", mbc_data,mbc_fit,a_Mbc,a_sMbc,"$M_{BC}$ / GeV", log=True)

de_data = mbcde_data.ProjectionY("de_data")
de_fit = mbcde_fit.ProjectionY("de_fit")
plot_dfs("de", de_data,de_fit,a_dE,a_sdE,"$\Delta E$ / GeV", log=True)


r2mpl.save_all("plots", png=False)
#pl.show()

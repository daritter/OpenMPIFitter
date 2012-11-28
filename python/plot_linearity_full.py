import os
import sys
import utils
import dspdsmks
import pyroot as pr
import r2mpl
import numpy as np
import matplotlib
import toymc
import linearity_full

infile = os.path.splitext(sys.argv[1])[0]
rfile = root.TFile(infile+".root")
colors = {
    "yield_signal_br":"r",
    "signal_dt_Jc":"r",
    "signal_dt_Js1":"g",
    "signal_dt_Js2":"b",
    "signal_dt_blifetime":"k"
}

for par, title in zip(linearity_full.lintest_params, linearity_full.lintest_pnames):
    result = rfile.Get("%s" % par)
    pull = rfile.Get("%s_pull" % par)

    if isinstance(result,root.TH2D):
        fig1, a1 = utils.get_plotaxes()
        fig2, a2 = utils.get_plotaxes()
        fig3, a3 = utils.get_plotaxes()
        a1.set_ylabel("fit result")
        a2.set_ylabel("mean of pull")
        a3.set_ylabel("sigma of pull")
        for a in a1,a2,a3:
            a.set_title("%s, linearity" % title)
            a.set_xlabel("input " + title)


        fit = root.TF1("gauss","gaus")
        lin_result = root.TH1D(
            "lin_%s" % par, "", result.GetNbinsX(),
            result.GetXaxis().GetXmin(), result.GetXaxis().GetXmax())
        lin_pull_mean = root.TH1D(
            "pull_mean_%s" % par, "", result.GetNbinsX(),
            result.GetXaxis().GetXmin(), result.GetXaxis().GetXmax())
        lin_pull_sigma = root.TH1D(
            "pull_sigma_%s" % par, "", result.GetNbinsX(),
            result.GetXaxis().GetXmin(), result.GetXaxis().GetXmax())

        for x in range(1,pull.GetNbinsX()+1):
            one_result = result.ProjectionY(result.GetName()+"_%d" % x,x,x)
            one_pull = pull.ProjectionY(pull.GetName()+"_%d" % x,x,x)

            one_result.Fit(fit,"Q")
            lin_result.SetBinContent(x,fit.GetParameter(1))
            lin_result.SetBinError(x,fit.GetParError(1))
            one_pull.Fit(fit,"Q")
            lin_pull_mean.SetBinContent(x,fit.GetParameter(1))
            lin_pull_mean.SetBinError(x,fit.GetParError(1))
            lin_pull_sigma.SetBinContent(x,fit.GetParameter(2))
            lin_pull_sigma.SetBinError(x,fit.GetParError(2))

        fit = root.TF1("line","pol1")
        for a,h in (a1,lin_result),(a2,lin_pull_mean),(a3,lin_pull_sigma):
            h.Fit(fit,"Q")
            r2mpl.plot(h, axes=a, errors=True, color=colors[par])
            r2mpl.plot(fit, axes=a, color=colors[par])

        ymin = result.GetYaxis().GetXmin()
        ymax = result.GetYaxis().GetXmax()
        a1.set_ylim(ymin,ymax)
        a2.set_ylim(-1,1)
        a3.set_ylim(0,2)

        lin_res = lin_result.Clone("tmp")
        for i in range(1,lin_res.GetNbinsX()+1):
            lin_res.SetBinContent(i,lin_res.GetBinContent(i)-lin_res.GetBinCenter(i))
        r2mpl.plot(lin_res, axes=a1, errors=True, color="k")

    elif isinstance(result,root.TH1D):
        toymc.draw_toyMC(result,title)
        toymc.draw_toyMC(pull,title + ", pull")

r2mpl.save_all(infile,png=False)

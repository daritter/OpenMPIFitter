import utils
import dspdsmks
import pyroot as pr
import r2mpl
import numpy as np
import matplotlib

infile = "toymc-lintest"

lin_params = ["Jc", "Js1", "Js2"]
cpv_params = ["signal_dt_Jc", "signal_dt_Js1", "signal_dt_Js2", "signal_dt_blifetime"]
cpv_names  = [r"$J_C/J_0$", r"$(2J_{s1}/J_0) \sin(2\phi_1)$", r"$(2J_{s2}/J_0) \cos(2\phi_1)$",r"$\tau$"]

rfile = root.TFile(infile+".root")

colors = {
    "signal_dt_Jc":"r",
    "signal_dt_Js1":"g",
    "signal_dt_Js2":"b",
    "signal_dt_blifetime":"k"
}

for lp, lp_name in zip(lin_params,cpv_names):
    fig1, a1 = utils.get_plotaxes()
    fig2, a2 = utils.get_plotaxes()
    fig3, a3 = utils.get_plotaxes()
    a1.set_ylabel("fit result")
    a2.set_ylabel("mean of pull")
    a3.set_ylabel("sigma of pull")
    for a in a1,a2,a3:
        a.set_title("%s, linearity" % lp_name)
        a.set_xlabel("input " + lp_name)

    for cpv, cpv_name in zip(cpv_params, cpv_names):
        result = rfile.Get("%s_%s" % (lp,cpv))
        pull = rfile.Get("%s_%s_pull" % (lp,cpv))

        fit = root.TF1("gauss","gaus")

        lin_result = root.TH1D(
            "lin_%s_%s" % (lp,cpv), "", result.GetNbinsX(),
            result.GetXaxis().GetXmin(), result.GetXaxis().GetXmax())
        lin_pull_mean = root.TH1D(
            "pull_mean_%s_%s" % (lp,cpv), "", result.GetNbinsX(),
            result.GetXaxis().GetXmin(), result.GetXaxis().GetXmax())
        lin_pull_sigma = root.TH1D(
            "pull_sigma_%s_%s" % (lp,cpv), "", result.GetNbinsX(),
            result.GetXaxis().GetXmin(), result.GetXaxis().GetXmax())

        for x in range(1,pull.GetNbinsX()+1):
            one_result = result.ProjectionY(result.GetName()+"_%d" % x,x,x)
            one_pull = pull.ProjectionY(pull.GetName()+"_%d" % x,x,x)

            one_result.Fit(fit,"LQ")
            lin_result.SetBinContent(x,fit.GetParameter(1))
            lin_result.SetBinError(x,fit.GetParError(1))
            one_pull.Fit(fit,"LQ")
            lin_pull_mean.SetBinContent(x,fit.GetParameter(1))
            lin_pull_mean.SetBinError(x,fit.GetParError(1))
            lin_pull_sigma.SetBinContent(x,fit.GetParameter(2))
            lin_pull_sigma.SetBinError(x,fit.GetParError(2))

        fit = root.TF1("line","pol1")
        r2mpl.plot(lin_result, axes=a1, errors=True, color=colors[cpv])#, label=cpv_name)
        lin_result.Fit(fit,"Q")
        r2mpl.plot(fit, axes=a1, color=colors[cpv], label=cpv_name)

        r2mpl.plot(lin_pull_mean, axes=a2, errors=True, color=colors[cpv])#, label=cpv_name)
        lin_pull_mean.Fit(fit,"Q")
        r2mpl.plot(fit, axes=a2, color=colors[cpv], label=cpv_name)

        r2mpl.plot(lin_pull_sigma, axes=a3, errors=True, color=colors[cpv])#, label=cpv_name)
        lin_pull_sigma.Fit(fit,"Q")
        r2mpl.plot(fit, axes=a3, color=colors[cpv], label=cpv_name)


    a1.set_ylim(-1.5,1.5)
    a2.set_ylim(-1.0,1.0)
    a3.set_ylim( 0.0,2.0)
    prop = matplotlib.font_manager.FontProperties(size=8)
    a1.legend(prop=prop,loc="upper left", numpoints=2)
    a2.legend(prop=prop,loc="upper left", numpoints=2)
    a3.legend(prop=prop,loc="upper left", numpoints=2)

r2mpl.save_all(infile,png=False)

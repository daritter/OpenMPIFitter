import os
import toymc
import dspdsmks
import r2mpl
import numpy as np
import utils
import subprocess

nexp = 500
infile = "ddk-out.par"
outfile = "toymc-lintest"
paramfile = "toymc/lintest/params.par"
templates = {
    "signal": "~/belle/DspDsmKs/skim/ddk-signal-correct.root",
    "misrecon": "~/belle/DspDsmKs/skim/ddk-signal-misrecon.root",
    "mixed": "~/belle/DspDsmKs/skim/ddk-mixed.root",
    "charged": "~/belle/DspDsmKs/skim/ddk-charged.root",
}
data = "~/belle/DspDsmKs/skim/ddk-on_resonance.root"
genflags = ["--fudge=2", "--gsim"]
fitflags = ["--fix=.*","--release=yield_signal_.*|signal_dt_J.*|signal_dt_blifetime"]

cpv_pars = ["Jc","Js1","Js2"]

rfile = root.TFile(outfile + ".root","RECREATE")
histograms = {}

def get_dir(par,val):
    directory = "toymc/lintest/_%s/_%+.1f_" % (par,val)
    paramfile = os.path.join(directory,"params.par")
    return directory, paramfile

def get_par(filename):
    dummy, par, val, dummy = [e.strip("/") for e in pfile.split("_")]
    return par, float(val)

cpv_params = ["signal_dt_Jc", "signal_dt_Js1", "signal_dt_Js2", "signal_dt_blifetime"]
cpv_names = [r"$J_C/J_0$", r"$(2J_{s1}/J_0) \sin(2\phi_1)$", r"$(2J_{s2}/J_0) \cos(2\phi_1)$",r"$\tau$"]
params = dspdsmks.Parameters()
for par in cpv_pars:
    histograms[par] = [
        root.TH2D(par+"_"+cpv_params[0],cpv_params[0],21,-1.05,1.05,50,-1.5,1.5),
        root.TH2D(par+"_"+cpv_params[1],cpv_params[1],21,-1.05,1.05,50,-1.5,1.5),
        root.TH2D(par+"_"+cpv_params[2],cpv_params[2],21,-1.05,1.05,50,-1.5,1.5),
        root.TH2D(par+"_"+cpv_params[3],cpv_params[3],21,-1.05,1.05,50,-0.5,0.5),
    ]
    histograms[par+"_pull"] = [
        root.TH2D(par+"_"+name+"_pull",name+"_pull",21,-1.05,1.05,50,-5,5) for name in cpv_params
    ]


    params.load(infile)
    params("yield_signal_svd1").value = utils.nevents[0]
    params("yield_signal_svd2").value = utils.nevents[1]
    params("yield_mixed_svd1").value /= 10
    params("yield_mixed_svd2").value /= 10

    for val in np.linspace(-1,1,21):
        directory, paramfile = get_dir(par,val)
        subprocess.call(["mkdir","-p",directory])
        params("signal_dt_" + par).value = val
        params.save(paramfile)

for i in range(nexp):
    for par in cpv_pars:
        for val in np.linspace(-1,1,21):
            directory, paramfile = get_dir(par,val)
            output = os.path.join(directory,"%03d" % i)
            toymc.add_job(output, paramfile, genflags, fitflags, data, templates)

results = toymc.run_jobs()

for pfile in results:
    params.load(pfile)
    par, val = get_par(pfile)
    pindex = cpv_params.index("signal_dt_"+par)
    cpv_input = [0]*4
    cpv_input[pindex] = val
    for i,(param,initial) in enumerate(zip(cpv_params,cpv_input)):
        p = params(param)
        histograms[par][i].Fill(val, p.value)
        histograms[par+"_pull"][i].Fill(val, (p.value-initial)/p.error)

rfile.Write()
rfile.Close()

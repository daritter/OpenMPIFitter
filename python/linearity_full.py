import sys
import os
import toymc
import dspdsmks
import r2mpl
import numpy as np
import utils
import subprocess
import random

lintest_params = ["signal_dt_Jc", "signal_dt_Js1", "signal_dt_Js2", "yield_signal_br"]#, "ratio_continuum_svd1", "ratio_continuum_svd2"]
lintest_pnames = [ r"$J_C/J_0$", r"$(2J_{s1}/J_0) \sin(2\phi_1)$", r"$(2J_{s2}/J_0) \cos(2\phi_1)$", r"$\mathcal{B}(\ddk)$"]#, r"$f_{q\bar{q}}$, SVD1", r"$f_{q\bar{q}}$, SVD2"]
lintest_irange = [(-1.0, 1.0)]*3 + [(2e-3, 6e-3)]#, (0, 0.3), (0, 0.3)]
lintest_orange = [(-2.5, 2.5)]*3 + [(0, utils.br*2)]#, (0, 0.5), (0, 0.5)]

def get_random_cp():
    while True:
        Jc = random.uniform(-1,1)
        Js1 = random.uniform(-1,1)
        Js2 = random.uniform(-1,1)
        if (Jc**2+Js1**2+Js2**2)**.5<=1:
            return Jc, Js1, Js2

def get_dir(i):
    subdir = i / 500
    directory = "%s/%03d" % (outdir, subdir)
    return directory

def get_par(filename):
    dummy, par, val, dummy = [e.strip("/") for e in pfile.split("_")]
    return par, float(val)

if __name__ == "__main__":
    mctype = sys.argv[1]
    nexp = int(sys.argv[2])
    offset = int(sys.argv[3])

    infile = "ddk-out.par"
    outfile = "toymc/lintest-full-%s-%05d" % (mctype,offset)
    outdir  = "toymc/lintest-full-%s" % mctype
    templates = {
        "signal": "~/belle/DspDsmKs/skim/ddk-signal-correct.root",
        "misrecon": "~/belle/DspDsmKs/skim/ddk-signal-misrecon.root",
        "bbar": "~/belle/DspDsmKs/skim/ddk-bbar.root",
        "continuum": None,
    }
    data = "~/belle/DspDsmKs/skim/ddk-on_resonance.root"
    genflags = ["--fudge=2"]
    if mctype.find("pdf")<0:
        genflags += ["--gsim"]
    if mctype.find("full")>=0:
        genflags += ["--fullgsim"]

    fitflags1 = ["--fix=.*", "--release=yield_signal_.*"]#|signal_dt_blifetime"]
    fitflags2 = ["--fix=.*", "--release=signal_dt_J.*"]#|signal_dt_blifetime"]
    if "bbar" in templates:
        fitflags1[1]+="|yield_bbar_.*"
    #if "continuum" in templates:
    #    fitflags[1]+="|ratio_continuum.*"

    rfile = root.TFile(outfile + ".root", "RECREATE")
    histograms = {}

    for par, irange, orange, title in zip(lintest_params, lintest_irange, lintest_orange, lintest_pnames):
        if irange is not None:
            hist_range = (40, ) + irange + (401, ) + orange
            histograms[par] = root.TH2D(par, title, *hist_range)
            histograms[par+"_pull"] = root.TH2D(par+"_pull", title+" pull", 40, irange[0], irange[1], 101, -5, 5)
        else:
            hist_range = (51, ) + orange
            histograms[par] = root.TH1D(par, title, *hist_range)
            histograms[par+"_pull"] = root.TH1D(par+"_pull", title+" pull", 100, -5, 5)

    params = dspdsmks.Parameters()
    params.load(infile)
    params("signal_svd1_nbb").value, params("signal_svd1_nbb").error = utils.nbb[0]
    params("signal_svd2_nbb").value, params("signal_svd2_nbb").error = utils.nbb[1]
    params("scale_bbar").value = 1
    params("scale_continuum").value = 1

    for i in range(nexp):
        exp_no = i+offset
        directory = get_dir(exp_no)
        output = os.path.join(directory, "%05d" % (exp_no))
        paramfile = output + ".in"
        if not os.path.exists(paramfile):
            #for name, irange in zip(lintest_params, lintest_irange):
            #    if irange is None: continue
            #    params(name).value = random.uniform(*irange)
            Jc,Js1,Js2 = get_random_cp()
            br = random.uniform(*lintest_irange[3])
            #f1, f2 = [random.uniform(*e) for e in lintest_irange[-2:]]
            for name, value in zip(lintest_params, [Jc,Js1,Js2,br]):#,f1,f2]):
                params(name).value = value

            subprocess.call(["mkdir", "-p", directory])
            params.save(paramfile)

        toymc.add_job(output, paramfile, genflags, fitflags1, fitflags2, data, templates)

    results = toymc.run_jobs(False, True)

    pinput = dspdsmks.Parameters()

    for pfile in results:
        params.load(pfile)
        pinput.load(os.path.splitext(pfile)[0] + ".in")
        for par in lintest_params:
            p = params(par)
            initial = pinput(par)
            if isinstance(histograms[par], root.TH2D):
                histograms[par].Fill(initial.value, p.value)
                histograms[par+"_pull"].Fill(initial.value, (p.value-initial.value)/p.error)
            else:
                histograms[par].Fill(p.value)
                histograms[par+"_pull"].Fill((p.value-initial.value)/p.error)

    rfile.Write()
    rfile.Close()

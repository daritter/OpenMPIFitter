import sys
import os
import time
import numpy as np
import dspdsmks
import subprocess
import pyroot as pr
import threading
import multiprocessing
import Queue as queue
import utils
import r2mpl

lock = threading.Lock()

def draw_toyMC(hist,title,xlabel="",ylabel="",exponent=None):
    fit = root.TF1("gauss","gaus")
    bbox_props = dict(boxstyle="round", fc="w", ec="k")
    textbox = r"\begin{align*}\mu&=%s\\\sigma&=%s\end{align*}"
    fig, a = utils.get_plotaxes()
    a.set_title(title)
    a.set_xlabel(xlabel)
    a.set_xlabel(ylabel)
    hist.Fit(fit,"Q")
    r2mpl.plot(hist,axes=a, errors=True, color="k", zorder=1)
    r2mpl.plot(fit,axes=a, color="r", zorder=0)
    a.annotate(textbox % (
        utils.format_error(fit.GetParameter(1), fit.GetParError(1), exponent=exponent),
        utils.format_error(fit.GetParameter(2), fit.GetParError(2), exponent=exponent)),
        xy=(1,1), xycoords="axes fraction", xytext=(-8,-8),
        textcoords="offset points", color="k", ha="right", va="top", bbox=bbox_props)
    return a, (fit.GetParameter(1), fit.GetParError(1)), (fit.GetParameter(2), fit.GetParError(2))


def generate_experiment(basename, params, components, flags, log, deltaT=True):
    did_something = False
    files = []
    deltaT = deltaT and ",deltat" or ""
    for component, events in components.items():
        rootfile = basename + "-%s.root" % component
        files.append(rootfile)
        if os.path.exists(rootfile): continue
        ret = subprocess.call(["./ddk-toymc","--cmp=%s%s" % (component,deltaT), "-i", params, "-o",rootfile, events] + flags, stdout = log)
        if ret!=0:
            if os.path.exists(rootfile): os.remove(rootfile)
            raise Exception("ToyMC generation failed for %s" % basename)
        did_something = True

    return did_something, files

def fit_experiment(basename, initial, files, flags, components, log):
    parfile  = basename + ".par"
    if os.path.exists(parfile): return parfile
    for run in range(10):
        ret = subprocess.call(["./ddk-fitter", "-i", initial, "-o", parfile, "--cmp=%s" % components] + files + flags, stdout = log)
        initial = parfile
        if ret == 0: return parfile

    os.remove(parfile)
    raise Exception("Experiment %s failed to converge" % basename)


def make_experiment(thread_id):
    while True:
        try:
            basename, params, genflags, fitflags, components = jobs.get(True,10)
        except queue.Empty:
            return

        with lock: print "Processing experiment %s on %d" % (basename, thread_id)

        try:
            logfile  = basename + ".log"
            parfile  = basename + ".par"
            log = open(logfile,"w")
            refit, rootfiles = generate_experiment(basename, params, components, genflags, log)
            #if refit and os.path.exists(parfile): os.remove(parfile)
            if os.path.exists(parfile): os.remove(parfile)
            outfile = fit_experiment(basename, params, rootfiles, fitflags, ",".join(components.keys()+["deltat"]), log)
            with lock: results.append(outfile)
        except Exception, e:
            with lock: print e
        finally:
            log.close()

        jobs.task_done()

jobs = queue.Queue()
results = []
def add_job(basename, params, genflags, fitflags, components):
    jobs.put((basename, params, genflags, fitflags, components))

def run_jobs():
    global results
    #Run all jobs
    threads = []
    results = []
    for i in range(multiprocessing.cpu_count()):
        thread = threading.Thread(target=make_experiment, args=(i,))
        thread.daemon=True
        thread.start()
        threads.append(thread)

    while not jobs.empty():
        time.sleep(1)

    jobs.join()
    return results


if __name__ == "__main__":
    nexp = 500
    toyname = "gsim"

    params = dspdsmks.Parameters()
    params.load("ddk-out.par")

    nevents = utils.nevents
    print "Signal events to generate per SVD: ", nevents
    params("yield_signal_svd1").value = nevents[0]
    params("yield_signal_svd2").value = nevents[1]
    params("signal_dt_Jc").value = utils.cpv[0,0]
    params("signal_dt_Js1").value = utils.cpv[1,0]
    params("signal_dt_Js2").value = utils.cpv[2,0]
    params("yield_mixed_svd1").value /= 10
    params("yield_mixed_svd2").value /= 10

    paramfile = "toymc/%s-params.par" % toyname
    params.save(paramfile)
    components = {
        "signal": "~/belle/DspDsmKs/skim/ddk-signal-correct.root",
        "misrecon": "~/belle/DspDsmKs/skim/ddk-signal-misrecon.root",
        "mixed": "~/belle/DspDsmKs/skim/ddk-mixed.root",
        "charged": "~/belle/DspDsmKs/skim/ddk-charged.root",
    }
    genflags = ["--fudge=2"]
    fitflags = ["--fix=.*","--release=yield_signal_.*|signal_dt_J.*"]
    if toyname.find("gsim")>=0:
        #if we generate from pdf take real data as base template
        genflags.append("--gsim")
        for k in components:
            components[k] = "~/belle/DspDsmKs/skim/ddk-on_resonance.root"

    for i in range(nexp):
        add_job("toymc/%s-%03d" % (toyname,i), paramfile, genflags, fitflags, components)

    br_result = root.TH1D("brresult","",100,0,2*utils.br)
    br_pull   = root.TH1D("brpull","",50,-5,5)
    jc_result = root.TH1D("jcresult","",50,-1.5,1.5)
    js1_result = root.TH1D("js1result","",50,-1.5,1.5)
    js2_result = root.TH1D("js2result","",50,-1.5,1.5)
    jc_pull = root.TH1D("jcpull","",50,-5,5)
    js1_pull = root.TH1D("js1pull","",50,-5,5)
    js2_pull = root.TH1D("js2pull","",50,-5,5)

    cpv_result = [jc_result, js1_result, js2_result]
    cpv_pull = [jc_pull, js1_pull, js2_pull]

    for parfile in run_jobs():
        if not os.path.exists(parfile): continue
        params.load(parfile)
        yields = np.array([params("yield_signal_svd1").value,params("yield_signal_svd2").value])
        yield_errors = np.array([params("yield_signal_svd1").error,params("yield_signal_svd2").error])
        result = yields / (utils.efficiency_mc[:,0] * utils.nbb[:,0])
        result_errors = yield_errors / (utils.efficiency_mc[:,0] * utils.nbb[:,0])
        combined = sum(yields / (utils.efficiency_mc[:,0])) / sum(utils.nbb[:,0])
        combined_error = sum((yield_errors / (utils.efficiency_mc[:,0]))**2)**.5 / sum(utils.nbb[:,0])
        br_result.Fill(combined)
        br_pull.Fill((combined-utils.br)/combined_error)

        for i,p in enumerate((params("signal_dt_Jc"), params("signal_dt_Js1"), params("signal_dt_Js2"))):
            cpv_result[i].Fill(p.value)
            cpv_pull[i].Fill((p.value-utils.cpv[i,0])/p.error)


    a,m,s = draw_toyMC(br_result, "Branching ratio, input=%.3g" % utils.br, "Fitted Br", exponent=-3)
    a.set_ylim(0, br_result.GetMaximum()*1.5)
    draw_toyMC(br_pull,"Branching ratio pull")

    names = [r"$J_C/J_0$", r"$(2J_{s1}/J_0) \sin(2\phi_1)$", r"$(2J_{s2}/J_0) \cos(2\phi_1)$"]

    for i in range(3):
        name = names[i]
        hist = cpv_result[i]
        pull = cpv_pull[i]
        draw_toyMC(hist, r"%s, input=%.3g" % (name, utils.cpv[i,0]))
        draw_toyMC(pull, r"Pull for %s" % (name))

    r2mpl.save_all("toymc-%s" % toyname, png=False)

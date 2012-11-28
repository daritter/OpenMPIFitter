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
import matplotlib

lock = threading.Lock()

def draw_toyMC(hist,title,xlabel="",ylabel="",exponent=None):
    fit = root.TF1("gauss","gaus")
    bbox_props = dict(boxstyle="round", fc="w", ec="k")
    font_props = matplotlib.font_manager.FontProperties(size="small")
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
        textcoords="offset points", color="k", ha="right", va="top", bbox=bbox_props, font_properties=font_props)
    return a, (fit.GetParameter(1), fit.GetParError(1)), (fit.GetParameter(2), fit.GetParError(2))


def generate_experiment(basename, params, data, components, flags, log, deltaT=True):
    did_something = False
    files = []
    deltaT = deltaT and ",deltat" or ""
    for component, events in components.items():
        rootfile = basename + "-%s.root" % component
        files.append(rootfile)
        if os.path.exists(rootfile): continue
        ret = subprocess.call(["./ddk-toymc","--cmp=%s%s" % (component,deltaT), "-i", params, "-o",rootfile, "--template=%s" % events, data] + flags, stdout = log)
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


def make_experiment(thread_id, force_refit):
    while True:
        try:
            basename, params, genflags, fitflags, data, components = jobs.get(True,10)
        except queue.Empty:
            return

        with lock: print "Processing experiment %s on %d" % (basename, thread_id)

        try:
            logfile  = basename + ".log"
            parfile  = basename + ".par"
            log = open(logfile,"w")
            refit, rootfiles = generate_experiment(basename, params, data, components, genflags, log)
            if (refit or force_refit) and os.path.exists(parfile): os.remove(parfile)
            outfile = fit_experiment(basename, params, rootfiles, fitflags, ",".join(components.keys()+["deltat"]), log)
            with lock: results.append(outfile)
        except Exception, e:
            with lock: print e
        finally:
            log.close()

        jobs.task_done()

jobs = queue.Queue()
results = []
def add_job(basename, params, genflags, fitflags, data, components):
    jobs.put((basename, params, genflags, fitflags, data, components))

def run_jobs(refit=False):
    global results
    #Run all jobs
    threads = []
    results = []
    for i in range(multiprocessing.cpu_count()-1):
        thread = threading.Thread(target=make_experiment, args=(i,refit))
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
    templates = {
        "signal": "~/belle/DspDsmKs/skim/ddk-signal-correct.root",
        "misrecon": "~/belle/DspDsmKs/skim/ddk-signal-misrecon.root",
        "mixed": "~/belle/DspDsmKs/skim/ddk-mixed.root",
        "charged": "~/belle/DspDsmKs/skim/ddk-charged.root",
    }
    data = "~/belle/DspDsmKs/skim/ddk-on_resonance.root"

    genflags = ["--fudge=2"]
    fitflags = ["--fix=.*","--release=yield_signal_.*|signal_dt_J.*|signal_dt_blifetime"]
    if toyname.find("gsim")<0:
        genflags.append("--gsim")

    for i in range(nexp):
        add_job("toymc/%s-%03d" % (toyname,i), paramfile, genflags, fitflags, data, templates)

    br_result = root.TH1D("brresult","",50,0,2*utils.br)
    br_pull   = root.TH1D("brpull","",50,-5,5)
    jc_result = root.TH1D("jc","",50,-1.5,1.5)
    js1_result = root.TH1D("js1","",50,-1.5,1.5)
    js2_result = root.TH1D("js2","",50,-1.5,1.5)
    bl_result = root.TH1D("blifetime","",50,-0.5,0.5)

    cpv_result = [jc_result, js1_result, js2_result, bl_result]
    cpv_pull = [root.TH1D(e.GetName()+"_pull","",50,-5,5) for e in cpv_result]
    cpv_params = ["signal_dt_Jc", "signal_dt_Js1", "signal_dt_Js2", "signal_dt_blifetime"]
    cpv_input = [params(e).value for e in cpv_params]

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

        for i,p in enumerate([params(e) for e in cpv_params]):
            cpv_result[i].Fill(p.value)
            cpv_pull[i].Fill((p.value-cpv_input[i])/p.error)


    a,m,s = draw_toyMC(br_result, "Branching ratio, input=%.3g" % utils.br, "Fitted Br", exponent=-3)
    a.set_ylim(0, br_result.GetMaximum()*1.5)
    draw_toyMC(br_pull,"Branching ratio pull")

    names = [r"$J_C/J_0$", r"$(2J_{s1}/J_0) \sin(2\phi_1)$", r"$(2J_{s2}/J_0) \cos(2\phi_1)$",r"$\tau$"]

    for i,name in enumerate(names):
        hist = cpv_result[i]
        pull = cpv_pull[i]
        draw_toyMC(hist, r"%s, input=%.3g" % (name, cpv_input[i]))
        draw_toyMC(pull, r"Pull for %s" % (name))

    r2mpl.save_all("toymc-%s" % toyname, png=False)

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

def draw_toyMC(hist,title,xlabel="",ylabel="",exponent=None, box="tl"):
    fit = root.TF1("gauss","gaus")
    textbox = r"\begin{align*}\mu&=%s\\\sigma&=%s\end{align*}"
    fig, a = utils.get_plotaxes((4,3.2))
    hist.Fit(fit,"LQ")
    r2mpl.plot(hist,axes=a, errors=True, color="k", zorder=1)
    r2mpl.plot(fit,axes=a, color="r", zorder=0)
    utils.text_box(a, box, textbox % (
        utils.format_error(fit.GetParameter(1), fit.GetParError(1), exponent=exponent),
        utils.format_error(fit.GetParameter(2), fit.GetParError(2), exponent=exponent)))

    a.set_title(title)
    a.set_xlabel(xlabel)
    a.set_ylabel(ylabel)
    return a, (fit.GetParameter(1), fit.GetParError(1)), (fit.GetParameter(2), fit.GetParError(2))


def generate_experiment(basename, params, data, components, flags, log, deltaT=True):
    did_something = False
    files = []
    deltaT = deltaT and ",deltat" or ""
    for component, events in components.items():
        rootfile = basename + "-%s.root" % component
        files.append(rootfile)
        if os.path.exists(rootfile): continue
        log.write("Generating events for %s\n" % component)
        log.flush()
        arguments = ["./ddk-toymc","--cmp=%s%s" % (component,deltaT), "-i", params, "-o",rootfile]
        if events is not None:
            arguments.append("--template=%s" % events)
        ret = subprocess.call(arguments + [data] + flags, stdout = log)
        if ret!=0:
            if os.path.exists(rootfile): os.remove(rootfile)
            raise Exception("ToyMC generation failed for %s" % basename)
        did_something = True

    return did_something, files

def fit_experiment(basename, initial, files, flags1, flags2, components, log):
    log.write("Fitting experiment with %s\n" % components)
    log.flush()
    parfile  = basename + ".par"
    if os.path.exists(parfile): return parfile
    ret = subprocess.call(["./ddk-fitter", "-i", initial, "-o", parfile, "--cmp=%s" % components] + files + flags1, stdout = log)
    if flags2 is not None:
        ret += subprocess.call(["./ddk-fitter", "-i", parfile, "-o", parfile, "--cmp=%s" % components] + files + flags2, stdout = log)
    if ret == 0: return parfile

    os.remove(parfile)
    raise Exception("Experiment %s failed to converge" % basename)


def make_experiment(thread_id, force_refit, deltaT):
    deltaTcmp = deltaT and ["deltat"] or []
    #if deltaT:
    #    dtc.append["deltat"]

    while True:
        try:
            basename, params, genflags, fitflags1, fitflags2, data, components = jobs.get(True,10)
        except queue.Empty:
            return

        with lock: print "Processing experiment %s on %d" % (basename, thread_id)

        try:
            logfile  = basename + ".log"
            parfile  = basename + ".par"
            log = open(logfile,"w")
            refit, rootfiles = generate_experiment(basename, params, data, components, genflags, log, deltaT)
            if (refit or force_refit) and os.path.exists(parfile): os.remove(parfile)
            outfile = fit_experiment(basename, params, rootfiles, fitflags1, fitflags2, ",".join(components.keys()+deltaTcmp), log)
            with lock: results.append(outfile)
        except Exception, e:
            with lock: print e
        finally:
            log.close()

        jobs.task_done()

jobs = queue.Queue()
results = []
def add_job(basename, params, genflags, fitflags1, fitflags2, data, components):
    jobs.put((basename, params, genflags, fitflags1, fitflags2, data, components))

def run_jobs(refit=False, deltaT=True):
    global results
    #Run all jobs
    threads = []
    results = []
    for i in range(multiprocessing.cpu_count()-1):
        thread = threading.Thread(target=make_experiment, args=(i,refit, deltaT))
        thread.daemon=True
        thread.start()
        threads.append(thread)

    while not jobs.empty():
        time.sleep(1)

    jobs.join()
    return results


if __name__ == "__main__":
    nexp = int(sys.argv[2]) #200
    toyname = sys.argv[1] #"gsim"

    params = dspdsmks.Parameters()
    params.load("ddk-out.par")

    params("signal_svd1_nbb").value, params("signal_svd1_nbb").error = utils.nbb[0] #* 10
    params("signal_svd2_nbb").value, params("signal_svd2_nbb").error = utils.nbb[1] #* 10
    params("signal_dt_Jc").value = utils.cpv[0,0]
    params("signal_dt_Js1").value = utils.cpv[1,0]
    params("signal_dt_Js2").value = utils.cpv[2,0]
    params("scale_bbar").value = 1
    params("scale_continuum").value = 1

    paramfile = "toymc/%s-params.par" % toyname
    params.save(paramfile)
    templates = {
        "signal": "~/belle/DspDsmKs/skim/ddk-signal-correct.root",
        "misrecon": "~/belle/DspDsmKs/skim/ddk-signal-misrecon.root",
        "bbar": "~/belle/DspDsmKs/skim/ddk-bbar.root",
        "continuum": None,
    }
    data = "~/belle/DspDsmKs/skim/ddk-on_resonance.root"

    genflags = ["--fudge=2"]
    fitflags1 = ["--fix=.*","--release=yield_signal_.*"]
    fitflags2 = ["--fix=.*","--release=signal_dt_J.*|signal_dt_blifetime"]
    if "bbar" in templates:
        fitflags1[-1] += "|yield_bbar.*"

    #fitflags = ["--fix=.*","--release=yield_.*|signal_dt_J.*"]

    if toyname.find("gsim")>=0:
        genflags.append("--gsim")
    if toyname.find("full")>=0:
        genflags.append("--fullgsim")
    if toyname.find("select")>=0:
        genflags.append("--select")


    for i in range(nexp):
        add_job("toymc/%s-%04d" % (toyname,i), paramfile, genflags, fitflags1, fitflags2, data, templates)

    br_result = root.TH1D("brresult","",100,0,2*utils.br)
    br_pull   = root.TH1D("brpull","",100,-5,5)
    jc_result = root.TH1D("jc","",100,-1.5,1.5)
    js1_result = root.TH1D("js1","",100,-1.5,1.5)
    js2_result = root.TH1D("js2","",100,-1.5,1.5)
    bl_result = root.TH1D("blifetime","",50,1.0,2.0)
    bbar_svd1_result = root.TH1D("bbar_svd1","",100,params("yield_bbar_svd1").value*0.5,params("yield_bbar_svd1").value*1.5)
    bbar_svd2_result = root.TH1D("bbar_svd2","",100,params("yield_bbar_svd2").value*0.5,params("yield_bbar_svd2").value*1.5)

    cpv_result = [jc_result, js1_result, js2_result, bl_result, bbar_svd1_result, bbar_svd2_result]
    cpv_pull = [root.TH1D(e.GetName()+"_pull","",100,-5,5) for e in cpv_result]
    cpv_params = ["signal_dt_Jc", "signal_dt_Js1", "signal_dt_Js2", "signal_dt_blifetime", "yield_bbar_svd1", "yield_bbar_svd2"]
    cpv_input = [params(e).value for e in cpv_params]

    for parfile in run_jobs(False, True):
        if not os.path.exists(parfile): continue
        params.load(parfile)
        br = params("yield_signal_br").value
        br_error = params("yield_signal_br").error
        br_result.Fill(br)
        br_pull.Fill((br-utils.br)/br_error)

        for i,p in enumerate([params(e) for e in cpv_params]):
            cpv_result[i].Fill(p.value)
            cpv_pull[i].Fill((p.value-cpv_input[i])/p.error)

    a,m,s = draw_toyMC(br_result, r"fit results, input=$%s$" % utils.format_error(utils.br, precision=1), "$\mathcal{B}(\ddk)$", exponent=-3, ylabel=r"Entries / \num{%s}" % br_result.GetBinWidth(1))
    a.set_ylim(0, br_result.GetMaximum()*1.5)
    draw_toyMC(br_pull,"pull distribution", xlabel=r"Pull($\mathcal{B}(\ddk)$)", ylabel="Entries / %s" % br_pull.GetBinWidth(1))

    names = [r"$J_C/J_0$", r"$(2J_{s1}/J_0) \sin(2\phi_1)$", r"$(2J_{s2}/J_0) \cos(2\phi_1)$",r"$\tau / ps$","m1","m2"]

    for i,name in enumerate(names):
        hist = cpv_result[i]
        pull = cpv_pull[i]
        draw_toyMC(hist, r"fit results, input=$%.3g$" % cpv_input[i], xlabel=name, ylabel="Entries / %s" % hist.GetBinWidth(1))
        draw_toyMC(pull, r"pull distribution", xlabel="Pull(%s)" % name, ylabel="Entries / %s" % pull.GetBinWidth(1))

    r2mpl.save_all("toymc/toymc-%s" % toyname, png=False, single_pdf=True)

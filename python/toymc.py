import sys
import os
import time
import numpy as np
import dspdsmks
import subprocess
import pyroot as pr
import threading
import Queue as queue
import matplotlib
matplotlib.use("Agg")
matplotlib.rc("path", simplify=False)
matplotlib.rc("font", family="serif")
matplotlib.rc("text", usetex=True)
matplotlib.rc("text.latex", unicode="true", preamble=r"\usepackage{sistyle},\SIthousandsep{},\usepackage{hepnames},\DeclareRobustCommand{\PDstpm}{\HepParticle{D}{}{\ast\pm}\xspace}")
from matplotlib import pyplot as pl
import utils
import r2mpl

lock = threading.Lock()

def make_experiment(thread_id):
    while True:
        try:
            i = jobs.get(True,10)
        except queue.Empty:
            return

        with lock:
            print "Processing experiment", i, "on", thread_id
        rootfile = "toymc/signal-%03d.root" % i
        parfile  = "toymc/signal-%03d.par" % i
        logfile  = "toymc/signal-%03d.log" % i
        log = open(logfile,"w")
        ret = subprocess.call(["./ddk-toymc","--cmp=signal,deltat","-i","toymc/params.par","-o",rootfile,
                        "~/belle/DspDsmKs/skim/ddk-signal-correct.root"], stdout = log)
        if ret!=0:
           with lock:
               print "ToyMC generation failed for %s" % i

        converged = False
        for run in range(10):
            ret = subprocess.call(["./ddk-fitter","--cmp","signal,deltat",
                                "-i","toymc/params.par","-o",parfile,"--fix=.*","--release=yield_signal.*|signal_dt_J.*",
                                rootfile], stdout = log)
            if ret == 0:
                converged = True
                break

        log.close()

        if not converged:
            with lock:
                print "Experiment %s failed to converge" % i

            os.remove(parfile)

        jobs.task_done()


nexp = 500
br = 4.1e-3

params = dspdsmks.Parameters()
params.load("ddk-out.par")

nevents = utils.nbb[:,0] * utils.efficiency_mc[:,0] * br
params("yield_signal_svd1").value = nevents[0]
params("yield_signal_svd2").value = nevents[1]
params("signal_dt_Jc").value = utils.cpv[0,0]
params("signal_dt_Js1").value = utils.cpv[1,0]
params("signal_dt_Js2").value = utils.cpv[2,0]

params.save("toymc/params.par")

jobs = queue.Queue()
for i in range(nexp):
    jobs.put(i)

threads = []
for i in range(16):
    thread = threading.Thread(target=make_experiment, args=(i,))
    thread.daemon=True
    thread.start()
    threads.append(thread)

while not jobs.empty():
    time.sleep(1)

jobs.join()

br_result = root.TH1D("brresult","",100,0,2*br)
br_pull   = root.TH1D("brpull","",50,-5,5)
jc_result = root.TH1D("jcresult","",50,-1.5,1.5)
js1_result = root.TH1D("js1result","",50,-1.5,1.5)
js2_result = root.TH1D("js2result","",50,-1.5,1.5)
jc_pull = root.TH1D("jcpull","",50,-5,5)
js1_pull = root.TH1D("js1pull","",50,-5,5)
js2_pull = root.TH1D("js2pull","",50,-5,5)

cpv_result = [jc_result, js1_result, js2_result]
cpv_pull = [jc_pull, js1_pull, js2_pull]

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
    return a


for i in range(nexp):
    parfile  = "toymc/signal-%03d.par" % i
    if not os.path.exists(parfile): continue
    params.load(parfile)
    yields = np.array([params("yield_signal_svd1").value,params("yield_signal_svd2").value])
    yield_errors = np.array([params("yield_signal_svd1").error,params("yield_signal_svd2").error])
    result = yields / (utils.efficiency_mc[:,0] * utils.nbb[:,0])
    result_errors = yield_errors / (utils.efficiency_mc[:,0] * utils.nbb[:,0])
    combined = sum(yields / (utils.efficiency_mc[:,0])) / sum(utils.nbb[:,0])
    combined_error = sum((yield_errors / (utils.efficiency_mc[:,0]))**2)**.5 / sum(utils.nbb[:,0])
    #br_result.Fill(result[0])
    #br_result.Fill(result[1])
    #br_pull.Fill((result[0]-br)/result_errors[0])
    #br_pull.Fill((result[1]-br)/result_errors[1])
    br_result.Fill(combined)
    br_pull.Fill((combined-br)/combined_error)

    for i,p in enumerate((params("signal_dt_Jc"), params("signal_dt_Js1"), params("signal_dt_Js2"))):
        cpv_result[i].Fill(p.value)
        cpv_pull[i].Fill((p.value-utils.cpv[i,0])/p.error)


a = draw_toyMC(br_result, "Branching ratio, input=%.3g" % br, "Fitted Br", exponent=-3)
a.set_ylim(0, br_result.GetMaximum()*1.5)
a = draw_toyMC(br_pull,"Branching ratio pull")

names = [r"$J_C/J_0$", r"$(2J_{s1}/J_0) \sin(2\phi_1)$", r"$(2J_{s2}/J_0) \cos(2\phi_1)$"]

for i in range(3):
    name = names[i]
    hist = cpv_result[i]
    pull = cpv_pull[i]
    draw_toyMC(hist, r"%s, input=%.3g" % (name, utils.cpv[i,0]))
    draw_toyMC(pull, r"Pull for %s" % (name))

r2mpl.save_all("toymc", png=False)

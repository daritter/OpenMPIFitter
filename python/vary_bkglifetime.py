import sys
import subprocess
import dspdsmks
import random
import multiprocessing
import os

cfg = "ctrl-run.d/10-box-lifetime.ini"
par = "ctrl-lifetime.par"

vary = [
    "bkg_svd1_dt_mean_tau",
    "bkg_svd1_dt_fraction_delta_sgl",
    "bkg_svd1_dt_fraction_delta_mul",
    "bkg_svd2_dt_mean_tau",
    "bkg_svd2_dt_fraction_delta_sgl",
    "bkg_svd2_dt_fraction_delta_mul",
]

start = int(sys.argv[1])
end = int(sys.argv[2])

params = dspdsmks.Parameters()


for i in range(start, end):
    params.load(par)
    print "Fit", i
    name = "toymc/ctrl-dt/%03d" % i
    if os.path.exists(name + ".out"): continue
    log = open(name + ".log", "w")
    for p in vary:
        params(p).value += random.gauss(0, params(p).error)
    params.save(name + ".in")
    ret = subprocess.call(["mpirun","-np", str(multiprocessing.cpu_count()-1), "./ddk-fitter","-c", cfg, "-i", name+".in", "-o", name+".out"], stdout=log)
    log.close()
    if ret!=0:
        print "failed"
        os.remove(name+".out")
    else:
        print "done"

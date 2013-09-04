import sys
import os
import glob
import dspdsmks
import utils
import matplotlib
from matplotlib import pyplot as pl
import pyroot as pr
import r2mpl
import toymc

params = dspdsmks.Parameters()

lifetime = root.TH1D("b_lifetime", "b_lifetime", 40, 1.20, 1.55)
pull_lifetime = root.TH1D("p_lifetime", "b_lifetime", 40, -3, 3)

for par in glob.glob("toymc/ctrl-dt/*.out"):
    params.load(par)
    v = params("signal_dt_blifetime").value
    e = params("signal_dt_blifetime").error
    lifetime.Fill(v)
    pull_lifetime.Fill((v-1.519)/(e**2+0.007**2)**.5)

toymc.draw_toyMC(lifetime, r"fit result", r"$\tau_{\PBz}$ / ps", ylabel=r"Entries / \num{%s}" % lifetime.GetBinWidth(1), box="tr")[0].axvline(1.519,c="b",ls="--")
toymc.draw_toyMC(pull_lifetime, r"pull distribution", r"Pull($\tau_{\PBz}$)", ylabel=r"Entries / \num{%s}" % pull_lifetime.GetBinWidth(1), box="tr")[0].axvline(0,c="b",ls="--")

r2mpl.save_all("vary_lifetime")

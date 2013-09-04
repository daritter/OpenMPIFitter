#!/usr/bin/env python
import sys
import os
import utils
import matplotlib
from matplotlib import pyplot as pl

#Calculate number three:
#s12 + s13 + s23 = m2B0 + m2D* + m2D0 + m2K-

#command handling here
filename = sys.argv[1]
sys.argv = sys.argv[:1]

import pyroot as pr
import r2mpl
control = True
if control:
    splus, sminus, s23 = ("$m(\PDstp\PKm)^2$ / GeV", "$m(\PDz\PKm)^2$ / GeV", "$m(\PDstp\PDz)^2$ / GeV")
else:
    splus, sminus, s23 = ("$m(\PDstp\PKzS)^2$ / GeV", "$m(\PDstm\PKzS)^2$ / GeV", "$m(\PDstp\PDstm)^2$ / GeV")

rootfile = root.TFile(filename + ".root")

def draw_dalitz(hist, splus, sminus):
    hist.Sumw2()
    hist.Rebin2D(2,2)
    f, a = utils.get_plotaxes()
    p = r2mpl.plot(hist, a)#, aspect="equal")
    label = r"Entries / $\num{%.3g} \times \num{%.3g}\text{ GeV}^2$" % \
        (hist.GetXaxis().GetBinWidth(1), hist.GetYaxis().GetBinWidth(1))

    formatter = matplotlib.ticker.ScalarFormatter(False,False)
    formatter.set_powerlimits((-3,3))
    cb = f.colorbar(p, format=formatter, fraction=0.2)
    cb.set_label(label)
    a.grid(False)
    a.set_xlabel(splus)
    a.set_ylabel(sminus)

    f, a = utils.get_plotaxes()
    h_x = hist.ProjectionX()
    #h_x.Rebin(2)
    r2mpl.plot(h_x, errors=True)
    a.set_xlabel(splus)
    a.set_ylabel(r"Entries / $\num{%.3g}\text{ GeV}$" % h_x.GetXaxis().GetBinWidth(1))

    f, a = utils.get_plotaxes()
    h_y = hist.ProjectionY()
    #h_y.Rebin(2)
    r2mpl.plot(h_y, errors=True)
    a.set_xlabel(sminus)
    a.set_ylabel(r"Entries / $\num{%.3g}\text{ GeV}$" % h_y.GetXaxis().GetBinWidth(1))

hist1 = rootfile.Get("sDalitz_svd2")
hist2 = rootfile.Get("sDalitz2_svd2")
hist3 = rootfile.Get("sDalitz3_svd2")
draw_dalitz(hist1, splus, sminus)
draw_dalitz(hist2, splus, s23)
draw_dalitz(hist3, sminus, s23)

r2mpl.save_all(filename, png=False, single_pdf=False)

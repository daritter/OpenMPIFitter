#!/usr/bin/env python
#coding:utf8
import sys
import os
import numpy as np
import math
import utils

#command handling here
#ngenerated = float(sys.argv[1])
filenames = sys.argv[1:]
sys.argv = sys.argv[:1]
print filenames

import pyroot as pr
import dspdsmks
import r2mpl

br_Dp = [0.0913]
br_D0 = [0.0387, 0.139, 0.0807]
br_DsDp = 0.307
br_DsD0 = 0.677
br_K0Ks = 0.500
br_KsPP = 0.692
br_DDK = 8.1e-3

eff_DsDp = br_DsDp * sum(br_Dp)
eff_DsD0 = br_DsD0 * sum(br_D0)
eff_DsDs = (eff_DsDp + eff_DsD0)**2;
eff_DDKs = eff_DsDs * br_KsPP
eff_DsD0Km = (eff_DsDp + eff_DsD0) * sum(br_D0)

#signal_yield_par = ("yield_svd1_signal", "yield_svd2_signal")
#params = dspdsmks.Parameters(sys.argv[1])
#signal_yield = [params(e) for e in signal_yield_par]

b0 = pr.chain(filenames[:1],"B0")
b0g = pr.chain(filenames[:1],"B0gen")
h_ngenerated = b0g.draw("svdVs", range=(2,0,2), option="goff")
#h_nevents = b0.draw("svdVs",cut="bestLHsig.mcInfo&1 && bestLHsig.flag==0 && bestLHsig.tag.flavour!=0", range=(2,0,2), option="goff")

h_ncorrect = b0.draw("bestLHsig.mcInfo&1", cut="", range=(2,0,2), option="goff")
h_bestB = b0.draw("bestLHsig.mcInfo&1", cut="nB0>1 && trueB0>0", range=(2,0,2), option="goff")
ncorrect = [h_ncorrect.GetBinContent(i+1) for i in range(2)]
bestB = [h_bestB.GetBinContent(i+1) for i in range(2)]

pdf = dspdsmks.DspDsmKsPDF()
pdf.set_filenames(filenames[1:])
pdf.load()

signal_yield = np.array([(pdf.size(i), math.sqrt(pdf.size(i))) for i in range(2)])

h_yield = b0.draw("svdVs", range=(2,0,2), cut=utils.cuts() + " && (bestLHsig.mcInfo&1)==1", option="goff")
s2 = np.array([(h_yield.GetBinContent(i+1), h_yield.GetBinError(i+1)) for i in range(2)])
if not np.array_equal(signal_yield,s2):
    print "utils.cuts seem to be not working"
    print "From fitter:", signal_yield
    print "From cuts:  ", s2

ngenerated = np.array([(h_ngenerated.GetBinContent(i+1), h_ngenerated.GetBinContent(i+1)) for i in range(2)])

raw_eff = signal_yield / ngenerated
rec_eff = raw_eff * eff_DDKs

ctl_eff = raw_eff * eff_DsD0Km

h_nB0 = b0.draw("nB0", range=(200,0,200), option="goff")
b_mult = h_nB0.GetMean()
f,a = utils.get_plotaxes()
r2mpl.plot(h_nB0, axes=a, log=True)
a.set_title("\PBz Multiplicity")
a.set_xlabel("Number of reconstructed \PBz")
a.set_ylabel("Number of events")
r2mpl.save_all("bmult", png=False, single_pdf=False)
print """ngenerated = np.array([
    %d,
    %d,
])

efficiency_mc = np.array([
   [%.7e, %.7e],
   [%.7e, %.7e],
])""" % (tuple(ngenerated[:,0]) + tuple(rec_eff.flat))

print """efficiency_ctrl = np.array([
   [%.7e, %.7e],
   [%.7e, %.7e],
])""" % tuple(ctl_eff.flat)

#print """efficency_mc = [
    #np.array([%.5e, %.5e]),
    #np.array([%.5e, %.5e]),
#]""" % tuple(rec_eff.flat)

print """
signal_svd1_eff                  %17.10e %17.10e      0      0     Y
signal_svd2_eff                  %17.10e %17.10e      0      0     Y
""" % tuple(rec_eff.flat)

print """
#Control channel
signal_svd1_eff                  %17.10e %17.10e      0      0     Y
signal_svd2_eff                  %17.10e %17.10e      0      0     Y
""" % tuple(ctl_eff.flat)

print r"""
\def\BMultiplicity{%.1f}
\def\FractionCorrectRec{\SI{%.1f}{\%%}}
\def\FractionCorrectBestB{\SI{%.1f}{\%%}}
\def\ReconstructedBR{\num{%.3e}}
\def\RawReconstructionEffSVDOne{%s}
\def\RawReconstructionEffSVDTwo{%s}
\def\ReconstructionEffSVDOne{%s}
\def\ReconstructionEffSVDTwo{%s}""" % (
    b_mult,
    100*ncorrect[1] / (ncorrect[0]+ncorrect[1]),
    100*bestB[1] / (bestB[0]+bestB[1]),
    eff_DDKs,
    utils.format_error(*raw_eff[0], exponent=-3),
    utils.format_error(*raw_eff[1], exponent=-3),
    utils.format_error(*rec_eff[0], exponent=-5),
    utils.format_error(*rec_eff[1], exponent=-5),
)

print r"""
\def\BMultiplicityCtrl{%.1f}
\def\FractionCorrectRecCtrl{\SI{%.1f}{\%%}}
\def\FractionCorrectBestBCtrl{\SI{%.1f}{\%%}}
\def\ReconstructedBRCtrl{\num{%.3e}}
\def\RawReconstructionEffSVDOneCtrl{%s}
\def\RawReconstructionEffSVDTwoCtrl{%s}
\def\ReconstructionEffSVDOneCtrl{%s}
\def\ReconstructionEffSVDTwoCtrl{%s}""" % (
    b_mult,
    100*ncorrect[1] / (ncorrect[0]+ncorrect[1]),
    100*bestB[1] / (bestB[0]+bestB[1]),
    eff_DsD0Km,
    utils.format_error(*raw_eff[0], exponent=-3),
    utils.format_error(*raw_eff[1], exponent=-3),
    utils.format_error(*ctl_eff[0], exponent=-5),
    utils.format_error(*ctl_eff[1], exponent=-5),
)

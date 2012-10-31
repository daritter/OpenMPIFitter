#!/usr/bin/env python
#coding:utf8
import sys
import os
import numpy as np
import math
import utils

#command handling here
ngenerated = float(sys.argv[1])
filenames = sys.argv[2:]
sys.argv = sys.argv[:1]

import pyroot as pr

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
#eff_DDK  = eff_DDKs * br_K0Ks

#signal_mc = (6e6,6e6)
#signal_yield_par = ("yield_svd1_signal", "yield_svd2_signal")
#sub_br = eff_DDKs
#params = dspdsmks.Parameters(sys.argv[1])
#signal_yield = [params(e) for e in signal_yield_par]
b0 = pr.chain(filenames,"B0")
h_nevents = b0.draw("svdVs",cut="bestLHsig.mcInfo&1 && bestLHsig.flag==0 && bestLHsig.tag.flavour!=0", range=(2,0,2), option="goff")

signal_yield = np.array([(h_nevents.GetBinContent(i+1), h_nevents.GetBinError(i+1)) for i in range(2)])

print signal_yield

raw_eff = signal_yield / ngenerated
rec_eff = raw_eff * eff_DDKs

print r"""
\def\ReconstructedBR{\num{%.3e}}
\def\RawReconstructionEffSVDOne{%s}
\def\RawReconstructionEffSVDTwo{%s}
\def\ReconstructionEffSVDOne{%s}
\def\ReconstructionEffSVDTwo{%s}""" % (
    eff_DDKs,
    utils.format_error(*raw_eff[0], exponent=-3),
    utils.format_error(*raw_eff[1], exponent=-3),
    utils.format_error(*rec_eff[0], exponent=-5),
    utils.format_error(*rec_eff[1], exponent=-5),
)

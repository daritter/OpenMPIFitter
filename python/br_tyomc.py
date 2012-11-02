#!/usr/bin/env python
#coding:utf8
import sys
import os
import numpy as np
import math
import utils
import random
import progressbar

#command handling here
filenames = sys.argv[1:]
sys.argv = sys.argv[:1]

import pyroot as pr
import r2mpl

channels =          [ "K-p+", "K-p+p0", "K-p+p+p-", "K-p+p+",  "D0p+",  "D+p0",    "Ks"]
br       = np.array([3.87e-2,  13.9e-2,    8.07e-2,  9.13e-2, 67.7e-2, 30.7e-2, 69.2e-2])
br_err   = np.array([0.05e-2,   0.5e-2,     0.2e-2,  0.19e-2,  0.5e-2,  0.5e-2, 0.05e-2])
whichDs  =          [      4,        4,          4,        5]
pdgIndex =          [     31,       49,         66,       40,       1,       2,       2]

correlations = {
    (1,0):  -4e-2,
    (2,0):  22e-2,
    (2,1):  55e-2,
    (5,4): -66e-2,
}

np.set_printoptions(precision=4, linewidth=120)
cor = np.matrix(np.identity(len(br)))
#note_cor = np.matrix(np.identity(len(br)))
#note_rep = {0:2,1:3,2:4,3:5,4:0,5:1,6:6}
for (i,j),v in correlations.items():
    cor[i,j] = v
    cor[j,i] = v
#    i = note_rep[i]
#    j = note_rep[j]
#    note_cor[i,j] = v
#    note_cor[j,i] = v

#print utils.latex_matrix(note_cor)

diag_err = np.matrix(np.diag(br_err))

cov = diag_err * cor *  diag_err

print "Branching Ratio Indicies:"
print channels
print "Branching Ratio Br (from PDG):"
print br
print "Branching Ratio errors δBr (from PDG):"
print br_err
print "Correlation matrix corr(Br) (from PDG):"
print cor
print "Covarianve matrix cov(Br) = diag(δBr)*corr(Br)*diag(δBr):"
print cov

def get_br(br, deriv, index, scale):
    if deriv != index: return br[index], scale
    if scale == 0:
        return 1,1
    else:
        return br[index],scale+1

def calc_eff(br, deriv=-1):
    n=4
    ks=6
    total_br = 0
    for d1 in range(n):
        for d2 in range(n):
            ds1 = whichDs[d1]
            ds2 = whichDs[d2]
            scale = deriv<0 and 1 or 0
            bd1, scale  = get_br(br, deriv, d1,  scale)
            bd2, scale  = get_br(br, deriv, d2,  scale)
            bds1, scale = get_br(br, deriv, ds1, scale)
            bds2, scale = get_br(br, deriv, ds2, scale)
            bks, scale  = get_br(br, deriv, ks,  scale)
            total_br += scale * bd1 * bd2 * bds1 * bds2 * bks
    return total_br

def calc_err(br, cov):
    n = 6
    total_err = 0
    deriv = [calc_eff(br,i) for i in range(n)]
    for i in range(n):
        total_err += deriv[i]**2 * cov[i,i]

    for i in range(n-1):
        for j in range(i+1,n):
            total_err += 2*deriv[i]*deriv[j]*cov[i,j]

    return math.sqrt(total_err)

eff_DsDp = br[4] * sum(br[0:3])
eff_DsD0 = br[5] * sum(br[3:4])
eff_DsDs = (eff_DsDp + eff_DsD0)**2;
eff_DDKs = eff_DsDs * br[6]

total_br = calc_eff(br)
total_err = calc_err(br, diag_err*diag_err)
print "Assuming uncorrelated: BR=%.3e, dBR=%.3e (%.2f%%)" % (total_br, total_err, total_err/total_br*100)

#def br_mod():
    #offset = np.zeros(br_err.shape)
    #for i in range(len(br_err)):
        #offset[i] = random.gauss(0,br_err[i])
    #return offset

#h = root.TH1D("br","br",100,0,0)
#h.SetBuffer(50000)
#for i in range(50000):
    #h.Fill(calc_eff(br + br_mod()))

#h.BufferEmpty()
#print "Uncorrelated, toymc: BR=%.3e, dBR=%.3e (%.2f%%)" % (h.GetMean(), h.GetRMS(), h.GetRMS()/h.GetMean()*100)

total_br = calc_eff(br)
total_err = calc_err(br, cov)
print "PDG correlations:  BR=%.3e, dBR=%.3e (%.2f%%)" % (total_br, total_err, total_err/total_br*100)

total_errp = 0
total_errm = 0
for i in range(len(br_err)):
    for j in range(-1,2,2):
        offset = np.zeros(br_err.shape)
        offset[i] += j*br_err[i]
        err = calc_eff(br+offset) - total_br
        if err<0: total_errm -= err
        else: total_errp += err
print "Conservative:  BR=%.3e, dBR= +%.3e (%.2f%%), -%.3e (%.2f%%)" % (total_br, total_errp, total_errp/total_br*100, total_errm, total_errm/total_br*100)

#h = root.TH1D("brerr","brerr",100,0,0)
#h.SetBuffer(100000)
#pbar = progressbar.ProgressBar(maxval=500000).start()
#for i in range(500000):
    #pbar.update(i)
    #sys.stdout.flush()
    #dcor = np.matrix(np.identity(7))
    #for (i,j) in correlations:
        #v = random.uniform(-1,1)
        #dcor[i,j] = v
        #dcor[j,i] = v

    #dcov = diag_err * dcor *  diag_err
    #err = calc_err(br, dcov)
    #h.Fill(100*err/total_br)
#pbar.finish()
#h.BufferEmpty()

#h.Draw()
#sys.stdin.readline()

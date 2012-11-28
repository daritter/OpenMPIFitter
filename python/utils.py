#coding: utf8
import math
import numpy as np
import matplotlib

matplotlib.use("Agg")
matplotlib.rc("path", simplify=False)
matplotlib.rc("font", family="serif")
matplotlib.rc("text", usetex=True)
matplotlib.rc("legend", fancybox=True, fontsize="small")
matplotlib.rc("text.latex", unicode="true", preamble=r"""
              \usepackage{sistyle},
              \SIthousandsep{},
              \usepackage{hepnames},
              \DeclareRobustCommand{\PDstpm}{\HepParticle{D}{}{\ast\pm}\xspace},
              \DeclareRobustCommand{\PDstp}{\HepParticle{D}{}{\ast+}\xspace},
              \DeclareRobustCommand{\PDstm}{\HepParticle{D}{}{\ast-}\xspace},
              \def\ddk{\ensuremath{\PBz \rightarrow \PDstp\PDstm\PKzS}\xspace}""")

from matplotlib import pyplot as pl

efficiency_mc = np.array([
   [8.15824e-05, 6.22585e-07],
   [3.19346e-04, 1.23177e-06],
])

#efficiency_mc = np.array([
    #[5.63679e-05, 5.17507e-07],
    #[2.17889e-04, 1.01746e-06],
#])

nbb = np.array([
    [151.961, 1.241],
    [619.620, 9.441],
])*1e6

cpv = np.array([
    [0.71,0.16],
    [0.03,0.21],
    [0.24,0.22],
])

br = 4.1e-3

nevents = nbb[:,0] * efficiency_mc[:,0] * br

def get_plotaxes():
    """Return axes"""
    f = pl.figure(figsize=(4,4))
    a = f.add_axes((0.20,0.13,0.73,0.75))
    formatter = matplotlib.ticker.ScalarFormatter(False,False)
    formatter.set_powerlimits((-3,3))
    a.get_xaxis().set_major_formatter(formatter)
    a.get_yaxis().set_major_formatter(formatter)
    a.grid()
    #a.get_yaxis().set_label_coords(-0.15,0.5)
    return f,a

def get_plotaxes_stacked():
    """Return axes suitable for fit + residuals or dT + asymmetry plots"""
    f = pl.figure(figsize=(4,4))
    a1  = f.add_axes((0.20,0.38,0.73,0.55))
    a2 = f.add_axes((0.20,0.13,0.73,0.25))
    for a,p in (a1,None), (a2,"upper"):
        formatter = matplotlib.ticker.ScalarFormatter(False,False)
        formatter.set_powerlimits((-3,3))
        a.get_xaxis().set_major_formatter(formatter)
        a.get_yaxis().set_major_formatter(formatter)
        if p is not None:
            locy = matplotlib.ticker.MaxNLocator(prune=p, nbins=5)
            a.get_yaxis().set_major_locator(locy)
        a.grid()
        a.get_yaxis().set_label_coords(-0.15,0.5)

    a1.set_xticklabels([])
    return f,a1,a2

def format_error(value,error, precision=3, exporder=3, exponent=None):
    if exponent is None: exponent = int(math.floor(math.log10(abs(value))))
    show_exp = False
    if abs(exponent)>=exporder:
        show_exp = True
        value /= math.pow(10,exponent)
        error /= math.pow(10,exponent)

    if(show_exp):
        return r"({0:.{precision}f} \pm {1:.{precision}f}) \times 10^{{{exp}}}".format(
            value, error, exp=exponent, precision=precision)
    else:
        return r"{0:.{precision}f} \pm {1:.{precision}f}".format(
            value, error, precision=precision)

def latex_matrix(matrix):
    lines = []
    matrix = np.array(matrix)
    for row in matrix:
        lines.append(" & ".join(str(e) for e in row))
    return  "\\begin{pmatrix}\n    " + "\\\\\n    ".join(lines) + "\n\\end{pmatrix}"

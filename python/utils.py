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

ngenerated = np.array([
    15000000,
    14999988,
])

efficiency_mc = np.array([
   [5.3093427e-05, 3.1765148e-07],
   [1.9792079e-04, 6.1330483e-07],
])
efficiency_ctrl = np.array([
   [9.7382357e-04, 5.0533889e-06],
   [1.7508169e-03, 6.7758386e-06],
])

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
br_ctrl = 2.47e-3

nevents = nbb[:,0] * efficiency_mc[:,0] * br
nevents_ctrl = nbb[:,0] * efficiency_ctrl[:,0] * br_ctrl

def cuts(min_Mbc=5.24, max_Mbc=5.30, min_dE=-0.15, max_dE=0.1):
    return "{bb}.flag==0 && {bb}.tag.flavour!=0 && \
    ({bb}.Mbc>{min_Mbc} && {bb}.Mbc<{max_Mbc}) &&\
    ({bb}.dE>{min_dE} && {bb}.dE<{max_dE}) &&\
    ({bb}.deltaZ*{inv_bgc}>{min_dT} && {bb}.deltaZ*{inv_bgc}<{max_dT}) && !(\
    ({bb}.vtx.ntrk>1    && {bb}.vtx.chi2/{bb}.vtx.ndf > {qual_cut}) ||\
    ({bb}.tag.ntrk>1    && {bb}.tag.chi2/{bb}.tag.ndf > {qual_cut}) ||\
    ({bb}.vtx.ntrk == 1 && {bb}.vtx.zerr > {sngl_cut}) ||\
    ({bb}.vtx.ntrk > 1  && {bb}.vtx.zerr > {mult_cut}) ||\
    ({bb}.tag.ntrk == 1 && {bb}.tag.zerr > {sngl_cut}) ||\
    ({bb}.tag.ntrk > 1  && {bb}.tag.zerr > {mult_cut}) )".format(
    bb="bestLHsig", qual_cut=50, mult_cut=0.02, sngl_cut=0.05,
    min_Mbc=min_Mbc, max_Mbc=max_Mbc,
    min_dE=min_dE, max_dE=max_dE,
    min_dT=-70, max_dT=70,
    inv_bgc=78.48566945838871754705,
    )

##Cuts to make sure we only use the correct Events
#cuts = "{bb}.flag==0 && {bb}.tag.flavour!=0 && \
    #({bb}.Mbc>{min_Mbc} && {bb}.Mbc<{max_Mbc}) &&\
    #({bb}.dE>{min_dE} && {bb}.dE<{max_dE}) &&\
    #({bb}.deltaZ*{inv_bgc}>{min_dT} && {bb}.deltaZ*{inv_bgc}<{max_dT}) && !(\
    #({bb}.vtx.ntrk>1    && {bb}.vtx.chi2/{bb}.vtx.ndf > {qual_cut}) ||\
    #({bb}.tag.ntrk>1    && {bb}.tag.chi2/{bb}.tag.ndf > {qual_cut}) ||\
    #({bb}.vtx.ntrk == 1 && {bb}.vtx.zerr > {sngl_cut}) ||\
    #({bb}.vtx.ntrk > 1  && {bb}.vtx.zerr > {mult_cut}) ||\
    #({bb}.tag.ntrk == 1 && {bb}.tag.zerr > {sngl_cut}) ||\
    #({bb}.tag.ntrk > 1  && {bb}.tag.zerr > {mult_cut}) )".format(
    #bb="bestLHsig", qual_cut=50, mult_cut=0.02, sngl_cut=0.05,
    #min_Mbc=5.24, max_Mbc=5.30,
    #min_dE=-0.15, max_dE=0.1,
    #min_dT=-70, max_dT=70,
    #inv_bgc=78.48566945838871754705,
#)

def get_plotaxes(figsize=(4,4)):
    """Return axes"""
    f = pl.figure(figsize=figsize)
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

def format_error(value, error=None, precision=3, exporder=3, exponent=None, align=False):
    if exponent is None:
        if value==0:
            exponent = 0
        else:
            exponent = int(math.floor(math.log10(abs(value))))
    show_exp = False
    if abs(exponent)>=exporder:
        show_exp = True
        value /= math.pow(10,exponent)
        if error is not None:
            error /= math.pow(10,exponent)

    operator = align and r'&\pm&' or r'\pm'

    if error is not None:
        if(show_exp):
            return r"({0:.{precision}f} {op} {1:.{precision}f}) \times 10^{{{exp}}}".format(
                value, error, exp=exponent, precision=precision, op=operator)
        else:
            return r"{0:.{precision}f} {op} {1:.{precision}f}".format(
                value, error, precision=precision, op=operator)
    else:
        if(show_exp):
            return r"{0:.{precision}f} \times 10^{{{exp}}}".format(
                value, exp=exponent, precision=precision)
        else:
            return r"{0:.{precision}f}".format(
                value, precision=precision)

def latex_matrix(matrix):
    lines = []
    matrix = np.array(matrix)
    for row in matrix:
        lines.append(" & ".join(str(e) for e in row))
    return  "\\begin{pmatrix}\n    " + "\\\\\n    ".join(lines) + "\n\\end{pmatrix}"

def text_box(axes, pos, text):
    bbox_props = dict(boxstyle="round", fc="w", ec="k")
    font_props = matplotlib.font_manager.FontProperties(size="small")
    if pos == "tr":
        pos_kwargs = {
            "xy":(1,1),
            "xytext":(-8,-8),
            "ha":"right",
            "va":"top",
        }
    elif pos == "tl":
        pos_kwargs = {
            "xy":(0,1),
            "xytext":(8,-8),
            "ha":"left",
            "va":"top",
        }
    elif pos == "br":
        pos_kwargs = {
            "xy":(1,0),
            "xytext":(-8,8),
            "ha":"right",
            "va":"bottom",
        }
    elif pos == "bl":
        pos_kwargs = {
            "xy":(0,0),
            "xytext":(8,8),
            "ha":"left",
            "va":"bottom",
        }

    axes.annotate(text, xycoords="axes fraction", textcoords="offset points",
               color="k", bbox=bbox_props, font_properties=font_props, **pos_kwargs)

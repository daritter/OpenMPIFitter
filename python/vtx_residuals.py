import sys
import pyroot as pr

b0 = pr.chain(sys.argv[1:],"B0")
cuts = "{bb}.flag==0 && {bb}.tag.flavour!=0 && \
    ({bb}.Mbc>{min_Mbc} && {bb}.Mbc<{max_Mbc}) &&\
    ({bb}.dE>{min_dE} && {bb}.dE<{max_dE}) &&\
    ({bb}.deltaZ*{inv_bgc}>{min_dT} && {bb}.deltaZ*{inv_bgc}<{max_dT}) && !(\
    ({bb}.vtx.ntrk>1    && {bb}.vtx.chi2/{bb}.vtx.ndf > {qual_cut}) ||\
    ({bb}.tag.ntrk>1    && {bb}.tag.chi2/{bb}.tag.ndf > {qual_cut}) ||\
    ({bb}.vtx.ntrk == 1 && {bb}.vtx.zerr > {sngl_cut}) ||\
    ({bb}.vtx.ntrk > 1  && {bb}.vtx.zerr > {mult_cut}) ||\
    ({bb}.tag.ntrk == 1 && {bb}.tag.zerr > {sngl_cut}) ||\
    ({bb}.tag.ntrk > 1  && {bb}.tag.zerr > {mult_cut}) )"

cuts = cuts.format(
    bb="bestLHsig", qual_cut=50, mult_cut=0.02, sngl_cut=0.05,
    min_Mbc=5.24, max_Mbc=5.30,
    min_dE=-0.15, max_dE=0.1,
    min_dT=-70, max_dT=70,
    inv_bgc=78.48566945838871754705,
)

b0.draw("(bestLHsig.tag.z - mcVtxTag)*1e4", range=(200,-1e3,1e3), cut=cuts)

sys.stdin.readline()

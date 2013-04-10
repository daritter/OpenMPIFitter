import sys
import utils
import pyroot as pr

b0 = pr.chain(sys.argv[1:],"B0")
b0.draw("(bestLHsig.tag.z - mcVtxTag)*1e4", range=(200,-1e3,1e3), cut=utils.cuts)

sys.stdin.readline()

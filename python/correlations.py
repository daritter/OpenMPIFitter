import sys
import utils
import dspdsmks
import pyroot as pr
import r2mpl
import numpy as np
import matplotlib

#bestB = "bestLHsig"
#components = ["Mbc","dE","deltaZ"]
#filenames = sys.argv[1:]

#data = pr.chain(filenames,"B0")

#corr = np.zeros([len(components)]*2)

#for i, ni in enumerate(components):
#    for j, nj in enumerate(components):
#        h = data.draw("{0}.{2}:{0}.{1}".format(bestB,ni,nj))
#        corr[i,j] = h.GetCorrelationFactor()
#print corr

filenames = sys.argv[1:]

axes = "xyz"

c = pr.get_canvas("projections",3,3)
sub = 1

for filename in sys.argv[1:]:
    rfile = root.TFile(filename)
    data = rfile.Get("mbcdedt_data")
    corr = np.zeros([len(axes)]*2)
    for i, ni in enumerate(axes):
        for j, nj in enumerate(axes):
            c.cd(sub)
            sub+=1
            if i==j:
                h = data.Project3D(ni)
                h.Draw()
                #corr[i,j]=1
            else:
                h = data.Project3D(ni+nj)
                h.Draw("colz")
                #corr[i,j] = h.GetCorrelationFactor()
                #print h.GetCorrelationFactor()

            corr[i,j] = data.GetCorrelationFactor(i+1,j+1)

    corr[1,0] = corr[0,1]

    print corr
    print utils.latex_matrix(corr)
    print r"""
\begin{tabular}{lrrr}
    \toprule
    &\mbc&\de&\dt\\
    \midrule
    \mbc & %.3f & %.3f & %.3f \\
    \de  & %.3f & %.3f & %.3f \\
    \dt  & %.3f & %.3f & %.3f \\
    \bottomrule
\end{tabular}""" % tuple(corr.flat)

sys.stdin.readline()

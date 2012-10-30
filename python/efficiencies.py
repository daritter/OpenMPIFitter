#!/usr/bin/env python
#coding:utf8
import sys
import os
sys.path.insert(0,os.path.expanduser("~/belle/DspDsmKs/python"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rc("path", simplify = False)
from matplotlib import pyplot as pl

#command handling here
filenames = sys.argv[1:]
sys.argv = sys.argv[:1]


import pyroot as pr
import r2mpl

bn835 = np.array(map(float, """
0.0068 ± 0.00026    0.0090 ± 0.00030
0.0014 ± 0.00012    0.0022 ± 0.00015
0.0025 ± 0.00016    0.0035 ± 0.00019
0.0016 ± 0.00013    0.0024 ± 0.00016
0.0060 ± 0.00025    0.0083 ± 0.00029
0.0013 ± 0.00011    0.0019 ± 0.00014
0.0006 ± 0.00008    0.0008 ± 0.00009
0.0003 ± 0.00005    0.0005 ± 0.00007
0.0010 ± 0.00010    0.0015 ± 0.00012
0.0026 ± 0.00016    0.0036 ± 0.00019
0.0005 ± 0.00007    0.0008 ± 0.00009
0.0010 ± 0.00010    0.0016 ± 0.00013
0.0005 ± 0.00007    0.0012 ± 0.00011
0.0024 ± 0.00016    0.0029 ± 0.00018
0.0017 ± 0.00013    0.0026 ± 0.00016
0.0002 ± 0.00005    0.0005 ± 0.00007
0.0006 ± 0.00008    0.0010 ± 0.00010
0.0003 ± 0.00006    0.0005 ± 0.00007
0.0016 ± 0.00013    0.0023 ± 0.00015
0.0062 ± 0.00025    0.0080 ± 0.00028
0.0011 ± 0.00010    0.0015 ± 0.00012
0.0022 ± 0.00015    0.0032 ± 0.00018
0.0014 ± 0.00012    0.0023 ± 0.00016
0.0055 ± 0.00024    0.0070 ± 0.00026

0.0065 ± 0.00026    0.0083 ± 0.00029
0.0012 ± 0.00011    0.0022 ± 0.00015
0.0024 ± 0.00016    0.0033 ± 0.00018
0.0019 ± 0.00014    0.0026 ± 0.00016
0.0066 ± 0.00026    0.0077 ± 0.00028
0.0054 ± 0.00023    0.0068 ± 0.00026
0.0010 ± 0.00010    0.0016 ± 0.00013
0.0022 ± 0.00015    0.0027 ± 0.00016
0.0014 ± 0.00012    0.0017 ± 0.00013
0.0046 ± 0.00021    0.0063 ± 0.00025

0.0052 ± 0.00023    0.0093 ± 0.00030
0.0055 ± 0.00023    0.0071 ± 0.00027
0.0012 ± 0.00011    0.0023 ± 0.00015
0.0011 ± 0.00010    0.0014 ± 0.00012
0.0027 ± 0.00017    0.0039 ± 0.00020
0.0021 ± 0.00015    0.0030 ± 0.00017
0.0016 ± 0.00013    0.0028 ± 0.00017
0.0015 ± 0.00012    0.0021 ± 0.00015
0.0065 ± 0.00025    0.0080 ± 0.00028
0.0050 ± 0.00022    0.0063 ± 0.00025
""".replace("±"," ").split()))

bn835.shape = (-1,4)

channels = ["K-p+", "K-p+p0", "K-p+p+p-", "K-p+p+"]
missing_channels = ["Ksp-p+", "K-K+", "K-K+p+"]
old_combinations = []
oldchannels_sp = [0,-1,1,2,-2]
oldchannels_sz = [3,-3]
old_eff_svd1 = np.zeros((4,4))
old_eff_svd2 = np.zeros((4,4))
old_err_svd1 = np.zeros((4,4))
old_err_svd2 = np.zeros((4,4))

for c1 in oldchannels_sp:
    for c2 in oldchannels_sp:
        if c1==-1 and c2==-1: continue
        old_combinations.append((c1,c2))
for c1 in oldchannels_sz:
    for c2 in oldchannels_sp:
        old_combinations.append((c1,c2))
for c1 in oldchannels_sp:
    for c2 in oldchannels_sz:
        old_combinations.append((c1,c2))

if len(bn835) != len(old_combinations):
    print "len not matching: ", len(bn835), len(old_combinations)
    #sys.exit(1)


for (c1,c2),(v1,e1,v2,e2) in zip(old_combinations, bn835):
    n1 = (c1>=0 and channels[c1]) or missing_channels[-1-c1]
    n2 = (c2>=0 and channels[c2]) or missing_channels[-1-c2]
    veto = c1>=0 and c2>=0 and " " or "*"
    print "%s%10s:%10s: %.4f %.4f" %(veto, n1, n2, v1, v2)
    if c1>=0 and c2>=0:
        old_eff_svd1[c1,c2] = v1
        old_eff_svd2[c1,c2] = v2
        old_err_svd1[c1,c2] = e1
        old_err_svd2[c1,c2] = e2

b0 = pr.chain(filenames,"B0")
b0g = pr.chain(filenames, "B0gen")

h_gen_svd1 = b0g.draw("channelP:channelM",cut="svdVs==0", option="text", range=(4,6,10,4,6,10))
h_gen_svd2 = b0g.draw("channelP:channelM",cut="svdVs==1", option="text", range=(4,6,10,4,6,10))

h_rec_svd1 = b0.draw("bestLHsig.channelP:bestLHsig.channelM",
                     cut="svdVs==0 && bestLHsig.mcInfo&1 && bestLHsig.flag==0 && bestLHsig.tag.flavour!=0", option="goff", range=(4,6,10,4,6,10))
h_rec_svd2 = b0.draw("bestLHsig.channelP:bestLHsig.channelM",
                     cut="svdVs==1 && bestLHsig.mcInfo&1 && bestLHsig.flag==0 && bestLHsig.tag.flavour!=0", option="goff", range=(4,6,10,4,6,10))


gen_svd1, dummy = r2mpl.getHistData(h_gen_svd1)
gen_svd2, dummy = r2mpl.getHistData(h_gen_svd2)
rec_svd1, dummy = r2mpl.getHistData(h_rec_svd1)
rec_svd2, dummy = r2mpl.getHistData(h_rec_svd2)

print "Generated:", gen_svd1.sum(), gen_svd2.sum()

eff_svd1 = rec_svd1/gen_svd1
eff_svd2 = rec_svd2/gen_svd2

diff_svd1 = (eff_svd1 - old_eff_svd1)/old_err_svd1
diff_svd2 = (eff_svd2 - old_eff_svd2)/old_err_svd2

for i,c1 in enumerate(channels):
    for j,c2 in enumerate(channels):
        print "%10s:%10s: %.4f %.4f vs %.4f %.4f = %+8.4f %+8.4f (%6d %6d)" % (
            c1, c2,
            eff_svd1[i,j], eff_svd2[i,j],
            old_eff_svd1[i,j], old_eff_svd2[i,j],
            diff_svd1[i,j], diff_svd2[i,j],
            gen_svd1[i,j], gen_svd2[i,j],
        )



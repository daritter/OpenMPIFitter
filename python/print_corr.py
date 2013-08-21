#!/usr/bin/python

import sys
import dspdsmks
import utils

fields = {
    "SVD1": [
        ("signal_corr_svd1_Mbc_mean",  "\Delta\mu(\mbc)"),
        ("signal_corr_svd1_Mbc_sigma", "\delta\sigma(\mbc)"),
        ("signal_corr_svd1_dE_mean",   "\Delta\mu(\de)"),
        ("signal_corr_svd1_dE_sigma",  "\delta\sigma(\de)"),
        ("signal_corr_svd1_rbin1",     "\delta f_Y(rbin0)"),
        ("signal_corr_svd1_rbin2",     "\delta f_Y(rbin1)"),
        ("signal_corr_svd1_rbin3",     "\delta f_Y(rbin2)"),
        ("signal_corr_svd1_rbin4",     "\delta f_Y(rbin3)"),
        ("signal_corr_svd1_rbin5",     "\delta f_Y(rbin4)"),
        ("signal_corr_svd1_rbin6",     "\delta f_Y(rbin5)"),
    ],
    "SVD2": [
        ("signal_corr_svd2_Mbc_mean",  "\Delta\mu(\mbc)"),
        ("signal_corr_svd2_Mbc_sigma", "\delta\sigma(\mbc)"),
        ("signal_corr_svd2_dE_mean",   "\Delta\mu(\de)"),
        ("signal_corr_svd2_dE_sigma",  "\delta\sigma(\de)"),
        ("signal_corr_svd2_rbin1",     "\delta f_Y(rbin0)"),
        ("signal_corr_svd2_rbin2",     "\delta f_Y(rbin1)"),
        ("signal_corr_svd2_rbin3",     "\delta f_Y(rbin2)"),
        ("signal_corr_svd2_rbin4",     "\delta f_Y(rbin3)"),
        ("signal_corr_svd2_rbin5",     "\delta f_Y(rbin4)"),
        ("signal_corr_svd2_rbin6",     "\delta f_Y(rbin5)"),
    ]
}

params = dspdsmks.Parameters()
params.load(sys.argv[1])


for svd, names in sorted(fields.items()):
    print r"""\begin{tabular}{LRCL}
    \toprule
    Name&Value\\
    \midrule"""
    for p,t in names:
        print r"    {0} & {1}\\" .format(t, utils.format_error(params(p).value, params(p).error, align=True))
    print r"""    \bottomrule
\end{tabular}"""



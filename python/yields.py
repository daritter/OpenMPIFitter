import sys
import dspdsmks
import utils
import ast
params = dspdsmks.Parameters()
params.load(sys.argv[1])
params("signal_svd1_nbb").value, params("signal_svd1_nbb").error = utils.nbb[0] #* 10
params("signal_svd2_nbb").value, params("signal_svd2_nbb").error = utils.nbb[1] #* 10
params("scale_bbar").value = 1
params("scale_continuum").value = 1
pdf = dspdsmks.DspDsmKsPDF()
yields = pdf.yields(params)
yields = ast.literal_eval(yields)

s_yield = yields["signal"]
del yields["signal"]

N1 = sum(e[0] for e in yields.values())
N2 = sum(e[1] for e in yields.values())

fractions1 = {e[0]:e[1][0]/N1 for e in yields.items()}
fractions2 = {e[0]:e[1][1]/N2 for e in yields.items()}

print """
yield_signal_svd1 %f %f -inf inf
yield_signal_svd2 %f %f -inf inf
yield_bkg_svd1  %f  %f -inf inf
yield_bkg_svd2  %f  %f -inf inf
f_bbar_svd1  %f 0 -inf inf
f_bbar_svd2  %f 0 -inf inf
f_misrecon_svd1 %f 0 -inf inf
f_misrecon_svd2 %f 0 -inf inf
f_continuum_svd1 %f 0 -inf inf
f_continuum_svd2 %f 0 -inf inf
""" % (
    s_yield[0], s_yield[0]/10, s_yield[1], s_yield[1]/10,
    N1, N1/10, N2, N2/10,
    fractions1["bbar"], fractions2["bbar"],
    fractions1["misrecon"], fractions2["misrecon"],
    fractions1["continuum"], fractions2["continuum"],
)

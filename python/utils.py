#coding: utf8
import math

def format_error(value,error, precision=3, exporder=3, exponent=None):
    if exponent is None: exponent = int(math.floor(math.log10(value)))
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

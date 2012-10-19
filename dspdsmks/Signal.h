#ifndef MPIFitter_Signal_h
#define MPIFitter_Signal_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"

namespace PAR {
    PARAM(sig_ratio);
    PARAM(sig_Mbc_mean_m0);
    PARAM(sig_Mbc_mean_m1);
    PARAM(sig_Mbc_sigma);
    PARAM(sig_Mbc_norm1);
    PARAM(sig_Mbc_meanshift1);
    PARAM(sig_Mbc_sigmascale1);
    PARAM(sig_Mbc_argusC);
    PARAM(sig_dE_mean);
    PARAM(sig_dE_sigma);
    PARAM(sig_dE_norm1);
    PARAM(sig_dE_meanshift1);
    PARAM(sig_dE_sigmascale1);
    PARAM(sig_dE_cheb1);
    PARAM(sig_dE_cheb2);
};


class Signal: public Component {
    public:
    Signal(double lowerMbc, double upperMbc, double lowerdE, double upperdE):
        Component(), signalPDF(lowerMbc,upperMbc, lowerdE, upperdE)
    {}

    virtual double operator()(const DspDsmKsEvent& e, const std::vector<double> &par) {
        //Set Parameters for signal component
        signalPDF.set(par[PAR::sig_ratio]);
        signalPDF.fcn1.fcnx.set(
                par[PAR::sig_Mbc_mean_m0] + e.dE*par[PAR::sig_Mbc_mean_m1],
                par[PAR::sig_Mbc_sigma], &par[PAR::sig_Mbc_norm1]
        );
        signalPDF.fcn1.fcny.set(&par[PAR::sig_dE_mean]);
        signalPDF.fcn2.fcnx.set(e.benergy, par[PAR::sig_Mbc_argusC]);
        signalPDF.fcn2.fcny.set(&par[PAR::sig_dE_cheb1]);

        return signalPDF(e.Mbc, e.dE);
    }

    private:

    /** PDF function components */
    Add2DFcn<
        CompoundFcn2D<MultiGauss<2>, MultiGauss<2> >,
        CompoundFcn2D<Argus, Chebychev<2> >
    > signalPDF;
};

#endif

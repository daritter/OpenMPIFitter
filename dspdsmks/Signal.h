#ifndef MPIFitter_Signal_h
#define MPIFitter_Signal_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"

namespace PAR {
    PARAM(yield_svd1_signal);
    PARAM(signal_svd1_ratio);
    PARAM(signal_svd1_Mbc_mean_m0);
    PARAM(signal_svd1_Mbc_mean_m1);
    PARAM(signal_svd1_Mbc_sigma);
    PARAM(signal_svd1_Mbc_norm1);
    PARAM(signal_svd1_Mbc_meanshift1);
    PARAM(signal_svd1_Mbc_sigmascale1);
    PARAM(signal_svd1_Mbc_argusC);
    PARAM(signal_svd1_dE_mean);
    PARAM(signal_svd1_dE_sigma);
    PARAM(signal_svd1_dE_norm1);
    PARAM(signal_svd1_dE_meanshift1);
    PARAM(signal_svd1_dE_sigmascale1);
    PARAM(signal_svd1_dE_cheb1);
    PARAM(signal_svd1_dE_cheb2);

    PARAM(yield_svd2_signal);
    PARAM(signal_svd2_ratio);
    PARAM(signal_svd2_Mbc_mean_m0);
    PARAM(signal_svd2_Mbc_mean_m1);
    PARAM(signal_svd2_Mbc_sigma);
    PARAM(signal_svd2_Mbc_norm1);
    PARAM(signal_svd2_Mbc_meanshift1);
    PARAM(signal_svd2_Mbc_sigmascale1);
    PARAM(signal_svd2_Mbc_argusC);
    PARAM(signal_svd2_dE_mean);
    PARAM(signal_svd2_dE_sigma);
    PARAM(signal_svd2_dE_norm1);
    PARAM(signal_svd2_dE_meanshift1);
    PARAM(signal_svd2_dE_sigmascale1);
    PARAM(signal_svd2_dE_cheb1);
    PARAM(signal_svd2_dE_cheb2);
};


class SignalPDF: public Component {
    public:
    SignalPDF(double lowerMbc, double upperMbc, double lowerdE, double upperdE):
        Component(),
        signalPDF_svd1(lowerMbc, upperMbc, lowerdE, upperdE),
        signalPDF_svd2(lowerMbc, upperMbc, lowerdE, upperdE)
    {}

    virtual double operator()(const DspDsmKsEvent& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            //Set Parameters for signal component
            signalPDF_svd1.set(par[PAR::signal_svd1_ratio]);
            signalPDF_svd1.fcn1.fcnx.set(
                    par[PAR::signal_svd1_Mbc_mean_m0] + e.dE*par[PAR::signal_svd1_Mbc_mean_m1],
                    par[PAR::signal_svd1_Mbc_sigma], &par[PAR::signal_svd1_Mbc_norm1]
            );
            signalPDF_svd1.fcn1.fcny.set(&par[PAR::signal_svd1_dE_mean]);
            signalPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::signal_svd1_Mbc_argusC]);
            signalPDF_svd1.fcn2.fcny.set(&par[PAR::signal_svd1_dE_cheb1]);

            return par[PAR::yield_svd1_signal] * signalPDF_svd1(e.Mbc, e.dE);
        }else{
            //Set Parameters for signal component
            signalPDF_svd2.set(par[PAR::signal_svd2_ratio]);
            signalPDF_svd2.fcn1.fcnx.set(
                    par[PAR::signal_svd2_Mbc_mean_m0] + e.dE*par[PAR::signal_svd2_Mbc_mean_m1],
                    par[PAR::signal_svd2_Mbc_sigma], &par[PAR::signal_svd2_Mbc_norm1]
            );
            signalPDF_svd2.fcn1.fcny.set(&par[PAR::signal_svd2_dE_mean]);
            signalPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::signal_svd2_Mbc_argusC]);
            signalPDF_svd2.fcn2.fcny.set(&par[PAR::signal_svd2_dE_cheb1]);

            return par[PAR::yield_svd2_signal] * signalPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH){
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::yield_svd1_signal];
        }
        if(svd & SVD2){
            yield += par[PAR::yield_svd2_signal];
        }
        return yield;
    }

    private:

    /** PDF function components */
    Add2DFcn<
        CompoundFcn2D<MultiGauss<2>, MultiGauss<2> >,
        CompoundFcn2D<Argus, Chebychev<2> >
    > signalPDF_svd1;

    Add2DFcn<
        CompoundFcn2D<MultiGauss<2>, MultiGauss<2> >,
        CompoundFcn2D<Argus, Chebychev<2> >
    > signalPDF_svd2;
};

#endif

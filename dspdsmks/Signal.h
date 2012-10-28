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

    PARAM(signal_dt_blifetime);
    PARAM(signal_dt_Jc);
    PARAM(signal_dt_Js1);
    PARAM(signal_dt_Js2);
    PARAM(signal_dt_fractionscale);
};


class SignalPDF: public Component {
    public:
    SignalPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        Component(range_dT, false, useDeltaT),
        signalPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        signalPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setParameters(
                PAR::signal_dt_blifetime, PAR::signal_dt_Jc, PAR::signal_dt_Js1, PAR::signal_dt_Js2,
                PAR::signal_dt_fractionscale, useDeltaT?Event::dt_signal:-1);
    }

    virtual double operator()(const Event& e, const std::vector<double> &par) {
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

            return get_deltaT(e,par)*par[PAR::yield_svd1_signal] * signalPDF_svd1(e.Mbc, e.dE);
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

            return get_deltaT(e,par)*par[PAR::yield_svd2_signal] * signalPDF_svd2(e.Mbc, e.dE);
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

    virtual Event get_maxEvent(const std::vector<double> &par, int svd){
        Event e = Component::get_maxEvent(par, svd);
        if(svd == SVD1){
            e.svdVs = 0;
            e.dE  = par[PAR::signal_svd1_dE_mean];
            e.Mbc = par[PAR::signal_svd1_Mbc_mean_m0] + e.dE*par[PAR::signal_svd1_Mbc_mean_m1];
        }else if(svd == SVD2){
            e.svdVs = 1;
            e.dE  = par[PAR::signal_svd2_dE_mean];
            e.Mbc = par[PAR::signal_svd2_Mbc_mean_m0] + e.dE*par[PAR::signal_svd2_Mbc_mean_m1];
        }
        return e;
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

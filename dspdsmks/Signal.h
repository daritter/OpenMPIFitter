#ifndef MPIFitter_Signal_h
#define MPIFitter_Signal_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"

namespace PAR {
    PARAM(yield_signal_br);

    PARAM(signal_svd1_rbin1);
    PARAM(signal_svd1_rbin2);
    PARAM(signal_svd1_rbin3);
    PARAM(signal_svd1_rbin4);
    PARAM(signal_svd1_rbin5);
    PARAM(signal_svd1_rbin6);

    PARAM(signal_svd1_eff);
    PARAM(signal_svd1_nbb);
    PARAM(signal_svd1_ratio);
    PARAM(signal_svd1_Mbc_mean_m0);
    PARAM(signal_svd1_Mbc_mean_m1);
    PARAM(signal_svd1_Mbc_mean_m2);
    PARAM(signal_svd1_Mbc_sigma_m0);
    PARAM(signal_svd1_Mbc_sigma_m1);
    PARAM(signal_svd1_Mbc_sigma_m2);
    PARAM(signal_svd1_Mbc_ratio);
    PARAM(signal_svd1_Mbc_ratio_m1);
    PARAM(signal_svd1_Mbc_ratio_m2);
    PARAM(signal_svd1_Mbc_meanshift);
    PARAM(signal_svd1_Mbc_meanshift_m1);
    PARAM(signal_svd1_Mbc_meanshift_m2);
    PARAM(signal_svd1_Mbc_sigmascale);
    PARAM(signal_svd1_Mbc_sigmascale_m1);
    PARAM(signal_svd1_Mbc_sigmascale_m2);
    PARAM(signal_svd1_Mbc_sigma2);
    PARAM(signal_svd1_Mbc_sigma2_m1);
    PARAM(signal_svd1_Mbc_sigma2_m2);
    PARAM(signal_svd1_Mbc_sigma2scale);
    PARAM(signal_svd1_Mbc_argusC);
    PARAM(signal_svd1_dE_mean);
    PARAM(signal_svd1_dE_sigma);
    PARAM(signal_svd1_dE_norm1);
    PARAM(signal_svd1_dE_meanshift1);
    PARAM(signal_svd1_dE_sigmascale1);
    PARAM(signal_svd1_dE_norm2);
    PARAM(signal_svd1_dE_meanshift2);
    PARAM(signal_svd1_dE_sigmascale2);
    PARAM(signal_svd1_dE_bkg_mean);
    PARAM(signal_svd1_dE_bkg_sigma);
    PARAM(signal_svd1_Mbc_g3_ratio);
    PARAM(signal_svd1_Mbc_g3_mean);
    PARAM(signal_svd1_Mbc_g3_sigma);
    PARAM(signal_svd1_Mbc_g3_sigma2);

    PARAM(signal_svd2_rbin1);
    PARAM(signal_svd2_rbin2);
    PARAM(signal_svd2_rbin3);
    PARAM(signal_svd2_rbin4);
    PARAM(signal_svd2_rbin5);
    PARAM(signal_svd2_rbin6);

    PARAM(signal_svd2_eff);
    PARAM(signal_svd2_nbb);
    PARAM(signal_svd2_ratio);
    PARAM(signal_svd2_Mbc_mean_m0);
    PARAM(signal_svd2_Mbc_mean_m1);
    PARAM(signal_svd2_Mbc_mean_m2);
    PARAM(signal_svd2_Mbc_sigma_m0);
    PARAM(signal_svd2_Mbc_sigma_m1);
    PARAM(signal_svd2_Mbc_sigma_m2);
    PARAM(signal_svd2_Mbc_ratio);
    PARAM(signal_svd2_Mbc_ratio_m1);
    PARAM(signal_svd2_Mbc_ratio_m2);
    PARAM(signal_svd2_Mbc_meanshift);
    PARAM(signal_svd2_Mbc_meanshift_m1);
    PARAM(signal_svd2_Mbc_meanshift_m2);
    PARAM(signal_svd2_Mbc_sigmascale);
    PARAM(signal_svd2_Mbc_sigmascale_m1);
    PARAM(signal_svd2_Mbc_sigmascale_m2);
    PARAM(signal_svd2_Mbc_sigma2);
    PARAM(signal_svd2_Mbc_sigma2_m1);
    PARAM(signal_svd2_Mbc_sigma2_m2);
    PARAM(signal_svd2_Mbc_sigma2scale);
    PARAM(signal_svd2_Mbc_argusC);
    PARAM(signal_svd2_dE_mean);
    PARAM(signal_svd2_dE_sigma);
    PARAM(signal_svd2_dE_norm1);
    PARAM(signal_svd2_dE_meanshift1);
    PARAM(signal_svd2_dE_sigmascale1);
    PARAM(signal_svd2_dE_norm2);
    PARAM(signal_svd2_dE_meanshift2);
    PARAM(signal_svd2_dE_sigmascale2);
    PARAM(signal_svd2_dE_bkg_mean);
    PARAM(signal_svd2_dE_bkg_sigma);
    PARAM(signal_svd2_Mbc_g3_ratio);
    PARAM(signal_svd2_Mbc_g3_mean);
    PARAM(signal_svd2_Mbc_g3_sigma);
    PARAM(signal_svd2_Mbc_g3_sigma2);

    PARAM(signal_dt_blifetime);
    PARAM(signal_dt_Jc);
    PARAM(signal_dt_Js1);
    PARAM(signal_dt_Js2);
    PARAM(signal_dt_fractionscale);
};


class SignalPDF: public DeltaTComponent<> {
    public:
    SignalPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        DeltaTComponent<>(range_dT, false, useDeltaT, false), range_mBC(range_mBC), range_dE(range_dE),
        signalPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        signalPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setParameters(
                PAR::signal_dt_blifetime, PAR::signal_dt_Jc, PAR::signal_dt_Js1, PAR::signal_dt_Js2,
                PAR::signal_dt_fractionscale, -1, useDeltaT?Event::dt_signal:-1);
    }
    virtual ~SignalPDF(){}

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            //Set Parameters for signal component
            signalPDF_svd1.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            signalPDF_svd1.set(par[PAR::signal_svd1_ratio]);

            signalPDF_svd1.fcn1.fcnx.set(par[PAR::signal_svd1_Mbc_g3_ratio]);
            signalPDF_svd1.fcn1.fcnx.fcn1.set(
                    std::max(0.0,std::min(1.0,par[PAR::signal_svd1_Mbc_ratio] + e.dE*par[PAR::signal_svd1_Mbc_ratio_m1] + e.dE*e.dE*par[PAR::signal_svd1_Mbc_ratio_m2])),
                    par[PAR::signal_svd1_Mbc_mean_m0] + e.dE*par[PAR::signal_svd1_Mbc_mean_m1] + e.dE*e.dE*par[PAR::signal_svd1_Mbc_mean_m2],
                    par[PAR::signal_svd1_Mbc_meanshift] + e.dE*par[PAR::signal_svd1_Mbc_meanshift_m1] + e.dE*e.dE*par[PAR::signal_svd1_Mbc_meanshift_m2],
                    std::max(0.0,par[PAR::signal_svd1_Mbc_sigma_m0] + e.dE*par[PAR::signal_svd1_Mbc_sigma_m1] + e.dE*e.dE*par[PAR::signal_svd1_Mbc_sigma_m2]),
                    std::max(0.0,par[PAR::signal_svd1_Mbc_sigmascale] + e.dE*par[PAR::signal_svd1_Mbc_sigmascale_m1] + e.dE*e.dE*par[PAR::signal_svd1_Mbc_sigmascale_m2]),
                    std::max(0.0,par[PAR::signal_svd1_Mbc_sigma2] + e.dE*par[PAR::signal_svd1_Mbc_sigma2_m1] + e.dE*e.dE*par[PAR::signal_svd1_Mbc_sigma2_m2]),
                    par[PAR::signal_svd1_Mbc_sigma2scale]
            );
            signalPDF_svd1.fcn1.fcnx.fcn2.set(par[PAR::signal_svd1_Mbc_g3_mean], par[PAR::signal_svd1_Mbc_g3_sigma], par[PAR::signal_svd1_Mbc_g3_sigma2]);

            signalPDF_svd1.fcn1.fcny.set(&par[PAR::signal_svd1_dE_mean]);
            signalPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::signal_svd1_Mbc_argusC]);
            signalPDF_svd1.fcn2.fcny.set(par[PAR::signal_svd1_dE_bkg_mean], par[PAR::signal_svd1_dE_bkg_sigma]);

            return get_deltaT(e,par)* get_yield(par, SVD1, e.rbin) * signalPDF_svd1(e.Mbc, e.dE);
        }else{
            //Set Parameters for signal component
            //signalPDF_svd2.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            signalPDF_svd2.set(par[PAR::signal_svd2_ratio]);

            signalPDF_svd2.fcn1.fcnx.set(par[PAR::signal_svd2_Mbc_g3_ratio]);
            signalPDF_svd2.fcn1.fcnx.fcn1.set(
                    std::max(0.0,std::min(1.0,par[PAR::signal_svd2_Mbc_ratio] + e.dE*par[PAR::signal_svd2_Mbc_ratio_m1] + e.dE*e.dE*par[PAR::signal_svd2_Mbc_ratio_m2])),
                    par[PAR::signal_svd2_Mbc_mean_m0] + e.dE*par[PAR::signal_svd2_Mbc_mean_m1] + e.dE*e.dE*par[PAR::signal_svd2_Mbc_mean_m2],
                    par[PAR::signal_svd2_Mbc_meanshift] + e.dE*par[PAR::signal_svd2_Mbc_meanshift_m1] + e.dE*e.dE*par[PAR::signal_svd2_Mbc_meanshift_m2],
                    std::max(0.0,par[PAR::signal_svd2_Mbc_sigma_m0] + e.dE*par[PAR::signal_svd2_Mbc_sigma_m1] + e.dE*e.dE*par[PAR::signal_svd2_Mbc_sigma_m2]),
                    std::max(0.0,par[PAR::signal_svd2_Mbc_sigmascale] + e.dE*par[PAR::signal_svd2_Mbc_sigmascale_m1] + e.dE*e.dE*par[PAR::signal_svd2_Mbc_sigmascale_m2]),
                    std::max(0.0,par[PAR::signal_svd2_Mbc_sigma2] + e.dE*par[PAR::signal_svd2_Mbc_sigma2_m1] + e.dE*e.dE*par[PAR::signal_svd2_Mbc_sigma2_m2]),
                    par[PAR::signal_svd2_Mbc_sigma2scale]
            );
            signalPDF_svd2.fcn1.fcnx.fcn2.set(par[PAR::signal_svd2_Mbc_g3_mean], par[PAR::signal_svd2_Mbc_g3_sigma], par[PAR::signal_svd2_Mbc_g3_sigma2]);

            signalPDF_svd2.fcn1.fcny.set(&par[PAR::signal_svd2_dE_mean]);
            signalPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::signal_svd2_Mbc_argusC]);
            signalPDF_svd2.fcn2.fcny.set(par[PAR::signal_svd2_dE_bkg_mean], par[PAR::signal_svd2_dE_bkg_sigma]);

            return get_deltaT(e,par)* get_yield(par, SVD2, e.rbin) * signalPDF_svd2(e.Mbc, e.dE);
        }
    }

    static double get_signal_yield(const std::vector<double> &par, EnabledSVD svd, int rbin=-1) {
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::yield_signal_br] * par[PAR::signal_svd1_nbb] * par[PAR::signal_svd1_eff] * get_rbinFraction(rbin, PAR::signal_svd1_rbin1, par);
        }
        if(svd & SVD2){
            yield += par[PAR::yield_signal_br] * par[PAR::signal_svd2_nbb] * par[PAR::signal_svd2_eff] * get_rbinFraction(rbin, PAR::signal_svd2_rbin1, par);
        }
        return yield;
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const {
        return SignalPDF::get_signal_yield(par,svd,rbin);
    }

    private:

    Range range_mBC;
    Range range_dE;

    /** PDF function components */
    Add2DFcn<
        CompoundFcn2D<Add1DFcn<DoubleGauss, Gauss>, MultiGauss<3> >,
        CompoundFcn2D<Argus, Gauss>
    > signalPDF_svd1;

    Add2DFcn<
        CompoundFcn2D<Add1DFcn<DoubleGauss, Gauss>, MultiGauss<3> >,
        CompoundFcn2D<Argus, Gauss>
    > signalPDF_svd2;
};

#endif

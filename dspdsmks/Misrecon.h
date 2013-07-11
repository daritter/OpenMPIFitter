#ifndef MPIFitter_Misrecon_h
#define MPIFitter_Misrecon_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"
#include "Signal.h"

namespace PAR {
    PARAM(ratio_misrecon_svd1);
    PARAM(misrecon_svd1_rbin1);
    PARAM(misrecon_svd1_rbin2);
    PARAM(misrecon_svd1_rbin3);
    PARAM(misrecon_svd1_rbin4);
    PARAM(misrecon_svd1_rbin5);
    PARAM(misrecon_svd1_rbin6);
    PARAM(misrecon_svd1_rbin7);

    PARAM(misrecon_svd1_ratio);
    PARAM(misrecon_svd1_Mbc_mean);
    PARAM(misrecon_svd1_Mbc_sigma);
    PARAM(misrecon_svd1_Mbc_sigma2);
    PARAM(misrecon_svd1_Mbc_argusC);
    PARAM(misrecon_svd1_dE_mean);
    PARAM(misrecon_svd1_dE_sigma);
    PARAM(misrecon_svd1_dE_bkg_mean);
    PARAM(misrecon_svd1_dE_bkg_sigma);

    PARAM(ratio_misrecon_svd2);
    PARAM(misrecon_svd2_rbin1);
    PARAM(misrecon_svd2_rbin2);
    PARAM(misrecon_svd2_rbin3);
    PARAM(misrecon_svd2_rbin4);
    PARAM(misrecon_svd2_rbin5);
    PARAM(misrecon_svd2_rbin6);
    PARAM(misrecon_svd2_rbin7);

    PARAM(misrecon_svd2_ratio);
    PARAM(misrecon_svd2_Mbc_mean);
    PARAM(misrecon_svd2_Mbc_sigma);
    PARAM(misrecon_svd2_Mbc_sigma2);
    PARAM(misrecon_svd2_Mbc_argusC);
    PARAM(misrecon_svd2_dE_mean);
    PARAM(misrecon_svd2_dE_sigma);
    PARAM(misrecon_svd2_dE_bkg_mean);
    PARAM(misrecon_svd2_dE_bkg_sigma);

    PARAM(misrecon_dt_blifetime);
    PARAM(misrecon_dt_fractionscale);
    PARAM(misrecon_dt_offset);
};


class MisreconPDF: public DeltaTComponent<> {
    public:
    MisreconPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT, bool eta_dependence):
        DeltaTComponent<>(range_dT, false, useDeltaT), range_mBC(range_mBC), range_dE(range_dE),
        misreconPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        misreconPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setParameters(
                PAR::misrecon_dt_blifetime, PAR::signal_dt_Jc, PAR::signal_dt_Js1, PAR::signal_dt_Js2,
                PAR::misrecon_dt_fractionscale, -1, Event::dt_misrecon);
        deltaT.setOffset(PAR::misrecon_dt_offset);
        deltaT.setEtaDependence(eta_dependence);
    }

    virtual ~MisreconPDF(){}

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            //Set Parameters for misrecon component
            misreconPDF_svd1.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            misreconPDF_svd1.set(par[PAR::misrecon_svd1_ratio]);
            misreconPDF_svd1.fcn1.fcnx.set(par[PAR::misrecon_svd1_Mbc_mean],par[PAR::misrecon_svd1_Mbc_sigma],par[PAR::misrecon_svd1_Mbc_sigma2]);
            misreconPDF_svd1.fcn1.fcny.set(&par[PAR::misrecon_svd1_dE_mean]);
            misreconPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::misrecon_svd1_Mbc_argusC]);
            misreconPDF_svd1.fcn2.fcny.set(par[PAR::misrecon_svd1_dE_bkg_mean], par[PAR::misrecon_svd1_dE_bkg_sigma]);

            return get_deltaT(e,par) * get_yield(par, SVD1, e.rbin) * misreconPDF_svd1(e.Mbc, e.dE);
        }else{
            //Set Parameters for misrecon component
            misreconPDF_svd2.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            misreconPDF_svd2.set(par[PAR::misrecon_svd2_ratio]);
            misreconPDF_svd2.fcn1.fcnx.set(par[PAR::misrecon_svd2_Mbc_mean],par[PAR::misrecon_svd2_Mbc_sigma],par[PAR::misrecon_svd2_Mbc_sigma2]);
            misreconPDF_svd2.fcn1.fcny.set(&par[PAR::misrecon_svd2_dE_mean]);
            misreconPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::misrecon_svd2_Mbc_argusC]);
            misreconPDF_svd2.fcn2.fcny.set(par[PAR::misrecon_svd2_dE_bkg_mean], par[PAR::misrecon_svd2_dE_bkg_sigma]);

            return get_deltaT(e,par) * get_yield(par, SVD2, e.rbin) * misreconPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const {
        double yield(0);
        if(svd & SVD1){
            yield += SignalPDF::get_signal_yield(par, SVD1) * par[PAR::ratio_misrecon_svd1] * get_rbinFraction(rbin, PAR::misrecon_svd1_rbin1, par);
        }
        if(svd & SVD2){
            yield += SignalPDF::get_signal_yield(par, SVD2) * par[PAR::ratio_misrecon_svd2] * get_rbinFraction(rbin, PAR::misrecon_svd2_rbin1, par);
        }
        return yield;
    }

    virtual void get_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, double scale=1.0){
        fill_rbinFractions(par,fractions, svd, PAR::misrecon_svd1_rbin1, PAR::misrecon_svd2_rbin1, scale);
    }

    private:

    Range range_mBC;
    Range range_dE;

    /** PDF function components */
    Add2DFcn<
        CompoundFcn2D<Gauss, MultiGauss<1> >,
        CompoundFcn2D<Argus, Gauss>
    > misreconPDF_svd1;

    Add2DFcn<
        CompoundFcn2D<Gauss, MultiGauss<1> >,
        CompoundFcn2D<Argus, Gauss>
    > misreconPDF_svd2;
};

#endif

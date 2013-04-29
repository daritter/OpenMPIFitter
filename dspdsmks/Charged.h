#ifndef MPIFitter_Charged_h
#define MPIFitter_Charged_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"
#include "Mixed.h"

namespace PAR {
    PARAM(scale_charged);
    PARAM(ratio_charged_svd1);
    PARAM(charged_svd1_ratio);
    PARAM(charged_svd1_Mbc_mean);
    PARAM(charged_svd1_Mbc_sigma);
    PARAM(charged_svd1_Mbc_argusC);
    PARAM(charged_svd1_dE_mean);
    PARAM(charged_svd1_dE_sigma);
    PARAM(charged_svd1_dE_cheb1);

    PARAM(ratio_charged_svd2);
    PARAM(charged_svd2_ratio);
    PARAM(charged_svd2_Mbc_mean);
    PARAM(charged_svd2_Mbc_sigma);
    PARAM(charged_svd2_Mbc_argusC);
    PARAM(charged_svd2_dE_mean);
    PARAM(charged_svd2_dE_sigma);
    PARAM(charged_svd2_dE_cheb1);
    PARAM(charged_svd2_dE_cheb2);
    PARAM(charged_svd2_ctrlpeak_ratio);
    PARAM(charged_svd2_ctrlpeak_Mbc_mean);
    PARAM(charged_svd2_ctrlpeak_Mbc_sigma);
    PARAM(charged_svd2_ctrlpeak_dE_mean);
    PARAM(charged_svd2_ctrlpeak_dE_sigma);

};


class ChargedPDF: public DeltaTComponent<BkgTPDF> {
    public:
    ChargedPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        DeltaTComponent<BkgTPDF>(range_dT, true, useDeltaT), range_mBC(range_mBC), range_dE(range_dE),
        chargedPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        chargedPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setCommonParameters(PAR::bkg_dt_mean_delta, PAR::bkg_dt_blifetime, PAR::bkg_dt_mean_tau,
                PAR::bkg_dt_outlier_fraction, PAR::bkg_dt_outlier_mean, PAR::bkg_dt_outlier_scale);
        deltaT.setParameters(PAR::bkg_dt_sigma_main_sgl, PAR::bkg_dt_sigma_tail_sgl,
                PAR::bkg_dt_fraction_delta_sgl, PAR::bkg_dt_fraction_tail_sgl, false);
        deltaT.setParameters(PAR::bkg_dt_sigma_main_mul, PAR::bkg_dt_sigma_tail_mul,
                PAR::bkg_dt_fraction_delta_mul, PAR::bkg_dt_fraction_tail_mul, true);
        deltaT.setAcp(PAR::bkg_dt_acp0);
    }

    virtual ~ChargedPDF(){}

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            //Set Parameters for charged component
            chargedPDF_svd1.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            chargedPDF_svd1.set(par[PAR::charged_svd1_ratio]);
            chargedPDF_svd1.fcn1.fcnx.set(par[PAR::charged_svd1_Mbc_mean], par[PAR::charged_svd1_Mbc_sigma]);
            chargedPDF_svd1.fcn1.fcny.set(par[PAR::charged_svd1_dE_mean], par[PAR::charged_svd1_dE_sigma]);
            chargedPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::charged_svd1_Mbc_argusC]);
            chargedPDF_svd1.fcn2.fcny.set(&par[PAR::charged_svd1_dE_cheb1]);

            return get_deltaT(e,par) * get_yield(par, SVD1, e.rbin) * chargedPDF_svd1(e.Mbc, e.dE);
        } else {
            //Set Parameters for charged component
            chargedPDF_svd2.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            chargedPDF_svd2.set(par[PAR::charged_svd2_ratio]);
            chargedPDF_svd2.fcn1.set(par[PAR::charged_svd2_ctrlpeak_ratio]);
            chargedPDF_svd2.fcn1.fcn1.fcnx.set(par[PAR::charged_svd2_Mbc_mean], par[PAR::charged_svd2_Mbc_sigma]);
            chargedPDF_svd2.fcn1.fcn1.fcny.set(par[PAR::charged_svd2_dE_mean], par[PAR::charged_svd2_dE_sigma]);
            if(par[PAR::charged_svd2_ctrlpeak_ratio]<1){
                chargedPDF_svd2.fcn1.fcn2.fcnx.set(par[PAR::charged_svd2_ctrlpeak_Mbc_mean], par[PAR::charged_svd2_ctrlpeak_Mbc_sigma]);
                chargedPDF_svd2.fcn1.fcn2.fcny.set(par[PAR::charged_svd2_ctrlpeak_dE_mean], par[PAR::charged_svd2_ctrlpeak_dE_sigma]);
            }
            chargedPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::charged_svd2_Mbc_argusC]);
            chargedPDF_svd2.fcn2.fcny.set(&par[PAR::charged_svd2_dE_cheb1]);

            return get_deltaT(e,par) * get_yield(par, SVD2, e.rbin) * chargedPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const {
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::ratio_charged_svd1] * par[PAR::yield_mixed_svd1] * get_rbinFraction(rbin, PAR::bkg_svd1_rbin1, par);
        }
        if(svd & SVD2){
            yield += par[PAR::ratio_charged_svd2] * par[PAR::yield_mixed_svd2] * get_rbinFraction(rbin, PAR::bkg_svd2_rbin1, par);
        }
        return par[PAR::scale_charged] * yield;
    }

    virtual void get_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, double scale=1.0){
        fill_rbinFractions(par, fractions, svd, PAR::bkg_svd1_rbin1, PAR::bkg_svd2_rbin1, scale);
    }


    private:

    Range range_mBC;
    Range range_dE;

    /** PDF function components */
    Add2DFcn<
        CompoundFcn2D<Gauss, Gauss>,
        CompoundFcn2D<Argus, Chebychev<1> >
    > chargedPDF_svd1;
    Add2DFcn<
        Add2DFcn< CompoundFcn2D<Gauss, Gauss>, CompoundFcn2D<Gauss, Gauss> >,
        CompoundFcn2D<Argus, Chebychev<2> >
    > chargedPDF_svd2;
};

#endif

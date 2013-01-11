#ifndef MPIFitter_Charged_h
#define MPIFitter_Charged_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"
#include "Mixed.h"

namespace PAR {
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

    PARAM(charged_dt_blifetime);
    PARAM(charged_dt_sgl_mean1);
    PARAM(charged_dt_sgl_mean2);
    PARAM(charged_dt_sgl_mean3);
    PARAM(charged_dt_sgl_mean4);
    PARAM(charged_dt_sgl_sigma1);
    PARAM(charged_dt_sgl_sigma2);
    PARAM(charged_dt_sgl_sigma3);
    PARAM(charged_dt_sgl_sigma4);
    PARAM(charged_dt_sgl_weight2);
    PARAM(charged_dt_sgl_weight3);
    PARAM(charged_dt_sgl_weight4);
    PARAM(charged_dt_mul_mean1);
    PARAM(charged_dt_mul_mean2);
    PARAM(charged_dt_mul_mean3);
    PARAM(charged_dt_mul_mean4);
    PARAM(charged_dt_mul_sigma1);
    PARAM(charged_dt_mul_sigma2);
    PARAM(charged_dt_mul_sigma3);
    PARAM(charged_dt_mul_sigma4);
    PARAM(charged_dt_mul_weight2);
    PARAM(charged_dt_mul_weight3);
    PARAM(charged_dt_mul_weight4);
    PARAM(charged_dt_fractionscale);
    PARAM(charged_dt_outliermean);
    PARAM(charged_dt_outlierscale);
    PARAM(charged_dt_outlierscale2);
    PARAM(charged_dt_outlierratio);
};


class ChargedPDF: public DeltaTComponent<GenericTPDF> {
    public:
    ChargedPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        DeltaTComponent<GenericTPDF>(range_dT, true, useDeltaT), range_mBC(range_mBC), range_dE(range_dE),
        chargedPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        chargedPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setCommonParameters(PAR::charged_dt_blifetime, PAR::charged_dt_fractionscale,
                PAR::charged_dt_outliermean, PAR::charged_dt_outlierscale, PAR::charged_dt_outlierscale2, PAR::charged_dt_outlierratio);
        deltaT.setParameters(false,
                PAR::charged_dt_sgl_mean1, PAR::charged_dt_sgl_mean2, PAR::charged_dt_sgl_mean3, PAR::charged_dt_sgl_mean4,
                PAR::charged_dt_sgl_sigma1, PAR::charged_dt_sgl_sigma2, PAR::charged_dt_sgl_sigma3,  PAR::charged_dt_sgl_sigma4,
                PAR::charged_dt_sgl_weight2, PAR::charged_dt_sgl_weight3, PAR::charged_dt_sgl_weight4);
        deltaT.setParameters(true,
                PAR::charged_dt_mul_mean1, PAR::charged_dt_mul_mean2, PAR::charged_dt_mul_mean3, PAR::charged_dt_mul_mean4,
                PAR::charged_dt_mul_sigma1, PAR::charged_dt_mul_sigma2, PAR::charged_dt_mul_sigma3, PAR::charged_dt_mul_sigma4,
                PAR::charged_dt_mul_weight2, PAR::charged_dt_mul_weight3, PAR::charged_dt_mul_weight4);
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

            return get_deltaT(e,par)*get_yield(par, Component::SVD1)*chargedPDF_svd1(e.Mbc, e.dE);
        } else {
            //Set Parameters for charged component
            chargedPDF_svd2.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            chargedPDF_svd2.set(par[PAR::charged_svd2_ratio]);
            chargedPDF_svd2.fcn1.fcnx.set(par[PAR::charged_svd2_Mbc_mean], par[PAR::charged_svd2_Mbc_sigma]);
            chargedPDF_svd2.fcn1.fcny.set(par[PAR::charged_svd2_dE_mean], par[PAR::charged_svd2_dE_sigma]);
            chargedPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::charged_svd2_Mbc_argusC]);
            chargedPDF_svd2.fcn2.fcny.set(&par[PAR::charged_svd2_dE_cheb1]);

            return get_deltaT(e,par)*get_yield(par, Component::SVD2)*chargedPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH){
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::ratio_charged_svd1]*par[PAR::yield_mixed_svd1];
        }
        if(svd & SVD2){
            yield += par[PAR::ratio_charged_svd2]*par[PAR::yield_mixed_svd2];
        }
        return yield;
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
        CompoundFcn2D<Gauss, Gauss>,
        CompoundFcn2D<Argus, Chebychev<1> >
    > chargedPDF_svd2;
};

#endif

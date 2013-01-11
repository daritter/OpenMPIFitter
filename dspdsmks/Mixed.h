#ifndef MPIFitter_Mixed_h
#define MPIFitter_Mixed_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"
#include "GenericT.h"

namespace PAR {
    PARAM(yield_mixed_svd1);
    PARAM(mixed_svd1_ratio);
    PARAM(mixed_svd1_Mbc_mean);
    PARAM(mixed_svd1_Mbc_sigma);
    PARAM(mixed_svd1_Mbc_argusC);
    PARAM(mixed_svd1_dE_mean);
    PARAM(mixed_svd1_dE_sigma);
    PARAM(mixed_svd1_dE_cheb1);

    PARAM(yield_mixed_svd2);
    PARAM(mixed_svd2_ratio);
    PARAM(mixed_svd2_Mbc_mean);
    PARAM(mixed_svd2_Mbc_sigma);
    PARAM(mixed_svd2_Mbc_argusC);
    PARAM(mixed_svd2_dE_mean);
    PARAM(mixed_svd2_dE_sigma);
    PARAM(mixed_svd2_dE_cheb1);

    PARAM(mixed_dt_blifetime);
    PARAM(mixed_dt_sgl_mean1);
    PARAM(mixed_dt_sgl_mean2);
    PARAM(mixed_dt_sgl_mean3);
    PARAM(mixed_dt_sgl_mean4);
    PARAM(mixed_dt_sgl_sigma1);
    PARAM(mixed_dt_sgl_sigma2);
    PARAM(mixed_dt_sgl_sigma3);
    PARAM(mixed_dt_sgl_sigma4);
    PARAM(mixed_dt_sgl_weight2);
    PARAM(mixed_dt_sgl_weight3);
    PARAM(mixed_dt_sgl_weight4);
    PARAM(mixed_dt_mul_mean1);
    PARAM(mixed_dt_mul_mean2);
    PARAM(mixed_dt_mul_mean3);
    PARAM(mixed_dt_mul_mean4);
    PARAM(mixed_dt_mul_sigma1);
    PARAM(mixed_dt_mul_sigma2);
    PARAM(mixed_dt_mul_sigma3);
    PARAM(mixed_dt_mul_sigma4);
    PARAM(mixed_dt_mul_weight2);
    PARAM(mixed_dt_mul_weight3);
    PARAM(mixed_dt_mul_weight4);
    PARAM(mixed_dt_fractionscale);
    PARAM(mixed_dt_outliermean);
    PARAM(mixed_dt_outlierscale);
    PARAM(mixed_dt_outlierscale2);
    PARAM(mixed_dt_outlierratio);
};


class MixedPDF: public DeltaTComponent<GenericTPDF> {
    public:
    MixedPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        DeltaTComponent<GenericTPDF>(range_dT, false, useDeltaT), range_mBC(range_mBC), range_dE(range_dE),
        mixedPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        mixedPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setCommonParameters(PAR::mixed_dt_blifetime, PAR::mixed_dt_fractionscale,
                PAR::mixed_dt_outliermean, PAR::mixed_dt_outlierscale, PAR::mixed_dt_outlierscale2, PAR::mixed_dt_outlierratio);
        deltaT.setParameters(false,
                PAR::mixed_dt_sgl_mean1, PAR::mixed_dt_sgl_mean2, PAR::mixed_dt_sgl_mean3, PAR::mixed_dt_sgl_mean4,
                PAR::mixed_dt_sgl_sigma1, PAR::mixed_dt_sgl_sigma2, PAR::mixed_dt_sgl_sigma3,  PAR::mixed_dt_sgl_sigma4,
                PAR::mixed_dt_sgl_weight2, PAR::mixed_dt_sgl_weight3, PAR::mixed_dt_sgl_weight4);
        deltaT.setParameters(true,
                PAR::mixed_dt_mul_mean1, PAR::mixed_dt_mul_mean2, PAR::mixed_dt_mul_mean3, PAR::mixed_dt_mul_mean4,
                PAR::mixed_dt_mul_sigma1, PAR::mixed_dt_mul_sigma2, PAR::mixed_dt_mul_sigma3, PAR::mixed_dt_mul_sigma4,
                PAR::mixed_dt_mul_weight2, PAR::mixed_dt_mul_weight3, PAR::mixed_dt_mul_weight4);
    }

    virtual ~MixedPDF(){}

    virtual double operator()(const Event& e, const std::vector<double> &par){
        if(e.svdVs==0){
            //Set Parameters for mixed component
            mixedPDF_svd1.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            mixedPDF_svd1.set(par[PAR::mixed_svd1_ratio]);
            mixedPDF_svd1.fcn1.fcnx.set(par[PAR::mixed_svd1_Mbc_mean], par[PAR::mixed_svd1_Mbc_sigma]);
            mixedPDF_svd1.fcn1.fcny.set(par[PAR::mixed_svd1_dE_mean], par[PAR::mixed_svd1_dE_sigma]);
            mixedPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::mixed_svd1_Mbc_argusC]);
            mixedPDF_svd1.fcn2.fcny.set(&par[PAR::mixed_svd1_dE_cheb1]);

            return get_deltaT(e,par)*par[PAR::yield_mixed_svd1] * mixedPDF_svd1(e.Mbc, e.dE);
        }else{
            //Set Parameters for mixed component
            mixedPDF_svd2.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            mixedPDF_svd2.set(par[PAR::mixed_svd2_ratio]);
            mixedPDF_svd2.fcn1.fcnx.set(par[PAR::mixed_svd2_Mbc_mean], par[PAR::mixed_svd2_Mbc_sigma]);
            mixedPDF_svd2.fcn1.fcny.set(par[PAR::mixed_svd2_dE_mean], par[PAR::mixed_svd2_dE_sigma]);
            mixedPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::mixed_svd2_Mbc_argusC]);
            mixedPDF_svd2.fcn2.fcny.set(&par[PAR::mixed_svd2_dE_cheb1]);

            return get_deltaT(e,par)*par[PAR::yield_mixed_svd2] * mixedPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH){
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::yield_mixed_svd1];
        }
        if(svd & SVD2){
            yield += par[PAR::yield_mixed_svd2];
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
    > mixedPDF_svd1;
    Add2DFcn<
        CompoundFcn2D<Gauss, Gauss>,
        CompoundFcn2D<Argus, Chebychev<1> >
    > mixedPDF_svd2;
};

#endif

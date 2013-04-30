#ifndef MPIFitter_Mixed_h
#define MPIFitter_Mixed_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"
#include "BkgT.h"

namespace PAR {
    PARAM(scale_mixed);
    PARAM(yield_mixed_svd1);
    PARAM(mixed_svd1_ratio);
    PARAM(mixed_svd1_Mbc_mean);
    PARAM(mixed_svd1_Mbc_sigma);
    PARAM(mixed_svd1_Mbc_argusC);
    PARAM(mixed_svd1_dE_mean);
    PARAM(mixed_svd1_dE_sigma);
    PARAM(mixed_svd1_dE_cheb1);
    PARAM(mixed_svd1_ctrlpeak_ratio);
    PARAM(mixed_svd1_ctrlpeak_Mbc_mean);
    PARAM(mixed_svd1_ctrlpeak_Mbc_sigma);
    PARAM(mixed_svd1_ctrlpeak_dE_mean);
    PARAM(mixed_svd1_ctrlpeak_dE_sigma);

    PARAM(yield_mixed_svd2);
    PARAM(mixed_svd2_ratio);
    PARAM(mixed_svd2_Mbc_mean);
    PARAM(mixed_svd2_Mbc_sigma);
    PARAM(mixed_svd2_Mbc_argusC);
    PARAM(mixed_svd2_dE_mean);
    PARAM(mixed_svd2_dE_sigma);
    PARAM(mixed_svd2_dE_cheb1);
    PARAM(mixed_svd2_dE_cheb2);
    PARAM(mixed_svd2_ctrlpeak_ratio);
    PARAM(mixed_svd2_ctrlpeak_Mbc_mean);
    PARAM(mixed_svd2_ctrlpeak_Mbc_sigma);
    PARAM(mixed_svd2_ctrlpeak_dE_mean);
    PARAM(mixed_svd2_ctrlpeak_dE_sigma);

    PARAM(bkg_svd2_rbin1);
    PARAM(bkg_svd2_rbin2);
    PARAM(bkg_svd2_rbin3);
    PARAM(bkg_svd2_rbin4);
    PARAM(bkg_svd2_rbin5);
    PARAM(bkg_svd2_rbin6);
    PARAM(bkg_svd2_rbin7);
    PARAM(bkg_svd1_rbin1);
    PARAM(bkg_svd1_rbin2);
    PARAM(bkg_svd1_rbin3);
    PARAM(bkg_svd1_rbin4);
    PARAM(bkg_svd1_rbin5);
    PARAM(bkg_svd1_rbin6);
    PARAM(bkg_svd1_rbin7);
    PARAM(bkg_dt_blifetime);
    PARAM(bkg_dt_mean_delta);
    PARAM(bkg_dt_mean_tau);
    PARAM(bkg_dt_sigma_main_sgl);
    PARAM(bkg_dt_sigma_tail_sgl);
    PARAM(bkg_dt_fraction_delta_sgl);
    PARAM(bkg_dt_fraction_tail_sgl);
    PARAM(bkg_dt_sigma_main_mul);
    PARAM(bkg_dt_sigma_tail_mul);
    PARAM(bkg_dt_fraction_delta_mul);
    PARAM(bkg_dt_fraction_tail_mul);
    PARAM(bkg_dt_outlier_fraction);
    PARAM(bkg_dt_outlier_mean);
    PARAM(bkg_dt_outlier_scale);
    PARAM(bkg_dt_acp0);
    PARAM(bkg_dt_acp1);
    PARAM(bkg_dt_acp2);
    PARAM(bkg_dt_acp3);
    PARAM(bkg_dt_acp4);
    PARAM(bkg_dt_acp5);
    PARAM(bkg_dt_acp6);

};


class MixedPDF: public DeltaTComponent<BkgTPDF> {
    public:
    MixedPDF(Range range_mBC, Range range_dE, Range range_dT, bool useMbcdE = true, bool useDeltaT=false):
        DeltaTComponent<BkgTPDF>(range_dT, false, useDeltaT), useMbcdE(useMbcdE), range_mBC(range_mBC), range_dE(range_dE),
        mixedPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        mixedPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setCommonParameters(PAR::bkg_dt_mean_delta, PAR::bkg_dt_blifetime, PAR::bkg_dt_mean_tau,
                PAR::bkg_dt_outlier_fraction, PAR::bkg_dt_outlier_mean, PAR::bkg_dt_outlier_scale);
        deltaT.setParameters(PAR::bkg_dt_sigma_main_sgl, PAR::bkg_dt_sigma_tail_sgl,
                PAR::bkg_dt_fraction_delta_sgl, PAR::bkg_dt_fraction_tail_sgl, false);
        deltaT.setParameters(PAR::bkg_dt_sigma_main_mul, PAR::bkg_dt_sigma_tail_mul,
                PAR::bkg_dt_fraction_delta_mul, PAR::bkg_dt_fraction_tail_mul, true);
        deltaT.setAcp(PAR::bkg_dt_acp0);
    }

    virtual ~MixedPDF(){}

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            if(useMbcdE){
                //Set Parameters for mixed component
                mixedPDF_svd1.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
                mixedPDF_svd1.set(par[PAR::mixed_svd1_ratio]);
                mixedPDF_svd1.fcn1.set(par[PAR::mixed_svd1_ctrlpeak_ratio]);
                mixedPDF_svd1.fcn1.fcn1.fcnx.set(par[PAR::mixed_svd1_Mbc_mean], par[PAR::mixed_svd1_Mbc_sigma]);
                mixedPDF_svd1.fcn1.fcn1.fcny.set(par[PAR::mixed_svd1_dE_mean], par[PAR::mixed_svd1_dE_sigma]);
                if(par[PAR::mixed_svd1_ctrlpeak_ratio]<1){
                    mixedPDF_svd1.fcn1.fcn2.fcnx.set(par[PAR::mixed_svd1_ctrlpeak_Mbc_mean], par[PAR::mixed_svd1_ctrlpeak_Mbc_sigma]);
                    mixedPDF_svd1.fcn1.fcn2.fcny.set(par[PAR::mixed_svd1_ctrlpeak_dE_mean], par[PAR::mixed_svd1_ctrlpeak_dE_sigma]);
                }
                mixedPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::mixed_svd1_Mbc_argusC]);
                mixedPDF_svd1.fcn2.fcny.set(&par[PAR::mixed_svd1_dE_cheb1]);
            }

            return get_deltaT(e,par) * get_yield(par, SVD1, e.rbin) * (useMbcdE?mixedPDF_svd1(e.Mbc, e.dE):1);
        } else {
            if(useMbcdE){
                //Set Parameters for mixed component
                mixedPDF_svd2.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
                mixedPDF_svd2.set(par[PAR::mixed_svd2_ratio]);
                mixedPDF_svd2.fcn1.set(par[PAR::mixed_svd2_ctrlpeak_ratio]);
                mixedPDF_svd2.fcn1.fcn1.fcnx.set(par[PAR::mixed_svd2_Mbc_mean], par[PAR::mixed_svd2_Mbc_sigma]);
                mixedPDF_svd2.fcn1.fcn1.fcny.set(par[PAR::mixed_svd2_dE_mean], par[PAR::mixed_svd2_dE_sigma]);
                if(par[PAR::mixed_svd2_ctrlpeak_ratio]<1){
                    mixedPDF_svd2.fcn1.fcn2.fcnx.set(par[PAR::mixed_svd2_ctrlpeak_Mbc_mean], par[PAR::mixed_svd2_ctrlpeak_Mbc_sigma]);
                    mixedPDF_svd2.fcn1.fcn2.fcny.set(par[PAR::mixed_svd2_ctrlpeak_dE_mean], par[PAR::mixed_svd2_ctrlpeak_dE_sigma]);
                }
                mixedPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::mixed_svd2_Mbc_argusC]);
                mixedPDF_svd2.fcn2.fcny.set(&par[PAR::mixed_svd2_dE_cheb1]);
            }

            return get_deltaT(e,par) * get_yield(par, SVD2, e.rbin) * (useMbcdE?mixedPDF_svd2(e.Mbc, e.dE):1);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const {
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::yield_mixed_svd1] * get_rbinFraction(rbin, PAR::bkg_svd1_rbin1, par);
        }
        if(svd & SVD2){
            yield += par[PAR::yield_mixed_svd2] * get_rbinFraction(rbin, PAR::bkg_svd2_rbin1, par);
        }
        return par[PAR::scale_mixed] * yield;
    }

    virtual void get_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, double scale=1.0){
        fill_rbinFractions(par, fractions, svd, PAR::bkg_svd1_rbin1, PAR::bkg_svd2_rbin1, scale);
    }

    private:

    bool useMbcdE;

    Range range_mBC;
    Range range_dE;

    /** PDF function components */
    Add2DFcn<
        Add2DFcn< CompoundFcn2D<Gauss, Gauss>, CompoundFcn2D<Gauss, Gauss> >,
        CompoundFcn2D<Argus, Chebychev<1> >
    > mixedPDF_svd1;
    Add2DFcn<
        Add2DFcn< CompoundFcn2D<Gauss, Gauss>, CompoundFcn2D<Gauss, Gauss> >,
        CompoundFcn2D<Argus, Chebychev<2> >
    > mixedPDF_svd2;
};

#endif

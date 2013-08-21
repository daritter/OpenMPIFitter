#ifndef MPIFitter_BBar_h
#define MPIFitter_BBar_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"
#include "BkgT.h"

namespace PAR {
    PARAM(scale_bbar);
    PARAM(yield_bbar_svd1);
    PARAM(bbar_svd1_ratio);
    PARAM(bbar_svd1_Mbc_mean);
    PARAM(bbar_svd1_Mbc_sigma);
    PARAM(bbar_svd1_Mbc_argusC);
    PARAM(bbar_svd1_dE_mean);
    PARAM(bbar_svd1_dE_sigma);
    PARAM(bbar_svd1_dE_cheb1);
    PARAM(bbar_svd1_ctrlpeak_ratio);
    PARAM(bbar_svd1_ctrlpeak_Mbc_mean);
    PARAM(bbar_svd1_ctrlpeak_Mbc_sigma);
    PARAM(bbar_svd1_ctrlpeak_dE_mean);
    PARAM(bbar_svd1_ctrlpeak_dE_sigma);

    PARAM(yield_bbar_svd2);
    PARAM(bbar_svd2_ratio);
    PARAM(bbar_svd2_Mbc_mean);
    PARAM(bbar_svd2_Mbc_sigma);
    PARAM(bbar_svd2_Mbc_argusC);
    PARAM(bbar_svd2_dE_mean);
    PARAM(bbar_svd2_dE_sigma);
    PARAM(bbar_svd2_dE_cheb1);
    PARAM(bbar_svd2_dE_cheb2);
    PARAM(bbar_svd2_ctrlpeak_ratio);
    PARAM(bbar_svd2_ctrlpeak_Mbc_mean);
    PARAM(bbar_svd2_ctrlpeak_Mbc_sigma);
    PARAM(bbar_svd2_ctrlpeak_dE_mean);
    PARAM(bbar_svd2_ctrlpeak_dE_sigma);

    PARAM(bkg_svd1_rbin1);
    PARAM(bkg_svd1_rbin2);
    PARAM(bkg_svd1_rbin3);
    PARAM(bkg_svd1_rbin4);
    PARAM(bkg_svd1_rbin5);
    PARAM(bkg_svd1_rbin6);
    PARAM(bkg_svd1_rbin7);

    PARAM(bkg_svd2_rbin1);
    PARAM(bkg_svd2_rbin2);
    PARAM(bkg_svd2_rbin3);
    PARAM(bkg_svd2_rbin4);
    PARAM(bkg_svd2_rbin5);
    PARAM(bkg_svd2_rbin6);
    PARAM(bkg_svd2_rbin7);

    PARAM(bkg_dt_blifetime);
    PARAM(bkg_svd1_dt_mean_delta);
    PARAM(bkg_svd1_dt_mean_tau);
    PARAM(bkg_svd1_dt_sigma_main_sgl);
    PARAM(bkg_svd1_dt_sigma_tail1_sgl);
    PARAM(bkg_svd1_dt_sigma_tail2_sgl);
    PARAM(bkg_svd1_dt_weight_tail1_sgl);
    PARAM(bkg_svd1_dt_weight_tail2_sgl);
    PARAM(bkg_svd1_dt_fraction_delta_sgl);
    PARAM(bkg_svd1_dt_sigma_main_mul);
    PARAM(bkg_svd1_dt_sigma_tail1_mul);
    PARAM(bkg_svd1_dt_sigma_tail2_mul);
    PARAM(bkg_svd1_dt_weight_tail1_mul);
    PARAM(bkg_svd1_dt_weight_tail2_mul);
    PARAM(bkg_svd1_dt_fraction_delta_mul);
    PARAM(bkg_svd1_dt_outlier_fraction);
    PARAM(bkg_svd1_dt_outlier_mean);
    PARAM(bkg_svd1_dt_outlier_scale);
    PARAM(bkg_svd1_dt_acp0);
    PARAM(bkg_svd1_dt_acp1);
    PARAM(bkg_svd1_dt_acp2);
    PARAM(bkg_svd1_dt_acp3);
    PARAM(bkg_svd1_dt_acp4);
    PARAM(bkg_svd1_dt_acp5);
    PARAM(bkg_svd1_dt_acp6);

    PARAM(bkg_svd2_dt_mean_delta);
    PARAM(bkg_svd2_dt_mean_tau);
    PARAM(bkg_svd2_dt_sigma_main_sgl);
    PARAM(bkg_svd2_dt_sigma_tail1_sgl);
    PARAM(bkg_svd2_dt_sigma_tail2_sgl);
    PARAM(bkg_svd2_dt_weight_tail1_sgl);
    PARAM(bkg_svd2_dt_weight_tail2_sgl);
    PARAM(bkg_svd2_dt_fraction_delta_sgl);
    PARAM(bkg_svd2_dt_sigma_main_mul);
    PARAM(bkg_svd2_dt_sigma_tail1_mul);
    PARAM(bkg_svd2_dt_sigma_tail2_mul);
    PARAM(bkg_svd2_dt_weight_tail1_mul);
    PARAM(bkg_svd2_dt_weight_tail2_mul);
    PARAM(bkg_svd2_dt_fraction_delta_mul);
    PARAM(bkg_svd2_dt_outlier_fraction);
    PARAM(bkg_svd2_dt_outlier_mean);
    PARAM(bkg_svd2_dt_outlier_scale);
    PARAM(bkg_svd2_dt_acp0);
    PARAM(bkg_svd2_dt_acp1);
    PARAM(bkg_svd2_dt_acp2);
    PARAM(bkg_svd2_dt_acp3);
    PARAM(bkg_svd2_dt_acp4);
    PARAM(bkg_svd2_dt_acp5);
    PARAM(bkg_svd2_dt_acp6);
};


class BBarPDF: public DeltaTComponent<BkgTPDF> {
    public:
    static void init_deltaT(BkgTPDF &deltaT, bool combined, bool eta_dependence){
        /*if(!combined){
            deltaT.setCommonParameters(0,PAR::bkg_svd1_dt_mean_delta, PAR::bkg_svd1_dt_blifetime, PAR::bkg_svd1_dt_mean_tau,
                    PAR::bkg_svd1_dt_outlier_fraction, PAR::bkg_svd1_dt_outlier_mean, PAR::bkg_svd1_dt_outlier_scale);
            deltaT.setParameters(0,PAR::bkg_svd1_dt_sigma_main_sgl, PAR::bkg_svd1_dt_sigma_tail_sgl,
                    PAR::bkg_svd1_dt_fraction_delta_sgl, PAR::bkg_svd1_dt_fraction_tail_sgl, false);
            deltaT.setParameters(0,PAR::bkg_svd1_dt_sigma_main_mul, PAR::bkg_svd1_dt_sigma_tail_mul,
                    PAR::bkg_svd1_dt_fraction_delta_mul, PAR::bkg_svd1_dt_fraction_tail_mul, true);
            deltaT.setAcp(0,PAR::bkg_svd1_dt_acp0);
        }else{
            deltaT.setCommonParameters(0,PAR::bkg_svd2_dt_mean_delta, PAR::bkg_svd2_dt_blifetime, PAR::bkg_svd2_dt_mean_tau,
                    PAR::bkg_svd2_dt_outlier_fraction, PAR::bkg_svd2_dt_outlier_mean, PAR::bkg_svd2_dt_outlier_scale);
            deltaT.setParameters(0,PAR::bkg_svd2_dt_sigma_main_sgl, PAR::bkg_svd2_dt_sigma_tail_sgl,
                    PAR::bkg_svd2_dt_fraction_delta_sgl, PAR::bkg_svd2_dt_fraction_tail_sgl, false);
            deltaT.setParameters(0,PAR::bkg_svd2_dt_sigma_main_mul, PAR::bkg_svd2_dt_sigma_tail_mul,
                    PAR::bkg_svd2_dt_fraction_delta_mul, PAR::bkg_svd2_dt_fraction_tail_mul, true);
            deltaT.setAcp(0,PAR::bkg_svd2_dt_acp0);
        }*/


        deltaT.setCommonParameters(0,PAR::bkg_svd1_dt_mean_delta, PAR::bkg_dt_blifetime, PAR::bkg_svd1_dt_mean_tau,
                PAR::bkg_svd1_dt_outlier_fraction, PAR::bkg_svd1_dt_outlier_mean, PAR::bkg_svd1_dt_outlier_scale);
        deltaT.setParameters(0,PAR::bkg_svd1_dt_sigma_main_sgl, PAR::bkg_svd1_dt_sigma_tail1_sgl, PAR::bkg_svd1_dt_sigma_tail2_sgl,
                PAR::bkg_svd1_dt_weight_tail1_sgl, PAR::bkg_svd1_dt_weight_tail2_sgl,
                PAR::bkg_svd1_dt_fraction_delta_sgl, false);
        if(!combined) {
            deltaT.setParameters(0,PAR::bkg_svd1_dt_sigma_main_mul, PAR::bkg_svd1_dt_sigma_tail1_mul, PAR::bkg_svd1_dt_sigma_tail2_mul,
                PAR::bkg_svd1_dt_weight_tail1_mul, PAR::bkg_svd1_dt_weight_tail2_mul,
                PAR::bkg_svd1_dt_fraction_delta_mul, true);
        }
        deltaT.setAcp(0,PAR::bkg_svd1_dt_acp0);

        deltaT.setCommonParameters(1,PAR::bkg_svd2_dt_mean_delta, PAR::bkg_dt_blifetime, PAR::bkg_svd2_dt_mean_tau,
                PAR::bkg_svd2_dt_outlier_fraction, PAR::bkg_svd2_dt_outlier_mean, PAR::bkg_svd2_dt_outlier_scale);
        deltaT.setParameters(1,PAR::bkg_svd2_dt_sigma_main_sgl, PAR::bkg_svd2_dt_sigma_tail1_sgl, PAR::bkg_svd2_dt_sigma_tail2_sgl,
                PAR::bkg_svd2_dt_weight_tail1_sgl, PAR::bkg_svd2_dt_weight_tail2_sgl,
                PAR::bkg_svd2_dt_fraction_delta_sgl, false);
        deltaT.setParameters(1,PAR::bkg_svd2_dt_sigma_main_mul, PAR::bkg_svd2_dt_sigma_tail1_mul, PAR::bkg_svd2_dt_sigma_tail2_mul,
                PAR::bkg_svd2_dt_weight_tail1_mul, PAR::bkg_svd2_dt_weight_tail2_mul,
                PAR::bkg_svd2_dt_fraction_delta_mul, true);
        deltaT.setAcp(1,PAR::bkg_svd2_dt_acp0);

        deltaT.setEtaDependence(eta_dependence);
    }

    BBarPDF(Range range_mBC, Range range_dE, Range range_dT, bool useMbcdE, bool useDeltaT, bool combinedDeltaT, bool eta_dependence):
        DeltaTComponent<BkgTPDF>(range_dT, false, useDeltaT), useMbcdE(useMbcdE), range_mBC(range_mBC), range_dE(range_dE),
        bbarPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        bbarPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        init_deltaT(deltaT,combinedDeltaT, eta_dependence);
    }

    virtual ~BBarPDF(){}

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            if(useMbcdE){
                //Set Parameters for bbar component
                bbarPDF_svd1.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
                bbarPDF_svd1.set(par[PAR::bbar_svd1_ratio]);
                bbarPDF_svd1.fcn1.set(par[PAR::bbar_svd1_ctrlpeak_ratio]);
                bbarPDF_svd1.fcn1.fcn1.fcnx.set(par[PAR::bbar_svd1_Mbc_mean], par[PAR::bbar_svd1_Mbc_sigma]);
                bbarPDF_svd1.fcn1.fcn1.fcny.set(par[PAR::bbar_svd1_dE_mean], par[PAR::bbar_svd1_dE_sigma]);
                if(par[PAR::bbar_svd1_ctrlpeak_ratio]<1){
                    bbarPDF_svd1.fcn1.fcn2.fcnx.set(par[PAR::bbar_svd1_ctrlpeak_Mbc_mean], par[PAR::bbar_svd1_ctrlpeak_Mbc_sigma]);
                    bbarPDF_svd1.fcn1.fcn2.fcny.set(par[PAR::bbar_svd1_ctrlpeak_dE_mean], par[PAR::bbar_svd1_ctrlpeak_dE_sigma]);
                }
                bbarPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::bbar_svd1_Mbc_argusC]);
                bbarPDF_svd1.fcn2.fcny.set(&par[PAR::bbar_svd1_dE_cheb1]);
            }

            return get_deltaT(e,par) * get_yield(par, SVD1, e.rbin) * (useMbcdE?bbarPDF_svd1(e.Mbc, e.dE):1);
        } else {
            if(useMbcdE){
                //Set Parameters for bbar component
                bbarPDF_svd2.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
                bbarPDF_svd2.set(par[PAR::bbar_svd2_ratio]);
                bbarPDF_svd2.fcn1.set(par[PAR::bbar_svd2_ctrlpeak_ratio]);
                bbarPDF_svd2.fcn1.fcn1.fcnx.set(par[PAR::bbar_svd2_Mbc_mean], par[PAR::bbar_svd2_Mbc_sigma]);
                bbarPDF_svd2.fcn1.fcn1.fcny.set(par[PAR::bbar_svd2_dE_mean], par[PAR::bbar_svd2_dE_sigma]);
                if(par[PAR::bbar_svd2_ctrlpeak_ratio]<1){
                    bbarPDF_svd2.fcn1.fcn2.fcnx.set(par[PAR::bbar_svd2_ctrlpeak_Mbc_mean], par[PAR::bbar_svd2_ctrlpeak_Mbc_sigma]);
                    bbarPDF_svd2.fcn1.fcn2.fcny.set(par[PAR::bbar_svd2_ctrlpeak_dE_mean], par[PAR::bbar_svd2_ctrlpeak_dE_sigma]);
                }
                bbarPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::bbar_svd2_Mbc_argusC]);
                bbarPDF_svd2.fcn2.fcny.set(&par[PAR::bbar_svd2_dE_cheb1]);
            }

            return get_deltaT(e,par) * get_yield(par, SVD2, e.rbin) * (useMbcdE?bbarPDF_svd2(e.Mbc, e.dE):1);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const {
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::yield_bbar_svd1] * get_rbinFraction(rbin, PAR::bkg_svd1_rbin1, par);
        }
        if(svd & SVD2){
            yield += par[PAR::yield_bbar_svd2] * get_rbinFraction(rbin, PAR::bkg_svd2_rbin1, par);
        }
        return par[PAR::scale_bbar] * yield;
    }

    virtual void get_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, double scale=1.0) const {
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
    > bbarPDF_svd1;
    Add2DFcn<
        Add2DFcn< CompoundFcn2D<Gauss, Gauss>, CompoundFcn2D<Gauss, Gauss> >,
        CompoundFcn2D<Argus, Chebychev<2> >
    > bbarPDF_svd2;
};

#endif

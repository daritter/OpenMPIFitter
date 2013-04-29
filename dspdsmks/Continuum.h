#ifndef MPIFitter_Continuum_h
#define MPIFitter_Continuum_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"
#include "Mixed.h"

namespace PAR {
    PARAM(scale_continuum);
    PARAM(ratio_continuum_svd1);
    PARAM(continuum_svd1_Mbc_argusC);
    PARAM(continuum_svd1_dE_cheb1);

    PARAM(ratio_continuum_svd2);
    PARAM(continuum_svd2_Mbc_argusC);
    PARAM(continuum_svd2_dE_cheb1);
};


class ContinuumPDF: public DeltaTComponent<BkgTPDF> {
    public:
    ContinuumPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        DeltaTComponent<BkgTPDF>(range_dT, true, useDeltaT), range_mBC(range_mBC), range_dE(range_dE),
        continuumPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        continuumPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setCommonParameters(PAR::bkg_dt_mean_delta, PAR::bkg_dt_blifetime, PAR::bkg_dt_mean_tau,
                PAR::bkg_dt_outlier_fraction, PAR::bkg_dt_outlier_mean, PAR::bkg_dt_outlier_scale);
        deltaT.setParameters(PAR::bkg_dt_sigma_main_sgl, PAR::bkg_dt_sigma_tail_sgl,
                PAR::bkg_dt_fraction_delta_sgl, PAR::bkg_dt_fraction_tail_sgl, false);
        deltaT.setParameters(PAR::bkg_dt_sigma_main_mul, PAR::bkg_dt_sigma_tail_mul,
                PAR::bkg_dt_fraction_delta_mul, PAR::bkg_dt_fraction_tail_mul, true);
        deltaT.setAcp(PAR::bkg_dt_acp0);
    }

    virtual ~ContinuumPDF(){}

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            //Set Parameters for continuum component
            continuumPDF_svd1.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            continuumPDF_svd1.fcnx.set(e.benergy, par[PAR::continuum_svd1_Mbc_argusC]);
            continuumPDF_svd1.fcny.set(&par[PAR::continuum_svd1_dE_cheb1]);

            return get_deltaT(e,par) * get_yield(par, SVD1, e.rbin) * continuumPDF_svd1(e.Mbc, e.dE);
        } else {
            //Set Parameters for continuum component
            continuumPDF_svd2.set_limits(range_mBC.vmin, std::min(e.benergy,(double) range_mBC.vmax), range_dE.vmin, range_dE.vmax);
            continuumPDF_svd2.fcnx.set(e.benergy, par[PAR::continuum_svd2_Mbc_argusC]);
            continuumPDF_svd2.fcny.set(&par[PAR::continuum_svd2_dE_cheb1]);

            return get_deltaT(e,par) * get_yield(par, SVD2, e.rbin) * continuumPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const {
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::ratio_continuum_svd1] * par[PAR::yield_mixed_svd1] * get_rbinFraction(rbin, PAR::bkg_svd1_rbin1, par);
        }
        if(svd & SVD2){
            yield += par[PAR::ratio_continuum_svd2] * par[PAR::yield_mixed_svd2] * get_rbinFraction(rbin, PAR::bkg_svd2_rbin1, par);
        }
        return par[PAR::scale_continuum] * yield;
    }

    virtual void get_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, double scale=1.0){
        fill_rbinFractions(par, fractions, svd, PAR::bkg_svd1_rbin1, PAR::bkg_svd2_rbin1, scale);
    }


    private:

    Range range_mBC;
    Range range_dE;

    /** PDF function components */
    CompoundFcn2D<Argus, Chebychev<1> > continuumPDF_svd1;
    CompoundFcn2D<Argus, Chebychev<1> > continuumPDF_svd2;
};

#endif

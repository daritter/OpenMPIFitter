#ifndef MPIFitter_Misrecon_h
#define MPIFitter_Misrecon_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"

namespace PAR {
    PARAM(yield_svd1_misrecon);
    PARAM(misrecon_svd1_ratio);
    PARAM(misrecon_svd1_Mbc_mean);
    //PARAM(misrecon_svd1_Mbc_mean_m1);
    PARAM(misrecon_svd1_Mbc_sigma);
    //PARAM(misrecon_svd1_Mbc_norm1);
    //PARAM(misrecon_svd1_Mbc_meanshift1);
    //PARAM(misrecon_svd1_Mbc_sigmascale1);
    PARAM(misrecon_svd1_Mbc_argusC);
    PARAM(misrecon_svd1_dE_mean);
    PARAM(misrecon_svd1_dE_sigma);
    //PARAM(misrecon_svd1_dE_norm1);
    //PARAM(misrecon_svd1_dE_meanshift1);
    //PARAM(misrecon_svd1_dE_sigmascale1);
    PARAM(misrecon_svd1_dE_bkg_mean);
    PARAM(misrecon_svd1_dE_bkg_sigma);

    PARAM(yield_svd2_misrecon);
    PARAM(misrecon_svd2_ratio);
    PARAM(misrecon_svd2_Mbc_mean);
    //PARAM(misrecon_svd2_Mbc_mean_m1);
    PARAM(misrecon_svd2_Mbc_sigma);
    //PARAM(misrecon_svd2_Mbc_norm1);
    //PARAM(misrecon_svd2_Mbc_meanshift1);
    //PARAM(misrecon_svd2_Mbc_sigmascale1);
    PARAM(misrecon_svd2_Mbc_argusC);
    PARAM(misrecon_svd2_dE_mean);
    PARAM(misrecon_svd2_dE_sigma);
    //PARAM(misrecon_svd2_dE_norm1);
    //PARAM(misrecon_svd2_dE_meanshift1);
    //PARAM(misrecon_svd2_dE_sigmascale1);
    PARAM(misrecon_svd2_dE_bkg_mean);
    PARAM(misrecon_svd2_dE_bkg_sigma);

    PARAM(misrecon_dt_blifetime);
    PARAM(misrecon_dt_Jc);
    PARAM(misrecon_dt_Js1);
    PARAM(misrecon_dt_Js2);
    PARAM(misrecon_dt_fractionscale);
};


class MisreconPDF: public Component {
    public:
    MisreconPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        Component(range_dT, false, useDeltaT),
        misreconPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        misreconPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setParameters(
                PAR::misrecon_dt_blifetime, PAR::misrecon_dt_Jc, PAR::misrecon_dt_Js1, PAR::misrecon_dt_Js2,
                PAR::misrecon_dt_fractionscale, useDeltaT?Event::dt_misrecon:-1);
    }

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            //Set Parameters for misrecon component
            misreconPDF_svd1.set(par[PAR::misrecon_svd1_ratio]);
            misreconPDF_svd1.fcn1.fcnx.set(&par[PAR::misrecon_svd1_Mbc_mean]);
            misreconPDF_svd1.fcn1.fcny.set(&par[PAR::misrecon_svd1_dE_mean]);
            misreconPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::misrecon_svd1_Mbc_argusC]);
            misreconPDF_svd1.fcn2.fcny.set(par[PAR::misrecon_svd1_dE_bkg_mean], par[PAR::misrecon_svd1_dE_bkg_sigma]);

            return get_deltaT(e,par)*par[PAR::yield_svd1_misrecon] * misreconPDF_svd1(e.Mbc, e.dE);
        }else{
            //Set Parameters for misrecon component
            misreconPDF_svd2.set(par[PAR::misrecon_svd2_ratio]);
            misreconPDF_svd2.fcn1.fcnx.set(&par[PAR::misrecon_svd2_Mbc_mean]);
            misreconPDF_svd2.fcn1.fcny.set(&par[PAR::misrecon_svd2_dE_mean]);
            misreconPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::misrecon_svd2_Mbc_argusC]);
            misreconPDF_svd2.fcn2.fcny.set(par[PAR::misrecon_svd2_dE_bkg_mean], par[PAR::misrecon_svd2_dE_bkg_sigma]);
            return get_deltaT(e,par)*par[PAR::yield_svd2_misrecon] * misreconPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH){
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::yield_svd1_misrecon];
        }
        if(svd & SVD2){
            yield += par[PAR::yield_svd2_misrecon];
        }
        return yield;
    }

    private:

    /** PDF function components */
    Add2DFcn<
        CompoundFcn2D<MultiGauss<1>, MultiGauss<1> >,
        CompoundFcn2D<Argus, Gauss>
    > misreconPDF_svd1;

    Add2DFcn<
        CompoundFcn2D<MultiGauss<1>, MultiGauss<1> >,
        CompoundFcn2D<Argus, Gauss>
    > misreconPDF_svd2;
};

#endif

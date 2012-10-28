#ifndef MPIFitter_Charged_h
#define MPIFitter_Charged_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"

namespace PAR {
    PARAM(yield_svd1_charged);
    PARAM(charged_svd1_ratio);
    PARAM(charged_svd1_Mbc_mean);
    PARAM(charged_svd1_Mbc_sigma);
    PARAM(charged_svd1_Mbc_argusC);
    PARAM(charged_svd1_dE_mean);
    PARAM(charged_svd1_dE_sigma);
    PARAM(charged_svd1_dE_cheb1);

    PARAM(yield_svd2_charged);
    PARAM(charged_svd2_ratio);
    PARAM(charged_svd2_Mbc_mean);
    PARAM(charged_svd2_Mbc_sigma);
    PARAM(charged_svd2_Mbc_argusC);
    PARAM(charged_svd2_dE_mean);
    PARAM(charged_svd2_dE_sigma);
    PARAM(charged_svd2_dE_cheb1);

    PARAM(charged_dt_blifetime);
    PARAM(charged_dt_Jc);
    PARAM(charged_dt_Js1);
    PARAM(charged_dt_Js2);
    PARAM(charged_dt_fractionscale);
};


class ChargedPDF: public Component {
    public:
    ChargedPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        Component(range_dT, true, useDeltaT),
        chargedPDF_svd1(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax),
        chargedPDF_svd2(range_mBC.vmin, range_mBC.vmax, range_dE.vmin, range_dE.vmax)
    {
        deltaT.setParameters(
                PAR::charged_dt_blifetime, PAR::charged_dt_Jc, PAR::charged_dt_Js1, PAR::charged_dt_Js2,
                PAR::charged_dt_fractionscale, useDeltaT?Event::dt_charged:-1);
    }

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            //Set Parameters for charged component
            chargedPDF_svd1.set(par[PAR::charged_svd1_ratio]);
            chargedPDF_svd1.fcn1.fcnx.set(par[PAR::charged_svd1_Mbc_mean], par[PAR::charged_svd1_Mbc_sigma]);
            chargedPDF_svd1.fcn1.fcny.set(par[PAR::charged_svd1_dE_mean], par[PAR::charged_svd1_dE_sigma]);
            chargedPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::charged_svd1_Mbc_argusC]);
            chargedPDF_svd1.fcn2.fcny.set(&par[PAR::charged_svd1_dE_cheb1]);

            return get_deltaT(e,par)*par[PAR::yield_svd1_charged]*chargedPDF_svd1(e.Mbc, e.dE);
        } else {
            //Set Parameters for charged component
            chargedPDF_svd2.set(par[PAR::charged_svd2_ratio]);
            chargedPDF_svd2.fcn1.fcnx.set(par[PAR::charged_svd2_Mbc_mean], par[PAR::charged_svd2_Mbc_sigma]);
            chargedPDF_svd2.fcn1.fcny.set(par[PAR::charged_svd2_dE_mean], par[PAR::charged_svd2_dE_sigma]);
            chargedPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::charged_svd2_Mbc_argusC]);
            chargedPDF_svd2.fcn2.fcny.set(&par[PAR::charged_svd2_dE_cheb1]);

            return get_deltaT(e,par)*par[PAR::yield_svd2_charged]*chargedPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH){
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::yield_svd1_charged];
        }
        if(svd & SVD2){
            yield += par[PAR::yield_svd2_charged];
        }
        return yield;
    }

    private:

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

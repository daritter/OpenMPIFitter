#ifndef MPIFitter_Mixed_h
#define MPIFitter_Mixed_h

#include <Functions.h>
#include "Event.h"
#include "Component.h"

namespace PAR {
    PARAM(yield_svd1_mixed);
    PARAM(mixed_svd1_ratio);
    PARAM(mixed_svd1_Mbc_mean);
    PARAM(mixed_svd1_Mbc_sigma);
    PARAM(mixed_svd1_Mbc_argusC);
    PARAM(mixed_svd1_dE_mean);
    PARAM(mixed_svd1_dE_sigma);
    PARAM(mixed_svd1_dE_cheb1);

    PARAM(yield_svd2_mixed);
    PARAM(mixed_svd2_ratio);
    PARAM(mixed_svd2_Mbc_mean);
    PARAM(mixed_svd2_Mbc_sigma);
    PARAM(mixed_svd2_Mbc_argusC);
    PARAM(mixed_svd2_dE_mean);
    PARAM(mixed_svd2_dE_sigma);
    PARAM(mixed_svd2_dE_cheb1);
};


class MixedPDF: public Component {
    public:
    MixedPDF(double lowerMbc, double upperMbc, double lowerdE, double upperdE):
        Component(),
        mixedPDF_svd1(lowerMbc, upperMbc, lowerdE, upperdE),
        mixedPDF_svd2(lowerMbc, upperMbc, lowerdE, upperdE)
    {}

    virtual double operator()(const DspDsmKsEvent& e, const std::vector<double> &par) {
        if(e.svdVs==0){
            //Set Parameters for mixed component
            mixedPDF_svd1.set(par[PAR::mixed_svd1_ratio]);
            mixedPDF_svd1.fcn1.fcnx.set(par[PAR::mixed_svd1_Mbc_mean], par[PAR::mixed_svd1_Mbc_sigma]);
            mixedPDF_svd1.fcn1.fcny.set(par[PAR::mixed_svd1_dE_mean], par[PAR::mixed_svd1_dE_sigma]);
            mixedPDF_svd1.fcn2.fcnx.set(e.benergy, par[PAR::mixed_svd1_Mbc_argusC]);
            mixedPDF_svd1.fcn2.fcny.set(&par[PAR::mixed_svd1_dE_cheb1]);

            return par[PAR::yield_svd1_mixed] * mixedPDF_svd1(e.Mbc, e.dE);
        }else{
            //Set Parameters for mixed component
            mixedPDF_svd2.set(par[PAR::mixed_svd2_ratio]);
            mixedPDF_svd2.fcn1.fcnx.set(par[PAR::mixed_svd2_Mbc_mean], par[PAR::mixed_svd2_Mbc_sigma]);
            mixedPDF_svd2.fcn1.fcny.set(par[PAR::mixed_svd2_dE_mean], par[PAR::mixed_svd2_dE_sigma]);
            mixedPDF_svd2.fcn2.fcnx.set(e.benergy, par[PAR::mixed_svd2_Mbc_argusC]);
            mixedPDF_svd2.fcn2.fcny.set(&par[PAR::mixed_svd2_dE_cheb1]);

            return par[PAR::yield_svd2_mixed] * mixedPDF_svd2(e.Mbc, e.dE);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH){
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::yield_svd1_mixed];
        }
        if(svd & SVD2){
            yield += par[PAR::yield_svd2_mixed];
        }
        return yield;
    }

    private:

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

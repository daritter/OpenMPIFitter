#ifndef MPIFitter_DeltaT_h
#define MPIFitter_DeltaT_h

#include <Functions.h>
#include "Event.h"

class DeltaTPDF {
    public:
    DeltaTPDF(Range range_dT, int isCharged=0):range_dT(range_dT), isCharged(isCharged), outlierPDF(range_dT.vmin, range_dT.vmax)
    {}

    void setParameters(int tau, int Jc, int Js1, int Js2, int cacheId=-1){
        this->tau = tau;
        this->Jc = Jc;
        this->Js1 = Js1;
        this->Js2 = Js2;
        this->cacheId=cacheId;
    }

    double operator()(const Event& e, const std::vector<double> &par) {
        double life_pdf(0), int_life_pdf(0), sin_pdf(0), cos_pdf(0);
        const Belle::dtres_param_t* const dtres_param = Belle::get_dtres_param( e.expNo, e.isMC );
        if(cacheId<0 || !e.dTcache[cacheId].get(tau,life_pdf, int_life_pdf, sin_pdf, cos_pdf)){
            //Calculate Lifetime components
            life_pdf = Belle::EfRkRdetRnp_fullrec(
                    e.deltaT, isCharged, par[tau], e.Ak, e.Ck,
                    e.vtx_ntrk, e.vtx_zerr, e.vtx_chi2, e.vtx_ndf,
                    e.tag_ntrk, e.tag_zerr, e.tag_chi2, e.tag_ndf, e.tag_isL,
                    dtres_param);

            int_life_pdf = Belle::norm_EfRkRdetRnp_fullrec(
                    range_dT.vmin, range_dT.vmax, isCharged, par[tau], e.Ak, e.Ck,
                    e.vtx_ntrk, e.vtx_zerr, e.vtx_chi2, e.vtx_ndf,
                    e.tag_ntrk, e.tag_zerr, e.tag_chi2, e.tag_ndf, e.tag_isL,
                    dtres_param);

            //Acp component
            cos_pdf = 0.5 / par[tau] * Belle::MfRkRdetRnp_fullrec(
                    e.deltaT, isCharged, par[tau], Event::deltaM, e.Ak, e.Ck,
                    e.vtx_ntrk, e.vtx_zerr, e.vtx_chi2, e.vtx_ndf,
                    e.tag_ntrk, e.tag_zerr, e.tag_chi2, e.tag_ndf, e.tag_isL,
                    dtres_param);

            //Scp component
            sin_pdf = 0.5 / par[tau] * Belle::AfRkRdetRnp_fullrec(
                    e.deltaT, isCharged, par[tau], Event::deltaM, e.Ak, e.Ck,
                    e.vtx_ntrk, e.vtx_zerr, e.vtx_chi2, e.vtx_ndf,
                    e.tag_ntrk, e.tag_zerr, e.tag_chi2, e.tag_ndf, e.tag_isL,
                    dtres_param);

            if(cacheId>=0) e.dTcache[cacheId].set(par[tau], life_pdf, int_life_pdf, sin_pdf, cos_pdf);
        }

        int eta = (e.m2DsmKs>e.m2DspKs)?1:-1;
        double sig_pdf = life_pdf*(1.0 - e.tag_q*e.wrongTag_dw)
            + e.tag_q*(1.0-2.0*e.wrongTag_w)*(eta*par[Jc]*cos_pdf - (par[Js1] + eta*par[Js2])*sin_pdf);

        outlierPDF.set(0.0,dtres_param->sig_ol);
        double fraction = (e.tag_ntrk>1 && e.vtx_ntrk>1)?dtres_param->fol_mul:dtres_param->fol_sgl;
        return (1-fraction)*sig_pdf + fraction*outlierPDF(e.deltaT);
    }

    private:
    Range range_dT;
    int tau;
    int Jc;
    int Js1;
    int Js2;
    int cacheId;
    int isCharged;

    /** PDF function components */
    Gauss outlierPDF;
};

#endif

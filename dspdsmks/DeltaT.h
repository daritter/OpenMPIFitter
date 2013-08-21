#ifndef MPIFitter_DeltaT_h
#define MPIFitter_DeltaT_h

#include <Functions.h>
#include "Event.h"

class DeltaTPDF {
    public:
        DeltaTPDF(Range range_dT, int isCharged):
            range_dT(range_dT), isCharged(isCharged), outlierPDF(range_dT.vmin, range_dT.vmax),
            offset(-1), eta_dependence(true),
            srec0(-1), srec1(-1), fol_sgl(-1), fol_mul(-1), lastExp(-1), lastMC(-1)
        {}

        void setParameters(int tau, int Jc, int Js1, int Js2,  int fractionScale, int outlierMean, int cacheId){
            this->tau = tau;
            this->Jc = Jc;
            this->Js1 = Js1;
            this->Js2 = Js2;
            this->fractionScale = fractionScale;
            this->outlierMean = outlierMean;
            this->cacheId=cacheId;
        }

        void setOffset(int offset){
            this->offset = offset;
        }

        void setCorrections(int srec0, int srec1, int fol_sgl, int fol_mul){
            this->srec0 = srec0;
            this->srec1 = srec1;
            this->fol_sgl = fol_sgl;
            this->fol_mul = fol_mul;
        }

        void setEtaDependence(bool eta){
            eta_dependence = eta;
        }

        double operator()(const Event& e, const std::vector<double> &par) {
            double life_pdf(0), int_life_pdf(0), sin_pdf(0), cos_pdf(0);
            //const Belle::dtres_param_t* const dtres_param = Belle::get_dtres_param( e.expNo, e.isMC );
            if(lastExp != e.expNo || lastMC != e.isMC){
                dtres_param = Belle::dtres_param_t(*Belle::get_dtres_param( e.expNo, e.isMC));
                lastExp = e.expNo;
                lastMC = e.isMC;
                if(srec0>=0) dtres_param.Srec[0] *= par[srec0+e.svdVs];
                if(srec1>=0) dtres_param.Srec[1] *= par[srec1+e.svdVs];
                if(fol_sgl>=0) dtres_param.fol_sgl *= par[fol_sgl+e.svdVs];
                if(fol_mul>=0) dtres_param.fol_mul *= par[fol_mul+e.svdVs];
            }
            const double abs_tau = par[tau];
            const double dT = e.deltaT + ((offset>=0)?par[offset]:0);
            if(cacheId<0 || !e.dTcache[cacheId].get(dT, abs_tau, life_pdf, int_life_pdf, sin_pdf, cos_pdf)){
                //Calculate Lifetime components
                life_pdf = Belle::EfRkRdetRnp_fullrec(
                        dT, isCharged, abs_tau, e.Ak, e.Ck,
                        e.vtx_ntrk, e.vtx_zerr, e.vtx_chi2, e.vtx_ndf,
                        e.tag_ntrk, e.tag_zerr, e.tag_chi2, (int)e.tag_ndf, e.tag_isL,
                        &dtres_param);

                int_life_pdf = Belle::norm_EfRkRdetRnp_fullrec(
                        range_dT.vmin, range_dT.vmax, isCharged, abs_tau, e.Ak, e.Ck,
                        e.vtx_ntrk, e.vtx_zerr, e.vtx_chi2, e.vtx_ndf,
                        e.tag_ntrk, e.tag_zerr, e.tag_chi2, (int)e.tag_ndf, e.tag_isL,
                        &dtres_param);

                //Acp component
                cos_pdf = 0.5 / abs_tau * Belle::MfRkRdetRnp_fullrec(
                        dT, isCharged, abs_tau, Event::deltaM, e.Ak, e.Ck,
                        e.vtx_ntrk, e.vtx_zerr, e.vtx_chi2, e.vtx_ndf,
                        e.tag_ntrk, e.tag_zerr, e.tag_chi2, (int)e.tag_ndf, e.tag_isL,
                        &dtres_param);

                //Scp component
                sin_pdf = 0.5 / abs_tau * Belle::AfRkRdetRnp_fullrec(
                        dT, isCharged, abs_tau, Event::deltaM, e.Ak, e.Ck,
                        e.vtx_ntrk, e.vtx_zerr, e.vtx_chi2, e.vtx_ndf,
                        e.tag_ntrk, e.tag_zerr, e.tag_chi2, (int)e.tag_ndf, e.tag_isL,
                        &dtres_param);

                if(cacheId>=0) e.dTcache[cacheId].set(dT, abs_tau, life_pdf, int_life_pdf, sin_pdf, cos_pdf);
            }

            double sig_pdf(0);
            if(eta_dependence) {
                sig_pdf = (
                        life_pdf * (1.0 - e.tag_q*e.wrongTag_dw) +
                        e.tag_q * (1.0-2.0*e.wrongTag_w) *
                        (e.eta * par[Jc] * cos_pdf - (par[Js1] + e.eta * par[Js2]) * sin_pdf)
                        )/int_life_pdf;

            }else{
                sig_pdf = (
                        life_pdf * (1.0 - e.tag_q*e.wrongTag_dw) +
                        e.tag_q * (1.0-2.0*e.wrongTag_w) *
                        (par[Jc] * cos_pdf + par[Js1] * sin_pdf)
                        )/int_life_pdf;
            }

            const double omean = (outlierMean>=0)?par[outlierMean]:0.0;
            outlierPDF.set(omean, dtres_param.sig_ol);
            double fraction = (e.tag_ntrk>1 && e.vtx_ntrk>1)?dtres_param.fol_mul:dtres_param.fol_sgl;
            if(fractionScale>=0) fraction*=par[fractionScale];

            return ((1-fraction)*sig_pdf + fraction*outlierPDF(dT))/(eta_dependence?4.0:2.0);
        }

    private:
        Range range_dT;
        int tau;
        int Jc;
        int Js1;
        int Js2;
        int fractionScale;
        int outlierMean;
        int cacheId;
        int isCharged;

        /** PDF function components */
        Gauss outlierPDF;

        int offset;
        bool eta_dependence;

        int srec0;
        int srec1;
        int fol_sgl;
        int fol_mul;
        int lastExp;
        int lastMC;

        Belle::dtres_param_t dtres_param;
};

#endif

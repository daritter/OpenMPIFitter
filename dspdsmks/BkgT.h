#ifndef MPIFitter_BkgT_h
#define MPIFitter_BkgT_h

#include <Functions.h>
#include "Event.h"

class BkgTPDF {
    public:
        BkgTPDF(Range range_dT, int charge):
            range_dT(range_dT), acp0(-1), promptPDF(range_dT.vmin, range_dT.vmax), outlierPDF(range_dT.vmin, range_dT.vmax)
        {
            s_main[1] = -1;
        }

        void setParameters(int sigma_main, int sigma_tail, int fraction_delta, int fraction_tail, bool multi=false) {
            int i = multi?1:0;
            this->s_main[i] = sigma_main;
            this->s_tail[i] = sigma_tail;
            this->f_delt[i] = fraction_delta;
            this->f_tail[i] = fraction_tail;
        }

        void setCommonParameters(int mu_delta, int tau, int mu_tau, int outlier_fraction=-1, int outlier_mean=-1, int outlier_scale=-1) {
            this->mu_delta = mu_delta;
            this->tau = tau;
            this->mu_tau = mu_tau;
            this->outlier_fraction = outlier_fraction;
            this->outlier_mean = outlier_mean;
            this->outlier_scale = outlier_scale;
        }

        void setAcp(int rbin0){
            this->acp0 = rbin0;
        }

        double Enp_conv_gauss(double weight, double dt, double tau, double mean, double sigma) const {
            if(weight==0.0) return 0;
            return weight * Belle::Enp_conv_gauss(dt, tau, tau, mean, sigma) /
                Belle::norm_Enp_conv_gauss(range_dT.vmin, range_dT.vmax, tau, tau, mean, sigma);
        }

        double operator()(const Event& e, const std::vector<double> &par) {
            const Belle::dtres_param_t* const dtres_param = Belle::get_dtres_param( e.expNo, e.isMC );
            const double abs_tau = par[tau];
            const double dT = e.deltaT;
            const double sigma = sqrt(e.tag_zerr*e.tag_zerr+e.vtx_zerr*e.vtx_zerr) * Belle :: dt_resol_global :: inv_bgc;
            const int i = (e.tag_ntrk>1 && e.vtx_ntrk>1 && s_main[1]>=0)?1:0;
            const double sigma_main = sigma * par[s_main[i]];
            const double sigma_tail = par[s_tail[i]];
            const double frac_tail = par[f_tail[i]];
            const double frac_delt = par[f_delt[i]];

            promptPDF.set((1-frac_tail), par[mu_delta], 0, sigma_main, sigma_tail);

            const double lifetime_part =
                Enp_conv_gauss(1-frac_tail, dT, abs_tau, par[mu_tau], sigma_main) +
                Enp_conv_gauss(frac_tail, dT, abs_tau, par[mu_tau], sigma_tail);

            const double pdf = ((acp0>=0)?(1+e.tag_q*par[acp0+e.rbin]):1.) * (frac_delt * promptPDF(dT) + (1-frac_delt)*lifetime_part);

            const double omean = (outlier_mean>=0)?par[outlier_mean]:0.0;
            outlierPDF.set(omean, dtres_param->sig_ol*(outlier_scale>=0?par[outlier_scale]:1));
            double fraction = (e.tag_ntrk>1 && e.vtx_ntrk>1)?dtres_param->fol_mul:dtres_param->fol_sgl;
            if(outlier_fraction>=0) fraction*=par[outlier_fraction];

            return ((1-fraction)*pdf + fraction*outlierPDF(dT))/4.0;
            //return Belle::AddOutlierWithBkg(expno, dT, 0, 0, pdf, e.vtx_ntrk, e.tag_ntrk, dtres_param, 1, 1, range_dT.vmin, range_dT.vmax, 1, 1);
        }

    private:
        Range range_dT;

        int s_main[2];
        int s_tail[2];
        int f_delt[2];
        int f_tail[2];
        int tau;
        int mu_tau;
        int mu_delta;
        int outlier_mean;
        int outlier_fraction;
        int outlier_scale;
        int acp0;

        /** PDF function components */
        DoubleGauss promptPDF;
        Gauss outlierPDF;
};

#endif

#ifndef MPIFitter_BkgT_h
#define MPIFitter_BkgT_h

#include <Functions.h>
#include "Event.h"

class BkgTPDF {
    public:
        BkgTPDF(Range range_dT, int charge):
            range_dT(range_dT), promptPDF(range_dT.vmin, range_dT.vmax), outlierPDF(range_dT.vmin, range_dT.vmax), eta_dependence(true)
        {
            s_main[0][1] = -1;
            s_main[1][1] = -1;
            acp0[0] = -1;
            acp0[1] = -1;
        }

        void setParameters(int svd, int sigma_main, int sigma_tail1, int sigma_tail2, int weight_tail1, int weight_tail2, int fraction_delta, bool multi=false) {
            int i = multi?1:0;
            s_main[svd][i] = sigma_main;
            s_tail1[svd][i] = sigma_tail1;
            s_tail2[svd][i] = sigma_tail2;
            w_tail1[svd][i] = weight_tail1;
            w_tail2[svd][i] = weight_tail2;
            f_delt[svd][i] = fraction_delta;
        }

        void setCommonParameters(int svd, int mu_delta, int tau, int mu_tau, int outlier_fraction=-1, int outlier_mean=-1, int outlier_scale=-1) {
            this->mu_delta[svd] = mu_delta;
            this->tau[svd] = tau;
            this->mu_tau[svd] = mu_tau;
            this->outlier_fraction[svd] = outlier_fraction;
            this->outlier_mean[svd] = outlier_mean;
            this->outlier_scale[svd] = outlier_scale;
        }

        void setAcp(int svd, int rbin0){
            this->acp0[svd] = rbin0;
        }

        void setEtaDependence(bool eta){
            eta_dependence = eta;
        }


        double Enp_conv_gauss(double weight, double dt, double tau, double mean, double sigma) const {
            if(weight==0.0) return 0;
            return weight * Belle::Enp_conv_gauss(dt, tau, tau, mean, sigma) /
                Belle::norm_Enp_conv_gauss(range_dT.vmin, range_dT.vmax, tau, tau, mean, sigma);
        }

        double operator()(const Event& e, const std::vector<double> &par) {
            const Belle::dtres_param_t* const dtres_param = Belle::get_dtres_param( e.expNo, e.isMC );
            const int svd = e.svdVs;
            const double abs_tau = par[tau[svd]];
            const double dT = e.deltaT;
            const double sigma = sqrt(e.tag_zerr*e.tag_zerr+e.vtx_zerr*e.vtx_zerr) * Belle::dt_resol_global::inv_bgc;
            const int i = (e.tag_ntrk>1 && e.vtx_ntrk>1 && s_main[svd][1]>=0)?1:0;
            const double sigma_main = sigma * par[s_main[svd][i]];
            const double sigma_tail1 = par[s_tail1[svd][i]];
            const double sigma_tail2 = par[s_tail2[svd][i]];
            const double weight_tail1 = par[w_tail1[svd][i]];
            const double weight_tail2 = par[w_tail2[svd][i]];
            const double frac_delt = par[f_delt[svd][i]];
            const double weight = 1 + weight_tail1 + weight_tail2;

            promptPDF.set(par[mu_delta[svd]], sigma_main, weight_tail1, 0, sigma_tail1, weight_tail2, 0, sigma_tail2);

            const double lifetime_part = (
                Enp_conv_gauss(           1, dT, abs_tau, par[mu_tau[svd]], sigma_main) +
                Enp_conv_gauss(weight_tail1, dT, abs_tau, par[mu_tau[svd]], sigma_main*sigma_tail1) +
                Enp_conv_gauss(weight_tail2, dT, abs_tau, par[mu_tau[svd]], sigma_main*sigma_tail1*sigma_tail2)
                ) / weight;

            const double pdf = ((acp0[svd]>=0)?(1+e.tag_q*par[acp0[svd]+e.rbin]):1.) * (frac_delt * promptPDF(dT) + (1-frac_delt)*lifetime_part);

            const double omean = (outlier_mean[svd]>=0)?par[outlier_mean[svd]]:0.0;
            outlierPDF.set(omean, dtres_param->sig_ol*(outlier_scale[svd]>=0?par[outlier_scale[svd]]:1));
            double fraction = (e.tag_ntrk>1 && e.vtx_ntrk>1)?dtres_param->fol_mul:dtres_param->fol_sgl;
            if(outlier_fraction[svd]>=0) fraction*=par[outlier_fraction[svd]];

            return ((1-fraction)*pdf + fraction*outlierPDF(dT))/(eta_dependence?4.0:2.0);
            //return Belle::AddOutlierWithBkg(expno, dT, 0, 0, pdf, e.vtx_ntrk, e.tag_ntrk, dtres_param, 1, 1, range_dT.vmin, range_dT.vmax, 1, 1);
        }

    private:
        Range range_dT;

        int s_main[2][2];
        int s_tail1[2][2];
        int s_tail2[2][2];
        int w_tail1[2][2];
        int w_tail2[2][2];
        int tau[2];
        int mu_tau[2];
        int mu_delta[2];
        int f_delt[2][2];
        int outlier_mean[2];
        int outlier_fraction[2];
        int outlier_scale[2];
        int acp0[2];

        /** PDF function components */
        MultiGauss<3> promptPDF;
        Gauss outlierPDF;

        bool eta_dependence;
};

#endif

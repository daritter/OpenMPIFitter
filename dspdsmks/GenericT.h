#ifndef MPIFitter_GenericT_h
#define MPIFitter_GenericT_h

#include <Functions.h>
#include "Event.h"

class GenericTPDF {
    public:
        GenericTPDF(Range range_dT, int isCharged=0):
            range_dT(range_dT), outlierPDF(range_dT.vmin, range_dT.vmax)
        {}

        void setParameters(bool multi, int mean1, int mean2, int mean3, int mean4, int sigma1, int sigma2, int sigma3, int sigma4, int weight2, int weight3, int weight4){
            int i = multi?1:0;
            this->mean1[i] = mean1;
            this->mean2[i] = mean2;
            this->mean3[i] = mean3;
            this->mean4[i] = mean4;
            this->sigma1[i] = sigma1;
            this->sigma2[i] = sigma2;
            this->sigma3[i] = sigma3;
            this->sigma4[i] = sigma4;
            this->weight2[i] = weight2;
            this->weight3[i] = weight3;
            this->weight4[i] = weight4;
        }

        void setCommonParameters(int tau, int fractionScale=-1, int outlierMean=-1, int outlierScale=-1){
            this->tau = tau;
            this->fractionScale = fractionScale;
            this->outlierMean = outlierMean;
            this->outlierScale = outlierScale;
        }

        double operator()(const Event& e, const std::vector<double> &par) {
            const Belle::dtres_param_t* const dtres_param = Belle::get_dtres_param( e.expNo, e.isMC );
            const double abs_tau = par[tau];
            const double sigma = sqrt(e.tag_zerr*e.tag_zerr+e.vtx_zerr*e.vtx_zerr);
            const int i = (e.tag_ntrk>1 && e.vtx_ntrk>1)?1:0;
            const double res_mean1 = par[mean1[i]];
            const double res_mean2 = res_mean1+par[mean2[i]];
            const double res_mean3 = res_mean2+par[mean3[i]];
            const double res_mean4 = res_mean3+par[mean4[i]];

            const double res_sigma1 = sigma*par[sigma1[i]]*Belle::dt_resol_global::inv_bgc;
            const double res_sigma2 = res_sigma1*par[sigma2[i]];
            const double res_sigma3 = res_sigma2*par[sigma3[i]];
            const double res_sigma4 = res_sigma3*par[sigma4[i]];

            const double res_weight2 = par[weight2[i]];
            const double res_weight3 = par[weight3[i]];
            const double res_weight4 = par[weight4[i]];
            const double res_sumw = 1.0/(1.0+res_weight2+res_weight3+res_weight4);

            //Calculate Lifetime components
            const double life_pdf = res_sumw * (
                            Belle::Ef_conv_gauss(e.deltaT, abs_tau, res_mean1, res_sigma1) +
                res_weight2*Belle::Ef_conv_gauss(e.deltaT, abs_tau, res_mean2, res_sigma2) +
                res_weight3*Belle::Ef_conv_gauss(e.deltaT, abs_tau, res_mean3, res_sigma3) +
                res_weight4*Belle::Ef_conv_gauss(e.deltaT, abs_tau, res_mean4, res_sigma4));

            const double int_life_pdf = res_sumw * (
                            Belle::norm_Ef_conv_gauss(range_dT.vmin, range_dT.vmax, abs_tau, res_mean1, res_sigma1) +
                res_weight2*Belle::norm_Ef_conv_gauss(range_dT.vmin, range_dT.vmax, abs_tau, res_mean2, res_sigma2) +
                res_weight3*Belle::norm_Ef_conv_gauss(range_dT.vmin, range_dT.vmax, abs_tau, res_mean3, res_sigma3) +
                res_weight4*Belle::norm_Ef_conv_gauss(range_dT.vmin, range_dT.vmax, abs_tau, res_mean4, res_sigma4));

            double omean = 0.0;
            double osig = dtres_param->sig_ol;
            if(outlierMean>=0){
                omean = par[outlierMean];
            }
            if(outlierScale>=0){
                osig *= par[outlierScale];
            }
            outlierPDF.set(omean, dtres_param->sig_ol);
            double fraction = (e.tag_ntrk>1 && e.vtx_ntrk>1)?dtres_param->fol_mul:dtres_param->fol_sgl;
            if(fractionScale>=0) fraction*=par[fractionScale];

            return (fraction*outlierPDF(e.deltaT) + (1.0-fraction)*(life_pdf/int_life_pdf))/4.0;
        }

    private:
        Range range_dT;
        int tau;
        int mean1[2];
        int mean2[2];
        int mean3[2];
        int mean4[2];
        int sigma1[2];
        int sigma2[2];
        int sigma3[2];
        int sigma4[2];
        int weight2[2];
        int weight3[2];
        int weight4[2];
        int fractionScale;
        int outlierMean;
        int outlierScale;

        Gauss outlierPDF;
};

#endif

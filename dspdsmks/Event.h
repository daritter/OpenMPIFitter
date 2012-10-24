#ifndef DDK_Event_h
#define DDK_Event_h

#include <string>
#include <limits>
#include <TTree.h>

#include <wtag.h>
#include "tatami/tatami.h"

class dTCache {
    public:
        dTCache():cache_lifetime(std::numeric_limits<double>::quiet_NaN()) {}

        bool get(double lifetime, double& life_pdf, double& int_life_pdf, double& cos_pdf, double& sin_pdf) const {
            if(cache_lifetime != lifetime) return false;
            lifetime = cache_lifetime;
            life_pdf = cache_life_pdf;
            int_life_pdf = cache_int_life_pdf;
            cos_pdf = cache_cos_pdf;
            sin_pdf = cache_sin_pdf;
            return true;
        }
        void set(double lifetime, double life_pdf, double int_life_pdf, double cos_pdf, double sin_pdf) {
            cache_lifetime = lifetime;
            cache_life_pdf = life_pdf;
            cache_int_life_pdf = int_life_pdf;
            cache_cos_pdf = cos_pdf;
            cache_sin_pdf = sin_pdf;
        }
        void reset() { cache_lifetime = std::numeric_limits<double>::quiet_NaN(); }

    protected:
        double cache_cos_pdf;
        double cache_sin_pdf;
        double cache_life_pdf;
        double cache_int_life_pdf;
        double cache_lifetime;
};

struct Event {
    static const double deltaM = 0.507;

    enum dtCachePosition {
        dt_signal  = 0,
        dt_mixed   = 1,
        dt_charged = 2
    };

    int expNo;
    int runNo;
    int evtNo;
    int svdVs;
    bool isMC;
    //double nB0;
    double benergy;
    int flag;
    double Mbc;
    double dE;
    //double channelP;
    //double channelM;
    int mcInfo;
    double m2DspKs;
    double m2DsmKs;
    double cosTheta;
    double deltaZ;
    //int q;
    //double r;

    int    vtx_ntrk;
    double vtx_zerr;
    double vtx_chi2;
    int    vtx_ndf;

    int    tag_ntrk;
    double tag_zerr;
    double tag_chi2;
    double tag_ndf;
    int    tag_q;
    double tag_r;
    int    tag_isL;

    mutable dTCache dTcache[3];

    int rbin;
    double wrongTag_w;
    double wrongTag_dw;
    double Ak;
    double Ck;
    double deltaT;

    void calculateValues(){
        deltaT = deltaZ*Belle::dt_resol_global::inv_bgc;
        tag_r = fabs(tag_r);
        rbin = Belle::set_rbin(tag_r);

        Ak = 0.0;
        Ck = 0.0;
        const double pb_cms_sq = (benergy*benergy) - (Mbc*Mbc);
        const double Eb_cms = sqrt((Belle::dt_resol_global::mbzero* Belle::dt_resol_global::mbzero) + pb_cms_sq);
        Belle::CalcAkCk( cosTheta, Eb_cms, &Ak, &Ck, Belle::dt_resol_global::mbzero );

        wrongTag_w = Belle::set_wtag( expNo, rbin, isMC );
        wrongTag_dw = Belle::set_dwtag( expNo, rbin, isMC );
    }



    void setBranches(TTree* tree, const std::string &bselection="bestLHsig"){
        std::string prefix("");
#define BADDRESS(x) tree->SetBranchAddress((prefix + #x).c_str(),&x)
        BADDRESS(expNo);
        BADDRESS(runNo);
        BADDRESS(evtNo);
        BADDRESS(svdVs);
        BADDRESS(isMC);
        BADDRESS(benergy);
        if(!bselection.empty()) prefix += bselection + ".";
        BADDRESS(flag);
        BADDRESS(Mbc);
        BADDRESS(dE);
        BADDRESS(mcInfo);
        BADDRESS(m2DspKs);
        BADDRESS(m2DsmKs);
        BADDRESS(cosTheta);
        BADDRESS(deltaZ);
        tree->SetBranchAddress((prefix + "tag.ntrk").c_str(),    &tag_ntrk);
        tree->SetBranchAddress((prefix + "tag.zerr").c_str(),    &tag_zerr);
        tree->SetBranchAddress((prefix + "tag.chi2").c_str(),    &tag_chi2);
        tree->SetBranchAddress((prefix + "tag.ndf").c_str(),     &tag_ndf);
        tree->SetBranchAddress((prefix + "tag.flavour").c_str(), &tag_q);
        tree->SetBranchAddress((prefix + "tag.qr").c_str(),      &tag_r);
        tree->SetBranchAddress((prefix + "tag.isL").c_str(),     &tag_isL);
        tree->SetBranchAddress((prefix + "vtx.ntrk").c_str(),    &vtx_ntrk);
        tree->SetBranchAddress((prefix + "vtx.zerr").c_str(),    &vtx_zerr);
        tree->SetBranchAddress((prefix + "vtx.chi2").c_str(),    &vtx_chi2);
        tree->SetBranchAddress((prefix + "vtx.ndf").c_str(),     &vtx_ndf);
    }
};


#endif

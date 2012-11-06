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
    static const double quality_cut = 50;
    static const double sigmaz_mult_cut = 0.02;
    static const double sigmaz_sngl_cut = 0.05;


    enum dtCachePosition {
        dt_signal,
        dt_misrecon,
        dt_mixed,
        dt_charged,
        MAX_DTCACHE
    };

    Event():
        expNo(0), svdVs(0), isMC(1), benergy(0),
        flag(0), Mbc(0), dE(0), m2DspKs(0), m2DsmKs(0), cosTheta(0), deltaZ(0),
        vtx_ntrk(0), vtx_zerr(0), vtx_chi2(0), vtx_ndf(0),
        tag_ntrk(0), tag_zerr(0), tag_chi2(0), tag_ndf(0), tag_q(1), tag_r(0), tag_isL(0),
        rbin(0), wrongTag_w(0), wrongTag_dw(0), Ak(0), Ck(0), deltaT(0), eta(0)
    {}

    int expNo;
    int svdVs;
    bool isMC;
    double benergy;
    int flag;
    double Mbc;
    double dE;
    double m2DspKs;
    double m2DsmKs;
    double cosTheta;
    double deltaZ;

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

    mutable dTCache dTcache[MAX_DTCACHE];

    //Calculated variables
    int rbin;
    double wrongTag_w;
    double wrongTag_dw;
    double Ak;
    double Ck;
    double deltaT;
    int eta;

    bool operator<(const Event& b) const {
        return benergy<b.benergy;
    }

    bool calculateValues(bool toyMC=false){
        if(toyMC){
            m2DspKs = 0;
            m2DsmKs = eta;
            deltaZ = deltaT / Belle::dt_resol_global::inv_bgc;
        }else{
            deltaT = deltaZ * Belle::dt_resol_global::inv_bgc;
            tag_r = fabs(tag_r);
            eta = (m2DsmKs>m2DspKs)?1:-1;
        }
        rbin = Belle::set_rbin(tag_r);

        //Check quality flag
        if(flag!=0 || tag_q==0) return false;
        //Random cuts imposed by ipcv for whatever reason
        if(vtx_ntrk > 1  && vtx_chi2/vtx_ndf > quality_cut) return false;
        if(tag_ntrk > 1  && tag_chi2/tag_ndf > quality_cut) return false;
        if(vtx_ntrk == 1 && vtx_zerr > sigmaz_sngl_cut) return false;
        if(vtx_ntrk > 1  && vtx_zerr > sigmaz_mult_cut) return false;
        if(tag_ntrk == 1 && tag_zerr > sigmaz_sngl_cut) return false;
        if(tag_ntrk > 1  && tag_zerr > sigmaz_mult_cut) return false;

        Ak = 0.0;
        Ck = 0.0;
        const double pb_cms_sq = (benergy*benergy) - (Mbc*Mbc);
        const double Eb_cms = sqrt((Belle::dt_resol_global::mbzero* Belle::dt_resol_global::mbzero) + pb_cms_sq);
        Belle::CalcAkCk( cosTheta, Eb_cms, &Ak, &Ck, Belle::dt_resol_global::mbzero );

        wrongTag_w = Belle::set_wtag( expNo, rbin, isMC );
        wrongTag_dw = Belle::set_dwtag( expNo, rbin, isMC );

        return true;
    }

    void setBranches(TTree* tree, const std::string &bselection){
        std::string prefix("");
#define BADDRESS__(x,var) tree->SetBranchAddress((prefix + x).c_str(),&var)
#define BADDRESS(x) BADDRESS__(#x,x)
#define BODDRESS(g,x) BADDRESS__(#g+"."+#x, g##_##x)
#define BUDDRESS(g,x,y) BADDRESS__(#g+"."+#x, g##_##y)
        BADDRESS(expNo);
        BADDRESS(svdVs);
        BADDRESS(isMC);
        BADDRESS(benergy);
        prefix += bselection + ".";
        BADDRESS(flag);
        BADDRESS(Mbc);
        BADDRESS(dE);
        BADDRESS(m2DspKs);
        BADDRESS(m2DsmKs);
        BADDRESS(cosTheta);
        BADDRESS(deltaZ);
        BODDRESS(tag,ntrk);
        BODDRESS(tag,zerr);
        BODDRESS(tag,chi2);
        BODDRESS(tag,ndf);
        BODDRESS(tag,isL);
        BUDDRESS(tag,flavour,q);
        BUDDRESS(tag,qr,r);
        BODDRESS(vtx,ntrk);
        BODDRESS(vtx,zerr);
        BODDRESS(vtx,chi2);
        BODDRESS(vtx,ndf);
    }

    void createBranches(TTree* tree, const std::string &bselection){
        std::string prefix("");
#define ADDBRANCH__(name,var,type)   tree->Branch((prefix + name).c_str(),&var,(prefix+name+"/"+#type).c_str())
#define ADDBRANCH(x,type)            ADDBRANCH__(#x,x,type)
#define ADDBRONCH(g,x,type)          ADDBRANCH__(#g+"."+#x,g##_##x,type)
#define ADDBRUNCH(g,x,y,type)        ADDBRANCH__(#g+"."+#x,g##_##y,type)
        ADDBRANCH(expNo,I);
        ADDBRANCH(svdVs,I);
        ADDBRANCH(isMC,O);
        ADDBRANCH(benergy,D);
        prefix += bselection + ".";
        ADDBRANCH(flag,I);
        ADDBRANCH(Mbc,D);
        ADDBRANCH(dE,D);
        ADDBRANCH(m2DspKs,D);
        ADDBRANCH(m2DsmKs,D);
        ADDBRANCH(cosTheta,D);
        ADDBRANCH(deltaZ,D);
        ADDBRONCH(tag,ntrk,I);
        ADDBRONCH(tag,zerr,D);
        ADDBRONCH(tag,chi2,D);
        ADDBRONCH(tag,ndf,D);
        ADDBRONCH(tag,isL,D);
        ADDBRUNCH(tag,flavour,q,I);
        ADDBRUNCH(tag,qr,r,D);
        ADDBRONCH(vtx,ntrk,I);
        ADDBRONCH(vtx,zerr,D);
        ADDBRONCH(vtx,chi2,D);
        ADDBRONCH(vtx,ndf,I);
    }

    void reset(){
        for(int i=0; i<MAX_DTCACHE; ++i) dTcache[i].reset();
    }
};


#endif

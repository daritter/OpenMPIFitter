#ifndef DDK_Event_h
#define DDK_Event_h

#include <string>
#include <limits>
#include <cassert>
#include <cmath>
#include <TTree.h>
#include "Range.h"

#include <wtag.h>

//Hacky dT inclusion to get rid of warning
//Include all dependencies of dt_resolution.h first
#include "tatami/libcnvl.h"
#include "tatami/dt_resolution_const.h"
//define typedef as empty token to get rid of the "typedef ignored" warning
#define typedef
//open Belle namespace as closing of namespace is out of include-guard in
//dt_resolution_const.h. as this is included first in dt_resolution.h this
//works. How wonderful
namespace Belle {
#include "tatami/dt_resolution.h"
//restore typedef and include the rest
#undef typedef
#include "tatami/tatami.h"

class dTCache {
    public:
        dTCache():cache_lifetime(std::numeric_limits<double>::quiet_NaN()) {}

        bool get(double dT, double lifetime, double& life_pdf, double& int_life_pdf, double& cos_pdf, double& sin_pdf) const {
            ++requests;
            if(cache_lifetime != lifetime || cache_dT != dT) return false;
            lifetime = cache_lifetime;
            life_pdf = cache_life_pdf;
            int_life_pdf = cache_int_life_pdf;
            cos_pdf = cache_cos_pdf;
            sin_pdf = cache_sin_pdf;
            ++hits;
            return true;
        }
        void set(double dT, double lifetime, double life_pdf, double int_life_pdf, double cos_pdf, double sin_pdf) {
            ++sets;
            cache_dT = dT;
            cache_lifetime = lifetime;
            cache_life_pdf = life_pdf;
            cache_int_life_pdf = int_life_pdf;
            cache_cos_pdf = cos_pdf;
            cache_sin_pdf = sin_pdf;
        }
        void reset() { cache_lifetime = std::numeric_limits<double>::quiet_NaN(); }

        static void print_stats();

    protected:
        double cache_dT;
        double cache_lifetime;
        double cache_cos_pdf;
        double cache_sin_pdf;
        double cache_life_pdf;
        double cache_int_life_pdf;

        static size_t requests;
        static size_t hits;
        static size_t sets;
};

struct Event {
    constexpr static double deltaM = 0.507;
    constexpr static double quality_cut = 50;
    constexpr static double sigmaz_mult_cut = 0.02;
    constexpr static double sigmaz_sngl_cut = 0.05;


    enum dtCachePosition {
        dt_signal,
        dt_misrecon,
        MAX_DTCACHE
    };

    Event():
        expNo(0), svdVs(0), isMC(1), benergy(0),
        flag(0), Mbc(0), dE(0), m2DspKs(0), m2DsmKs(0), cosTheta(0), deltaZ(0),
        vtx_ntrk(0), vtx_zerr(0), vtx_chi2(0), vtx_ndf(0),
        tag_ntrk(0), tag_zerr(0), tag_chi2(0), tag_ndf(0), tag_q(1), tag_r(0), tag_isL(0),
        rbin(-1), wrongTag_w(0), wrongTag_dw(0), Ak(0), Ck(0), deltaT(0), eta(0)
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

    //Calculated variables
    int rbin;
    double wrongTag_w;
    double wrongTag_dw;
    double Ak;
    double Ck;
    double deltaT;
    int eta;

    mutable dTCache dTcache[MAX_DTCACHE];

    bool operator<(const Event& b) const {
        return benergy<b.benergy;
    }

    bool calculateValues(bool qualityCuts=true, bool toyMC=false){
        if(toyMC){
            m2DspKs = 0;
            m2DsmKs = eta;
            tag_r = Belle::set_r(rbin);
            assert(Belle::set_rbin(tag_r) == rbin);
            deltaZ = deltaT / Belle::dt_resol_global::inv_bgc;
        }else{
            deltaT = deltaZ * Belle::dt_resol_global::inv_bgc;
            tag_r = fabs(tag_r);
            eta = (m2DsmKs>m2DspKs)?1:-1;
        }

        //Check quality flag
        if(flag!=0 || tag_q==0) return false;

        rbin = Belle::set_rbin(tag_r);
        assert(rbin>=0);
        assert(rbin<7);

        if(qualityCuts){
            //Random cuts imposed by ipcv for whatever reason
            if(vtx_ntrk > 1  && vtx_chi2/vtx_ndf > quality_cut) return false;
            if(tag_ntrk > 1  && tag_chi2/tag_ndf > quality_cut) return false;
            if(vtx_ntrk == 1 && vtx_zerr > sigmaz_sngl_cut) return false;
            if(vtx_ntrk > 1  && vtx_zerr > sigmaz_mult_cut) return false;
            if(tag_ntrk == 1 && tag_zerr > sigmaz_sngl_cut) return false;
            if(tag_ntrk > 1  && tag_zerr > sigmaz_mult_cut) return false;
        }

        Ak = 0.0;
        Ck = 0.0;
        const double pb_cms_sq = (benergy*benergy) - (Mbc*Mbc);
        const double Eb_cms = sqrt((Belle::dt_resol_global::mbzero* Belle::dt_resol_global::mbzero) + pb_cms_sq);
        Belle::CalcAkCk( cosTheta, Eb_cms, &Ak, &Ck, Belle::dt_resol_global::mbzero );

        wrongTag_w = Belle::set_wtag( expNo, rbin, isMC );
        wrongTag_dw = Belle::set_dwtag( expNo, rbin, isMC );

        return true;
    }

    void setBranches(TTree* tree, const std::string &bselection);
    void createBranches(TTree* tree, const std::string &bselection);
    void reset(){
        for(int i=0; i<MAX_DTCACHE; ++i) dTcache[i].reset();
    }
};


#endif

#ifndef DDK_Event_h
#define DDK_Event_h

#include <string>
#include <limits>

#include <TTree.h>

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

struct DspDsmKsEvent {
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
    int q;
    double r;

    //mutable dTCache dTsignal;

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
        tree->SetBranchAddress((prefix + "tag.flavour").c_str(), &q);
        tree->SetBranchAddress((prefix + "tag.qr").c_str(), &r);
    }
};


#endif

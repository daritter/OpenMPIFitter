#ifndef MPIFitter_Component_h
#define MPIFitter_Component_h

#include "Event.h"
#include "Range.h"
#include "DeltaT.h"

class Component {
    public:
    enum EnabledSVD {
        SVD1 = 1<<0,
        SVD2 = 1<<1,
        BOTH = SVD1 | SVD2
    };

    Component() {}

    virtual ~Component(){}
    virtual double operator()(const Event& e, const std::vector<double> &par, bool normalized=false) = 0;
    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const = 0;
    virtual double get_deltaT(const Event& e, const std::vector<double> &par, bool anyway=false) = 0;
    virtual double get_cosTheta(const Event &e) = 0;

    virtual double get_fraction(const std::vector<double> &par, EnabledSVD svd=BOTH) {
        return 1;
    }

    virtual void get_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, double scale=1.0) const = 0;

    static double get_rbinFraction(int rbin, int rbin1, const std::vector<double> &par, int corr_rbin1=-1) {
        if(rbin<0) return 1.0;
        if(rbin>6) return 0.0;
        double sum(0);
        for(int i=0; i<7; ++i) sum += par[rbin1+i];
        if(corr_rbin1>=0){
            if(rbin<6) return par[rbin1+rbin]/sum * par[corr_rbin1+rbin];
            double corr_rbin6(1);
            for(int i=0; i<6; ++i) corr_rbin6 -= par[rbin1+i]/sum * par[corr_rbin1+i];
            return std::max(0., corr_rbin6);
        }
        return par[rbin1+rbin]/sum;
    }

    protected:

    void fill_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, int svd1_rbin1, int svd2_rbin1, double scale, int svd1_cbin1=-1, int svd2_cbin1=-1) const {
        if((svd & SVD1) && (svd & SVD2)){
            std::cerr << "Cannot give rbin fractions for both svds";
        }
        const int rbin1 = (svd & SVD1)?svd1_rbin1:svd2_rbin1;
        const int cbin1 = (svd & SVD1)?svd1_cbin1:svd2_cbin1;
        for(int i=0; i<7; ++i){
            fractions[i] = get_rbinFraction(i, rbin1, par, cbin1)*scale;
        }
    }
};

template<class T=DeltaTPDF> class DeltaTComponent: public Component {
    public:

    DeltaTComponent(Range range_dT, bool isCharged, bool useDeltaT, bool flatCosTheta=true):useDeltaT(useDeltaT), flatCosTheta(flatCosTheta), deltaT(range_dT, isCharged?1:0)
    {}

    virtual ~DeltaTComponent(){}

    virtual double get_deltaT(const Event& e, const std::vector<double> &par, bool anyway=false){
        if(!useDeltaT && !anyway) return 1.0;
        return deltaT(e, par);
    }

    virtual double get_cosTheta(const Event &e){
        return flatCosTheta?1.0:(1.0-e.cosTheta*e.cosTheta);
    }

    protected:
    bool useDeltaT;
    bool flatCosTheta;
    T deltaT;
};

#endif

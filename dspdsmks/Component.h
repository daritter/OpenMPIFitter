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
    virtual double operator()(const Event& e, const std::vector<double> &par) = 0;
    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const = 0;
    virtual double get_deltaT(const Event& e, const std::vector<double> &par, bool anyway=false) = 0;
    virtual double get_cosTheta(const Event &e) = 0;

    virtual void get_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, double scale=1.0) = 0;

    protected:

    void fill_rbinFractions(const std::vector<double> &par, std::vector<double> &fractions, EnabledSVD svd, int svd1_rbin1, int svd2_rbin1, double scale){
        if((svd & SVD1) && (svd & SVD2)){
            std::cerr << "Cannot give rbin fractions for both svds";
        }
        const int rbin1 = (svd & SVD1)?svd1_rbin1:svd2_rbin1;
        double sum(0);
        for(int i=0; i<7; ++i) {
            sum += par[rbin1+i];
        }
        for(int i=0; i<7; ++i) {
            fractions[i] = par[rbin1+i]/sum*scale;
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

    static double get_rbinFraction(int rbin, int rbin1, const std::vector<double> &par) {
        if(rbin<0) return 1.0;
        if(rbin>6) return 0.0;
        double sum(0);
        for(int i=0; i<7; ++i) sum += par[rbin1+i];
        return par[rbin1+rbin]/sum;
    }

    protected:
    bool useDeltaT;
    bool flatCosTheta;
    T deltaT;
};

#endif

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
};

template<class T=DeltaTPDF> class DeltaTComponent: public Component {
    public:

    DeltaTComponent(Range range_dT, bool isCharged, bool useDeltaT, bool flatCosTheta=true):useDeltaT(useDeltaT), flatCosTheta(flatCosTheta), deltaT(range_dT, isCharged?1:0)
    {}

    virtual ~DeltaTComponent(){}
    virtual double operator()(const Event& e, const std::vector<double> &par) = 0;
    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH, int rbin=-1) const = 0;

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
        if(rbin==6){
            double fraction(1.);
            for(int i=0; i<6; ++i) fraction -= par[rbin1+i];
            //if(fraction<0){
            //    std::cerr << "rbin Fraction for bin 7 out of bounds: " << fraction << std::endl;
            //    std::abort();
            //}
            //std::cout << rbin1 << ", " << "6: " << fraction << std::endl;
            return fraction;
        }
        //std::cout << rbin1 << ", " << rbin << ": " << par[rbin1+rbin] << std::endl;
        return par[rbin1+rbin];
    }

    protected:
    bool useDeltaT;
    bool flatCosTheta;
    T deltaT;
};

#endif

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
    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH) = 0;
    virtual double get_deltaT(const Event& e, const std::vector<double> &par, bool anyway=false) = 0;
    virtual double get_cosTheta(const Event &e) = 0;
};

template<class T=DeltaTPDF> class DeltaTComponent: public Component {
    public:

    DeltaTComponent(Range range_dT, bool isCharged, bool useDeltaT, bool flatCosTheta=true):useDeltaT(useDeltaT), flatCosTheta(flatCosTheta), deltaT(range_dT, isCharged?1:0)
    {}

    virtual ~DeltaTComponent(){}
    virtual double operator()(const Event& e, const std::vector<double> &par) = 0;
    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH) = 0;

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

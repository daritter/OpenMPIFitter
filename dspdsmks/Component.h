#ifndef MPIFitter_Component_h
#define MPIFitter_Component_h


class Component {
    public:
    enum EnabledSVD {
        SVD1 = 1<<0,
        SVD2 = 1<<1,
        BOTH = SVD1 | SVD2
    };

    Component()    {}

    virtual double operator()(const DspDsmKsEvent& e, const std::vector<double> &par) = 0;
    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH) = 0;
};

#endif

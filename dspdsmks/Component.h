#ifndef MPIFitter_Component_h
#define MPIFitter_Component_h


class Component {
    public:
    Component()    {}

    virtual double operator()(const DspDsmKsEvent& e, const std::vector<double> &par) = 0;
};

#endif

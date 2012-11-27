
class DummyPDF: public DeltaTComponent<> {
    public:
    DummyPDF(Range range_mBC, Range range_dE, Range range_dT, bool useDeltaT=false):
        DeltaTComponent<>(range_dT, false, false), range_mBC(range_mBC),
        dummyPDF(range_mBC.vmin, range_mBC.vmax)
    {}
    virtual ~DummyPDF(){}

    virtual double operator()(const Event& e, const std::vector<double> &par) {
        if(e.svdVs == 0){
            dummyPDF.set_limits(range_mBC.vmin,std::min(e.benergy,(double) range_mBC.vmax));
            dummyPDF.set(par[PAR::misrecon_svd1_ratio]);
            dummyPDF.fcn1.set(&par[PAR::misrecon_svd1_Mbc_mean]);
            dummyPDF.fcn2.set(e.benergy, par[PAR::misrecon_svd1_Mbc_argusC]);
            return get_yield(par,Component::SVD1) * dummyPDF(e.Mbc);
        }else{
            dummyPDF.set_limits(range_mBC.vmin,std::min(e.benergy,(double) range_mBC.vmax));
            dummyPDF.set(par[PAR::misrecon_svd2_ratio]);
            dummyPDF.fcn1.set(&par[PAR::misrecon_svd2_Mbc_mean]);
            dummyPDF.fcn2.set(e.benergy, par[PAR::misrecon_svd2_Mbc_argusC]);
            return get_yield(par,Component::SVD2) * dummyPDF(e.Mbc);
        }
    }

    virtual double get_yield(const std::vector<double> &par, EnabledSVD svd=BOTH){
        double yield(0);
        if(svd & SVD1){
            yield += par[PAR::ratio_misrecon_svd1]*par[PAR::yield_signal_svd1];
        }
        if(svd & SVD2){
            yield += par[PAR::ratio_misrecon_svd2]*par[PAR::yield_signal_svd2];
        }
        return yield;
    }

    private:

    Range range_mBC;

    /** PDF function components */
    //Gauss dummyPDF;
    Add1DFcn<MultiGauss<1>, Argus> dummyPDF;
};


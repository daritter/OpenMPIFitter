#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <Parameters.h>
#include "../dspdsmks/DspDsmKs.h"
#include <fstream>

struct ParameterHelper {
    static const std::string tostring(Parameter &self){
        return (boost::format("<%1% at %2% value=%3% error=%4%>") % self.name() % &self % self.value() % self.error()).str();
    }
};

struct RangeHelper {
    static const std::string tostring(Range &self){
        return (boost::format("<Range vmin=%1% vmax=%2%>") % self.vmin % self.vmax).str();
    }
};

struct PDFHelper {
    static void set_files(DspDsmKsPDF &self, boost::python::list& ns){
        std::vector<std::string> &filenames = self.getFiles();
        filenames.clear();
        for (int i = 0; i < len(ns); ++i) {
            filenames.push_back(boost::python::extract<std::string>(ns[i]));
        }
    }
    static void print_yields(DspDsmKsPDF &self, const Parameters &p){
        const std::vector<double> par = p.getValues();
        std::string names[] = {"signal","misrecon","bbar","continuum"};
        int components[] = {DspDsmKsPDF::CMP_signal, DspDsmKsPDF::CMP_misrecon, DspDsmKsPDF::CMP_bbar, DspDsmKsPDF::CMP_continuum};
        for(int i=0; i<4; ++i){
                int cmp = components[i];
                self.setComponents((DspDsmKsPDF::EnabledComponents)cmp);
                std::string name = names[i];
                std::cout << name << ": SVD1 = " << self.get_yield(par, Component::SVD1) << ", SVD2 = " << self.get_yield(par, Component::SVD2) << std::endl;
        }
    }
};

bool operator==(const Event &a, const Event &b){ return false;}
bool operator!=(const Event &a, const Event &b){ return true;}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(LoadParamsFromFile, Parameters::load, 1, 4);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(LoadEvents, DspDsmKsPDF::load, 0, 3);

BOOST_PYTHON_MODULE(dspdsmks)
{
    using namespace boost::python;

    double (Parameter::*value_get)() const = &Parameter::value;
    void   (Parameter::*value_set)(double) = &Parameter::value;
    double (Parameter::*error_get)() const = &Parameter::error;
    void   (Parameter::*error_set)(double) = &Parameter::error;
    class_<Parameter>("Parameter", no_init)
        .add_property("name", make_function(&Parameter::name, return_value_policy<copy_const_reference>()))
        .add_property("value", value_get, value_set)
        .add_property("error", error_get, error_set)
        .add_property("fixed", &Parameter::isFixed)
        .add_property("min", &Parameter::min)
        .add_property("max", &Parameter::max)
        .add_property("has_min", &Parameter::hasMin)
        .add_property("has_max", &Parameter::hasMax)
        .def("__repr__", &ParameterHelper::tostring)
        ;

    Parameter& (Parameters::*get_id)(int) = &Parameters::operator[];
    Parameter& (Parameters::*get_name)(const std::string&) = &Parameters::operator[];
    class_<Parameters>("Parameters")
        .def("__len__", &Parameters::size)
        .def("__getitem__", get_id, return_value_policy<reference_existing_object>())
        .def("__call__", get_name, return_value_policy<reference_existing_object>())
        .def("__iter__", iterator<Parameters>())
        .def("load", &Parameters::load, LoadParamsFromFile())
        .def("save", &Parameters::save)
        .def("print", &Parameters::print)
        ;

    class_<Range>("Range", init<double, double>())
        .def("__call__", &Range::operator())
        .def_readwrite("vmin",&Range::vmin)
        .def_readwrite("vmax",&Range::vmax)
        .def("__repr__", &RangeHelper::tostring)
        ;

    class_<Event>("Event")
        .def_readonly("expNo", &Event::expNo)
        .def_readonly("svdVs", &Event::svdVs)
        .def_readonly("isMC", &Event::isMC)
        .def_readonly("benergy", &Event::benergy)
        .def_readonly("Mbc", &Event::Mbc)
        .def_readonly("dE", &Event::dE)
        .def_readonly("m2DspKs", &Event::m2DspKs)
        .def_readonly("m2DsmKs", &Event::m2DsmKs)
        .def_readonly("cosTheta", &Event::cosTheta)
        .def_readonly("deltaZ", &Event::deltaZ)
        .def_readonly("vtx_ntrk", &Event::vtx_ntrk)
        .def_readonly("vtx_zerr", &Event::vtx_zerr)
        .def_readonly("vtx_chi2", &Event::vtx_chi2)
        .def_readonly("vtx_ndf", &Event::vtx_ndf)
        .def_readonly("tag_ntrk", &Event::tag_ntrk)
        .def_readonly("tag_zerr", &Event::tag_zerr)
        .def_readonly("tag_chi2", &Event::tag_chi2)
        .def_readonly("tag_ndf", &Event::tag_ndf)
        .def_readonly("tag_q", &Event::tag_q)
        .def_readonly("tag_r", &Event::tag_r)
        .def_readonly("tag_isL", &Event::tag_isL)
        .def_readonly("rbin", &Event::rbin)
        .def_readonly("deltaT", &Event::deltaT)
        .def_readonly("eta", &Event::eta)
        .def_readonly("wrongTag_w", &Event::wrongTag_w)
        .def_readonly("wrongTag_dw", &Event::wrongTag_dw)
        .def_readonly("Ak", &Event::Ak)
        .def_readonly("Ck", &Event::Ck)
        ;

    class_< std::vector<Event> >("EventList")
        .def(vector_indexing_suite< std::vector<Event> >())
        ;

    Range& (DspDsmKsPDF::*range_mBC)() = &DspDsmKsPDF::getRange_mBC;
    Range& (DspDsmKsPDF::*range_dE)() = &DspDsmKsPDF::getRange_dE;
    Range& (DspDsmKsPDF::*range_dT)() = &DspDsmKsPDF::getRange_dT;
    class_<DspDsmKsPDF>("DspDsmKsPDF")
        .def("size", &DspDsmKsPDF::size)
        .def("load", &DspDsmKsPDF::load, LoadEvents())
        .def("set_filenames", &PDFHelper::set_files)
        .def("range_mBC", range_mBC, return_value_policy<reference_existing_object>())
        .def("range_dE", range_dE, return_value_policy<reference_existing_object>())
        .def("range_dT", range_dT, return_value_policy<reference_existing_object>())
        .def("data", &DspDsmKsPDF::getData, return_value_policy<reference_existing_object>())
        .def("yields", &PDFHelper::print_yields)
        ;
}

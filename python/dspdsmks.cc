#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>
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
        return (boost::format("<Range vmin=%2% vmax=%3%>") % self.vmin % self.vmax).str();
    }
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(LoadParamsFromFile, Parameters::load, 1, 4);

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
}

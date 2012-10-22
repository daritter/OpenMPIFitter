#include "Parameters.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

using namespace std;
using namespace boost;

Parameters::Parameters(){
    const ParameterList::list_type &parameterList = ParameterList::getParameterList();
    m_parameters.reserve(parameterList.size());
    BOOST_FOREACH(const std::string &name, parameterList){
        m_parameters.push_back(Parameter(name));
    }
}

void Parameters::load(istream &in){
    int lineNr=0;
    while(!in.eof()){
        ++lineNr;
        string line;
        getline(in,line);
        size_t comment_pos = line.find('#');
        if(comment_pos!=string::npos){
            line = line.substr(0,comment_pos);
        }
        trim(line);
        if(line=="") continue;

        string name;
        stringstream buffer(line);
        buffer >> name;
        try{
            int index = ParameterList::getIndex(name);
            m_parameters[index].load(buffer);
        }catch(bad_lexical_cast &e){
            cerr << "ERROR reading parameter '" << name << "' (line " << lineNr << "): " << e.what() << endl;
            throw e;
        }catch(invalid_argument &e){
            cerr << "WARNING: Unknown parameter '" << name << "' (line " << lineNr << ")\n";
        }
    }
}

void Parameters::save(ostream &out) const {
    out << boost::format("#Name %|32t| %=17s %=17s %=17s %=17s %5s\n")
        % "value" % "error" % "min" % "max" % "fixed";
    BOOST_FOREACH(const Parameter &p, m_parameters){
        out << p;
    }
}

std::vector<double> Parameters::getValues() const {
    std::vector<double> result(m_parameters.size());
    for(size_t i=0; i<m_parameters.size(); ++i){
        result[i] = m_parameters[i].value();
    }
    return result;
}

void Parameters::fixParameters(const std::string& fixParameters){
    boost::regex fixed(fixParameters);
    BOOST_FOREACH(Parameter& p, m_parameters){
        p.setDynamicFix(boost::regex_match(p.name(),fixed));
    }
}

ROOT::Minuit2::MnUserParameters Parameters::getMnParams() const {
    ROOT::Minuit2::MnUserParameters mnParams;

    BOOST_FOREACH(const Parameter& p, m_parameters){
        if(!mnParams.Add(p.name(), p.value(), p.error())){
            std::cerr << "ERROR adding parameter" << std::endl;
        }
        int index = mnParams.Index(p.name());
        assert(index == ParameterList::getIndex(p.name()));
        if(p.isFixed()) mnParams.Fix(index);
        else if(p.hasMin() && p.hasMax()) mnParams.SetLimits(index, p.min(), p.max());
        else if(p.hasMin()) mnParams.SetLowerLimit(index, p.min());
        else if(p.hasMax()) mnParams.SetUpperLimit(index, p.max());
    }

    return mnParams;
}

void Parameters::update(ROOT::Minuit2::MnUserParameters mnParams){
    BOOST_FOREACH(const ROOT::Minuit2::MinuitParameter &mnp, mnParams.Parameters()){
        int index = ParameterList::getIndex(mnp.GetName());
        Parameter &p = m_parameters[index];
        if(!p.isFixed()) {
            p.value(mnp.Value());
            p.error(mnp.Error());
        }
    }
}

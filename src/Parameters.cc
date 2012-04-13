#include "Parameters.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

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

//MnUserParameters Parameters::getMnParams() const {
    //MnUserParameters mnParams;

    //BOOST_FOREACH(const Parameter& p, m_parameters){
        //mnParams.Add(p.name(), p.value(), p.error());
        //if(p.hasMin()) mnParams.SetLowerLimit(p.name(), p.min());
        //if(p.hasMax()) mnParams.SetUpperLimit(p.name(), p.max());
        //if(p.isFixed()) mnParams.Fix(p.name());
    //}

    //return mnParams;
//}

//Parameters::update(MnUserParameters mnParams){
    //BOOST_FOREACH(const MinuitParameter &mnp, mnParams.Parameters()){
        //int index = ParameterList::getIndex(mnp.GetName());
        //Parameter &p = m_parameters[index];
        //p.value(mnp.Value());
        //p.error(mnp.Error());
    //}
//}

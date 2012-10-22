#include "Parameter.h"

#include <limits>
#include <stdexcept>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

namespace {
    inline void readValue(istream &in, double &val){
        if(in.eof()) return;
        string tmp;
        in >> tmp;
        to_lower(tmp);
        if(tmp == "inf")       val =  numeric_limits<double>::infinity();
        else if(tmp == "-inf") val = -numeric_limits<double>::infinity();
        else if(tmp == "nan")  val =  numeric_limits<double>::quiet_NaN();
        else                   val =  lexical_cast<double>(tmp);
    }

    inline void readBool(istream &in, bool &val){
        if(in.eof()) return;
        string tmp;
        in >> tmp;
        trim_if(tmp,is_any_of("*"));
        if(tmp == "" || tmp=="0" || tmp=="false" || tmp=="no" || tmp=="N" || tmp=="n"){
            val = false;
        }else{
            val = true;
        }
    }
}

void Parameter::load(istream& in){
    readValue(in,m_value);
    readValue(in,m_error);
    readValue(in,m_min);
    readValue(in,m_max);
    readBool(in,m_fixed);
}

void Parameter::save(ostream& out) const {
    string fixed=m_dynfix?"*":"";
    fixed += m_fixed?"Y":"N";
    out << format("%s %|32t| %17.10e %17.10e %17.10e %17.10e %5s\n")
        % m_name % m_value % m_error % m_min % m_max % fixed;
}

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
        //Read floating point from stream and convert nan and inf accordingly
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
        //Read a bool value from stream
        if(in.eof()) return;
        string tmp;
        in >> tmp;
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

bool Parameter::save(ostream& out, bool istty) const {
    if(istty && isFixed()) return false;
    bool end=false;
    if(istty){
        double significance = std::fabs(m_value / m_error);
        bool loose = (significance <= 1.0);
        if(m_changed && loose) out << ANSI_PURPLE;
        else if(loose) out << ANSI_RED;
        else if(m_changed) out << ANSI_GREEN;
        end = loose | m_changed;
        out << format("%s %|32t| %17.10e %17.10e %6.6g %6.6g %12.2f")
            % m_name % m_value % m_error % m_min % m_max % significance;
    }else{
        out << format("%s %|32t| %17.10e %17.10e %6.6g %6.6g %5s")
            % m_name % m_value % m_error % m_min % m_max % (m_fixed?"Y":"N");
    }

    if(end) out << ANSI_END;
    return true;
}

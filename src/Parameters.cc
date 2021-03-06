#include "Parameters.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <fstream>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>

using namespace std;
using namespace boost;

Parameters::Parameters(){
    const ParameterList::list_type &parameterList = ParameterList::getParameterList();
    m_parameters.reserve(parameterList.size());
    for(const std::string &name: parameterList){
        m_parameters.push_back(Parameter(name));
    }
}

bool Parameters::load(const std::string &filename, const std::string &overrides, const std::string &fixes, const std::string &releases, bool ignore_unknown){
    std::ifstream input(filename.c_str());
    if(!input){
        std::cerr << "ARRRRRRRRR: Thy parrrrrameter file could not be opened, abandoning ship" << std::endl;
        return false;
    }else{
        std::cout << "Aye, loading thy parrrameters from " << filename << std::endl;
    }
    try{
        load_stream(input, ignore_unknown);
        input.close();
        if(!overrides.empty()) overrideParameters(overrides);
        if(!fixes.empty()) fixParameters(fixes);
        if(!releases.empty()) releaseParameters(releases);
    }catch(std::invalid_argument &e){
        std::cerr << e.what() << std::endl;
        return false;
    }
    return true;
}

void Parameters::save(const std::string& filename) const {
    ofstream output(filename.c_str());
    if(!output) {
        std::cerr << "Could not open parameter file '" << filename << "'" << std::endl;
        return;
    }
    save_stream(output);
}


void Parameters::load_stream(istream &in, bool ignore_unknown){
    int lineNr=0;
    //Read in the whole file and store enough information that we can reproduce it completely (modulo white space change)
    while(!in.eof()){
        //Keep track of the line number for error messages
        ++lineNr;
        //Get the line, if end of file we are done
        string line;
        getline(in,line);
        if(line.empty() && in.eof()) return;
        //Check if there is a comment in the line. If so split comment from rest of the line
        std::string comment;
        size_t comment_pos = line.find('#');
        if(comment_pos!=string::npos){
            //Ignore first line because it is normally just the header which we put in ourselfes
            if(line[0]=='#' && lineNr==1) continue;
            comment = line.substr(comment_pos);
            line = line.substr(0,comment_pos);
            trim(comment);
        }
        //So, remove remaining whitespace and check if there is a parameter left in the line. If so remember which parameter
        trim(line);
        int index(-1);
        if(!line.empty()) {
            string name;
            stringstream buffer(line);
            buffer >> name;
            try{
                index = ParameterList::getIndex(name);
                m_parameters[index].load(buffer);
            }catch(bad_lexical_cast &e){
                cerr << "ERROR reading parameter '" << name << "' (line " << lineNr << "): " << e.what() << endl;
                throw e;
            }catch(invalid_argument &e){
                cerr << "ERROR: Unknown parameter '" << name << "' (line " << lineNr << ")\n";
                if(!ignore_unknown) throw e;
            }
        }
        //Remember which parameter (-1 for none) and the comment for this line
        m_originalLines.push_back(std::make_pair(index,comment));
    }
}

void Parameters::save_stream(ostream &out, bool istty, bool onlyfree) const {
    //If istty is set we want nice color output of the non fixed parameters only
    if(istty){
        out << ANSI_BLUE << boost::format("# Name %|32t| %=13s %=13s %6s %6s %9s %9s\n")
            % "value" % "error" % "min" % "max" % "sig" % "change" << ANSI_END;
    }else{
        out << boost::format("# Name %|32t| %=13s %=13s %6s %6s %5s\n")
            % "value" % "error" % "min" % "max" % "fixed";
    }
    //Bookkeeping to see if we have saved all parameters
    std::vector<bool> saved(m_parameters.size(),false);
    size_t nsaved(0);
    //Go over all the info saved when loading and print the corresponding
    //parameter and comment for each line
    for(const line_info &l: m_originalLines){
        //Write parameter if present
        if(l.first>=0){
            const Parameter &p = m_parameters[l.first];
            //Remember we saved this parameter
            if(!saved[l.first]) ++nsaved;
            saved[l.first] = true;
            //If we only show fre parameters: skip fixed ones
            if(onlyfree && p.isFixed()) continue;
            //If the parameter is fixed we don't show it or its comment in tty
            //mode, go to the next line
            if(!p.save(out, istty)) continue;
        }
        //Write comment if present
        if(!l.second.empty()){
            if(l.first>=0) out << " ";
            if(istty) out << ANSI_BLUE;
            out << l.second;
            if(istty) out << ANSI_END;
        }
        //Write newline when appropriate
        if(!(istty || onlyfree) || l.first>=0 || !l.second.empty()) out << std::endl;
    }

    //Print all parameters which are not saved yet
    if(nsaved < m_parameters.size()){
        if(istty) out << ANSI_BLUE;
        out << std::endl << "# Parameters not present in initial file" << std::endl;
        if(istty) out << ANSI_END;
        for(size_t i=0; i<m_parameters.size(); ++i){
            if(!saved[i]) out << m_parameters[i] << std::endl;
        }
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
    for(Parameter& p: m_parameters){
        if(boost::regex_match(p.name(),fixed)){
            p.setDynamicFix(+1);
        }
    }
}

void Parameters::releaseParameters(const std::string& fixParameters){
    boost::regex fixed(fixParameters);
    for(Parameter& p: m_parameters){
        if(boost::regex_match(p.name(),fixed)){
            p.setDynamicFix(-1);
        }
    }
}

void Parameters::overrideParameters(const std::string &overrides){
    boost::char_separator<char> sep(", ");
    boost::tokenizer<boost::char_separator<char> > tokens(overrides,sep);
    for(const string & tok: tokens) {
        size_t pos = tok.find(':');
        if(pos==string::npos)
            throw std::invalid_argument("Override '" + tok + "' not well formed, expected <name>:<value>");
        Parameter &p = (*this)[tok.substr(0,pos)];
        try{
            const double value = boost::lexical_cast<double>(tok.substr(pos+1));
            p.override(value);
        }catch(bad_lexical_cast &e){
            throw std::invalid_argument("Override '" + tok + "' not well formed, value cannot be parsed");
        }
    }
}

ROOT::Minuit2::MnUserParameters Parameters::getMnParams() const {
    ROOT::Minuit2::MnUserParameters mnParams;

    for(const Parameter& p: m_parameters){
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

void Parameters::update(const ROOT::Minuit2::MnUserParameters &mnParams){
    for(const ROOT::Minuit2::MinuitParameter &mnp: mnParams.Parameters()){
        int index = ParameterList::getIndex(mnp.GetName());
        Parameter &p = m_parameters[index];
        if(!p.isFixed()) {
            p.error(mnp.Error());
            p.value(mnp.Value());
        }
    }
}

void Parameters::update(const std::vector<double> &params){
    assert(params.size() == m_parameters.size());
    for(size_t i=0; i<params.size(); ++i){
        m_parameters[i].value(params[i]);
    }
}

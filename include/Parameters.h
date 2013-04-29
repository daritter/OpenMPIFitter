#ifndef MPIFitter_Parameters_H
#define MPIFitter_Parameters_H

#include <vector>
#include <string>

#include "Parameter.h"
#include "ParameterList.h"

#include "Minuit2/MnUserParameters.h"

class Parameters {
    public:
        typedef std::vector<Parameter>::iterator iterator;
        typedef std::pair<int, std::string> line_info;

        Parameters();

        bool load(const std::string &filename, const std::string &overrides="", const std::string &fixes="", const std::string &releases="");
        //Load the parameters from a stream, usually used with the >> operator
        void load_stream(std::istream &in);
        void save(const std::string &filename) const;
        //Save the parameters to a stream, is istty is true, use color and show only non-fixed parameters
        void save_stream(std::ostream &out, bool istty=false, bool onlyfree=false) const;
        //Short hand for printing with color
        void print(bool istty=true, bool onlyfree=false) const { save_stream(std::cout, istty, onlyfree); }
        //Dynamically fix the parameters which match the given regular expression
        void fixParameters(const std::string& fixParameters = "");
        //Dynamically release the parameters which match the given regular expression
        void releaseParameters(const std::string& releaseParameters = "");
        //Dynamically override parameters
        void overrideParameters(const std::string &overrides);
        //Return Minut2 user parameter object
        ROOT::Minuit2::MnUserParameters getMnParams() const;
        //Update the parameters from the Minuit2 user parameter object
        void update(const ROOT::Minuit2::MnUserParameters mnParams);
        //Return vector containing the values for all parameters
        std::vector<double> getValues() const;
        //Return the number of parameters
        size_t size() const { return m_parameters.size(); }
        //Access the parameter by id
        Parameter& operator[](int id){ return m_parameters[id]; }
        //Access the parameter by name
        Parameter& operator[](const std::string &name){ return m_parameters[ParameterList::getIndex(name)]; }
        //Iterator to the first parameter
        iterator begin() { return m_parameters.begin(); }
        //Iterator behind the last parameter
        iterator end() { return m_parameters.end(); }
    protected:
        //List of all parameters
        std::vector<Parameter> m_parameters;
        //Information needed to save parameters in same order as they were loaded and restore any comments in the file
        std::vector<line_info> m_originalLines;
};

inline std::istream& operator>>(std::istream &in,  Parameters &p){ p.load_stream(in);  return in; }
inline std::ostream& operator<<(std::ostream &out, const Parameters &p){ p.save_stream(out); return out;}

#endif //MPIFitter_Parameter_H

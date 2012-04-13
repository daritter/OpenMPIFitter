#ifndef DsDsKsFitter_Parameters_H
#define DsDsKsFitter_Parameters_H

#include <string>
#include <boost/array.hpp>

#include "Parameter.h"
#include "ParameterList.h"

//#include "Minuit2/MnUserParameters.h"

class Parameters {
    public:

        Parameters();

        void load(std::istream &in);
        void save(std::ostream &out) const;

        //MnUserParameters getMnParams() const;
        //update(const MnUserParameters mnParams);

        Parameter& operator[](int id){ return m_parameters[id]; }
        Parameter& operator[](const std::string &name){ return m_parameters[ParameterList::getIndex(name)]; }
    protected:
        std::vector<Parameter> m_parameters;
};

inline std::istream& operator>>(std::istream &in,  Parameters &p){ p.load(in);  return in; }
inline std::ostream& operator<<(std::ostream &out, const Parameters &p){ p.save(out); return out;}

#endif //DsDsKsFitter_Parameter_H

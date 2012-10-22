#ifndef MPIFitter_Parameter_H
#define MPIFitter_Parameter_H

#include <limits>
#include <string>
#include <iostream>

class Parameter {
    public:
        Parameter(const std::string &name="", double value=0, double error=0,
                double min=-std::numeric_limits<double>::infinity(),
                double max=std::numeric_limits<double>::infinity(), bool fixed=false):
            m_name(name),m_value(value),m_error(error),m_min(min),m_max(max),m_fixed(fixed), m_dynfix(false){};

        const  std::string &name() const { return m_name; }
        void   value(double value) { m_value = value; }
        double value() const { return m_value; }
        double error() const { return m_value; }
        void   error(double error)  { m_error = error; }
        double min() const { return m_min; }
        double max() const { return m_max; }
        bool   isFixed() const { return m_fixed || m_min>=m_max || m_dynfix; }
        bool   hasMin() const { return m_min != -std::numeric_limits<double>::infinity(); }
        bool   hasMax() const { return m_max !=  std::numeric_limits<double>::infinity(); }

        void load(std::istream& in);
        void save(std::ostream& out) const;
        void setDynamicFix(bool fix){ m_dynfix = fix; }
    protected:
        std::string m_name;
        double m_value;
        double m_error;
        double m_min;
        double m_max;
        bool   m_fixed;
        bool   m_dynfix;
};
inline std::istream& operator>>(std::istream &in,  Parameter &p){ p.load(in);  return in; }
inline std::ostream& operator<<(std::ostream &out, const Parameter &p){ p.save(out); return out;}

#endif //MPIFitter_Parameter_H

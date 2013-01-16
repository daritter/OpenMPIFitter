#ifndef MPIFitter_Parameter_H
#define MPIFitter_Parameter_H

#include <limits>
#include <string>
#include <iostream>
#include <cmath>

#define ANSI_RED    "\033[91m"
#define ANSI_GREEN  "\033[92m"
#define ANSI_YELLOW "\033[93m"
#define ANSI_BLUE   "\033[94m"
#define ANSI_PURPLE "\033[95m"
#define ANSI_CYAN   "\033[96m"
#define ANSI_WHITE  "\033[97m"
#define ANSI_END    "\033[0m"

class Parameter {
    public:
        //Create a new parameter, also works as default constructor
        Parameter(const std::string &name="", double value=0, double error=0,
                double min=-std::numeric_limits<double>::infinity(),
                double max=std::numeric_limits<double>::infinity(), bool fixed=false):
            m_name(name),m_value(value),m_error(error),m_min(min),m_max(max),m_fixed(fixed), m_dynfix(0), m_oldvalue(0){};

        //Return the name of the parameter
        const  std::string &name() const { return m_name; }
        //Set the value
        void   value(double value) { m_value = value;}
        //Override the value, also setting it start value
        void   override(double value) { m_value = value; m_oldvalue = value;}
        //Get the value
        double value() const { return m_value; }
        //Set the error
        void   error(double error)  { m_error = error; }
        //Get the error
        double error() const { return m_error; }
        //Get the lower limit
        double min() const { return m_min; }
        //Get the upper limit
        double max() const { return m_max; }
        //Check if the parameter is fixed
        bool   isFixed() const { return m_dynfix>0 || (m_fixed && m_dynfix>=0) || m_min>=m_max; }
        //Check if ther is a lower limit
        bool   hasMin() const { return m_min != -std::numeric_limits<double>::infinity(); }
        //Check if ther is a upper limit
        bool   hasMax() const { return m_max !=  std::numeric_limits<double>::infinity(); }
        //Load the parameter from the given stream
        void load(std::istream& in);
        //Save the parameter to the given stream, if istty=true, only save non fixed and use color
        bool save(std::ostream& out, bool istty=false) const;
        //Set a dynamic fix: Parameter will be fixed in fit but this will not be saved to file
        void setDynamicFix(int fix){ m_dynfix = fix; }

    protected:
        //name of the parameter
        std::string m_name;
        //value of the parameter
        double m_value;
        //error of the parameter
        double m_error;
        //lower limit of the parameter
        double m_min;
        //upper limit of the parameter
        double m_max;
        //indicates wether the parameter is permanently fixed
        bool   m_fixed;
        //indicates wether the parameter is temporarily fixed (>0) or released (<0)
        int   m_dynfix;
        //indicates wether the parameter was changed
        double m_oldvalue;
};
inline std::istream& operator>>(std::istream &in,  Parameter &p){ p.load(in);  return in; }
inline std::ostream& operator<<(std::ostream &out, const Parameter &p){ p.save(out); return out;}

#endif //MPIFitter_Parameter_H

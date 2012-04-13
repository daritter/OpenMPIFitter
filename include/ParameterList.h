#ifndef DsDsKsFitter_ParameterList_h
#define DsDsKsFitter_ParameterList_h

#include <map>
#include <vector>
#include <string>
#include <stdexcept>

struct ParameterList {
    public:
        typedef std::map<std::string, int> map_type;
        typedef std::vector<std::string> list_type;

        /** Add a parameter to both lists */
        static int addParameter(const std::string &name) {
            int index = m_parameterList.size();
            if(!m_parameterMap.insert(make_pair(name,index)).second){
                throw std::logic_error("Duplicate Parameter: " + name);
            }
            m_parameterList.push_back(name);
            return index;
        }

        /** Get parameter index for a given name */
        static int getIndex(const std::string &name) {
            map_type::const_iterator it = m_parameterMap.find(name);
            if(it == m_parameterMap.end()) throw std::invalid_argument("Unknown parameter: " + name);
            return it->second;
        }

        /** Get parameter name for a given index */
        static std::string getName(int index) {
            return m_parameterList.at(index);
        }

        static const map_type & getParameterMap() { return m_parameterMap; }
        static const list_type & getParameterList() { return m_parameterList; }

    private:
        /** Static class, hide constructor */
        ParameterList();
        /** Static class, hide destructor */
        ~ParameterList();
        /** Static class, hide copy constructor */
        ParameterList(const ParameterList&);

        /** Map providing mapping between parameter name and parameter index */
        static map_type  m_parameterMap;
        /** List providing mapping between parameter index and parameter name */
        static list_type m_parameterList;
};

/** Macro to define a parameter easily */
#define PARAM(x) const size_t x = ParameterList::addParameter(#x)

#endif //DsDsKsFitter_ParameterList_H


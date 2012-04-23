#ifndef MPIFitter_ParameterList_h
#define MPIFitter_ParameterList_h

#include <map>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>

/** This class keeps a name<->index matching for all known parameters */
struct ParameterList {
    public:
        typedef std::map<std::string, int> map_type;
        typedef std::vector<std::string> list_type;

        /** Add a parameter to both lists */
        static int addParameter(const std::string &name) {
            ParameterList &pl = ParameterList::getInstance();
            int index = pl.m_parameterList.size();
            if(!pl.m_parameterMap.insert(make_pair(name,index)).second){
                throw std::logic_error("Duplicate Parameter: " + name);
            }
            pl.m_parameterList.push_back(name);
            return index;
        }

        /** Get parameter index for a given name */
        static int getIndex(const std::string &name) {
            ParameterList &pl = ParameterList::getInstance();
            map_type::const_iterator it = pl.m_parameterMap.find(name);
            if(it == pl.m_parameterMap.end()) throw std::invalid_argument("Unknown parameter: " + name);
            return it->second;
        }

        /** Get parameter name for a given index */
        static std::string getName(int index) {
            ParameterList &pl = ParameterList::getInstance();
            return pl.m_parameterList.at(index);
        }

        static const map_type& getParameterMap() { return  ParameterList::getInstance().m_parameterMap; }
        static const list_type& getParameterList() { return  ParameterList::getInstance().m_parameterList; }

    private:
        static ParameterList& getInstance();

        /** Static class, hide constructor */
        ParameterList() {}
        /** Static class, hide destructor */
        ~ParameterList() {}
        /** Static class, hide copy constructor */
        ParameterList(const ParameterList&);

        /** Map providing mapping between parameter name and parameter index */
        map_type  m_parameterMap;
        /** List providing mapping between parameter index and parameter name */
        list_type m_parameterList;

        friend class std::auto_ptr<ParameterList>;
};

/** Macro to define a parameter easily */
#define PARAM(x) static const size_t x = ParameterList::addParameter(#x)

#endif //MPIFitter_ParameterList_H


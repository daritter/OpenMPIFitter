#include "ParameterList.h"
#include <memory>

ParameterList& ParameterList::getInstance() {
    static std::auto_ptr<ParameterList> instance(new ParameterList());
    return *instance;
}

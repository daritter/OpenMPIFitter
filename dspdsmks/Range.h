#ifndef MPIFitter_Range_h
#define MPIFitter_Range_h

#include <iostream>
#include <string>

struct Range {
    Range(const std::string &name, double vmin, double vmax):name(name),vmin(vmin),vmax(vmax) {}
    bool operator()(double value){ return value>=vmin &&  value<=vmax; }

    std::string name;
    float vmin;
    float vmax;

};

inline std::ostream& operator<<(std::ostream& out, const Range &r) {
    out << r.vmin << " <= " << r.name << " <= " << r.vmax;
    return out;
}

#endif

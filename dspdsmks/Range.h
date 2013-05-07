#ifndef MPIFitter_Range_h
#define MPIFitter_Range_h

#include <iostream>
#include <string>

struct Range {
    Range(double vmin, double vmax):vmin(vmin),vmax(vmax) {}
    bool operator()(double value){ return value>=vmin &&  value<=vmax; }

    float vmin;
    float vmax;

    void print(std::ostream& out, const std::string &name) {
        out << vmin << " <= " << name << " <= " << vmax;
    }
};

#endif

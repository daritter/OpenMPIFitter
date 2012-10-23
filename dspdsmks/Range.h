#ifndef MPIFitter_Range_h
#define MPIFitter_Range_h

struct Range {
    Range(double vmin, double vmax):vmin(vmin),vmax(vmax) {}
    bool operator()(double value){ return value>=vmin &&  value<=vmax; }

    float vmin;
    float vmax;
};

#endif

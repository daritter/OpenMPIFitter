#ifndef PROGRESS_H
#define PROGRESS_H

#include <cmath>

/** Progressbar class
 *
 * Usage:
 * ProgressBar pbar(100);
 * for(int i=0; i<100; ++i){
 *     ... do something ...
 *
 *     std::cout << pbar++;
 * }
 */
class ProgressBar {
    protected:
        int maxval;
        int lastpos;
        int curpos;
        int lastVal;
    public:
        ProgressBar(int max_value):maxval(max_value),lastpos(-1),curpos(0), lastVal(0){};
        ProgressBar &  operator() (int value){
            curpos = (int)std::min(round(value*60.0/maxval),60.0);
            return *this;
        }
        ProgressBar & operator++() {
            return (*this)(++lastVal);
        }
        void show(std::ostream &out){
            if(curpos>lastpos){
                if(lastpos < 0) out<<"[";
                for(int i=lastpos+1;i<=curpos;i++){
                    if(i%15==0){
                        out << (i/60.*100) << '%';
                    }else{
                        out << (i%5==0?'+':'-');
                    }
                }
                if(curpos>=60) out<<"]"<< std::endl;
                out << std::flush;
            }
            lastpos = curpos;
        };
};

std::ostream & operator<<(std::ostream &out, ProgressBar &pb){
    pb.show(out);
    return out;
}

#endif

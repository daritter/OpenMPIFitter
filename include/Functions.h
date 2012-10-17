#ifndef OpenMPIFitter_Functions_h
#define OpenMPIFitter_Functions_h

#include <cmath>
#include <limits>

class Gauss {
    public:
        static const double sqrt_2pi_2 = 1.2533141373155002512078826424055226265034;

        Gauss(double lower, double upper):
            lower(lower),upper(upper),mean(std::numeric_limits<double>::quiet_NaN()),sigma(-1),sigma2(-1),norm(0) {}

        /** Setup the parameters */
        void set(double mean, double sigma, double sigma2 = -1){
            //std::cout << "params: " << mean << ", " << sigma << ", " << sigma2 << std::endl;
            //We only have to recalculate if the parameters changed since the last call
            if(mean == this->mean && sigma == this->sigma && sigma2 == this->sigma2) return;

            // If one of the sigmas is 0 we have a delta function
            if(sigma == 0 || sigma2 == 0) {
                norm = 0;
            } else {
                //Otherwise we calculate the normalization factor in the given bounds
                const double x1 = (upper - mean)*M_SQRT1_2;
                const double x2 = (lower - mean)*M_SQRT1_2;
                if(sigma2 == -1 || upper < mean) {
                    norm = sigma * ( erf(x1/sigma) - erf(x2/sigma) );
                }else if( lower > mean) {
                    norm = sigma2 * ( erf(x1/sigma2) - erf(x2/sigma2) );
                }else{
                    norm = sigma2 * erf(x1/sigma2) - sigma * erf(x2/sigma);
                }
                norm *= sqrt_2pi_2;
            }

            //std::cout << "norm: " << norm << ", " << lower << ", " << upper << std::endl;
            //And save the parameters
            this->mean = mean;
            this->sigma = sigma;
            this->sigma2 = sigma2;
        }

        /** Calculate the fcn for a given value */
        double operator()(double x) const {
            //std::cout << "gauss: " << sigma << ", " << norm << std::endl;
            if(sigma==0.0 || sigma2==0.0) return (x==mean)?std::numeric_limits<double>::max():0;
            if(!finite(sigma) || norm==0) return 0.0;
            const double s = (sigma2!=-1 && x>mean)?sigma2:sigma;
            const double e = (x-mean)/s;
            return exp(-0.5*e*e)/norm;
        }

        double getNorm() const { return norm; }

    protected:
        double lower;
        double upper;
        double mean;
        double sigma;
        double sigma2;
        double norm;
};

class Argus {
    public:
        Argus(double lower, double upper):
            lower(lower),upper(upper),a(std::numeric_limits<double>::quiet_NaN()), benergy(0), norm(0) {}

        void set(double benergy, double a){
            if(a == this->a && benergy == this->benergy) return;

            const double e2 = benergy * benergy;
            const double ll2 = 1.0 - (lower*lower) / e2;
            const double ul2 = std::max(0.0, 1.0 - (upper*upper) / e2);
            if(ll2 < 0.0) {
                norm = 0.0;
                return;
            }

            norm = 0.5*e2/a*( (sqrt(ll2)*exp(a*ll2) - sqrt(ul2)*exp(a*ul2))
                 + 0.5*sqrt(M_PI/fabs(a)) * (erfc(sqrt(fabs(a)*ll2)) - erfc(sqrt(fabs(a)*ul2))));

            this->a = a;
            this->benergy = benergy;
        }

        /** Calculate the fcn for a given value */
        double operator()(double x) const {
            if ( x > benergy ) return 0.0;
            const double dz2 = 1.0 - ( x/benergy ) * ( x/benergy );
            return (x * sqrt(dz2) * exp(a*dz2))/norm;
        }

        double getNorm() const { return norm; }

    protected:
        double lower;
        double upper;
        double a;
        double benergy;
        double norm;
};

class Chebychev1 {
    public:
        Chebychev1(double lower, double upper):
            lower(lower),upper(upper),c(std::numeric_limits<double>::quiet_NaN()),norm(0) {}

        void set(double c){
            if(c == this->c) return;

            const double ll2 = lower*lower;
            const double ul2 = upper*upper;
            norm = upper - lower + (c*0.5*ul2) - (c*0.5*ll2);
            this->c = c;
        }

        /** Calculate the fcn for a given value */
        double operator()(double x) const {
            return (1.0 + (c*x))/norm;
        }

        double getNorm() const { return norm; }

    protected:
        double lower;
        double upper;
        double c;
        double norm;
};

/** Compound FCN for adding to generic 1D FCN functions*/
template<class FCN1, class FCN2> class Add1DFcn {
    public:
        Add1DFcn(double lower, double upper): fcn1(lower,upper), fcn2(lower,upper), ratio(1) {}

        /** Set the parameters. The parameters of the components have to be set separately by calling
         * fcn1.set and fcn2.set
         */
        void set(double ratio) {
            this->ratio = ratio;
        }

        /** Calculate the fcn for a given value */
        double operator()(double x) const {
            return ratio*fcn1(x) + (1.0-ratio)*fcn2(x);
        }

        FCN1 fcn1;
        FCN2 fcn2;
    protected:
        double ratio;
};


template<class FCN1, class FCN2> class Add2DFcn {
    public:
        Add2DFcn(double lowerX, double upperX, double lowerY, double upperY):
            fcn1(lowerX,upperX,lowerY,upperY), fcn2(lowerX,upperX,lowerY,upperY), ratio(0.5) {}

        /** Set the parameters. The parameters of the components have to be set separately by calling
         * fcn1.set and fcn2.set
         */
        void set(double ratio) {
            this->ratio = ratio;
        }

        /** Calculate the fcn for a given value */
        double operator()(double x, double y) const {
            return ratio*fcn1(x,y) + (1.0-ratio)*fcn2(x,y);
        }

        FCN1 fcn1;
        FCN2 fcn2;
    protected:
        double ratio;
};

template<class FCNX, class FCNY> class CompoundFcn2D {
    public:
        CompoundFcn2D(double lowerX, double upperX, double lowerY, double upperY): fcnx(lowerX,upperX), fcny(lowerY,upperY) {}

        /** Calculate the fcn for a given value */
        double operator()(double x, double y) const {
            return fcnx(x) * fcny(y);
        }

        FCNX fcnx;
        FCNY fcny;
};


class DoubleGauss: public Add1DFcn<Gauss, Gauss> {
    public:
        DoubleGauss(double lower, double upper): Add1DFcn<Gauss, Gauss>(lower,upper) {}

        void set(double ratio, double mean, double meanshift, double sigma, double sigmascale, double sigma2=-1, double sigma2scale=-1){
            Add1DFcn<Gauss, Gauss>::set(ratio);
            fcn1.set(mean, sigma, sigma2);
            if(sigma2 != -1){
                if(sigma2scale == -1) sigma2scale = sigmascale;
                sigma2 *= sigma2scale;
            }
            fcn2.set(mean+meanshift, sigma*sigmascale, sigma2);
        }
};

#endif //OpenMPIFitter_Functions_h

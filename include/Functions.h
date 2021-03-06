#ifndef OpenMPIFitter_Functions_h
#define OpenMPIFitter_Functions_h

#include <cmath>
#include <limits>
#include <stdexcept>

class Gauss {
    public:
        constexpr static double sqrt_pi_2 = 1.2533141373155002512078826424055226265034;

        Gauss(double lower=std::numeric_limits<double>::quiet_NaN(), double upper=0):
            lower(lower),upper(upper),mean(std::numeric_limits<double>::quiet_NaN()),sigma(-1),sigma2(-1),norm(0) {}

        void set_limits(double lower, double upper){
            if(lower == this->lower && upper == this->upper) return;
            this->lower = lower;
            this->upper = upper;
            mean = std::numeric_limits<double>::quiet_NaN();
        }

        /** Setup the parameters */
        void set(double mean, double sigma, double sigma2 = -1){
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
                norm *= sqrt_pi_2;
            }

            //And save the parameters
            this->mean = mean;
            this->sigma = sigma;
            this->sigma2 = sigma2;
        }

        /** Calculate the fcn for a given value */
        double operator()(double x) const {
            if(x<lower || x>upper) return 0.0;
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

        void set_limits(double lower, double upper){
            if(lower == this->lower && upper == this->upper) return;
            this->lower = lower;
            this->upper = upper;
            a = std::numeric_limits<double>::quiet_NaN();
        }

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
            if(x<lower || x>upper) return 0.0;
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

template<int N> class Chebychev {
    public:
        Chebychev(double lower, double upper): norm(0) {
            if(N<0) throw std::logic_error("Chebychev with N<0 makes no sense");
            if(N>4) throw std::logic_error("Chebychev with N>4 not supported");
            u[0] = std::numeric_limits<double>::quiet_NaN();
            set_limits(lower,upper);
        }

        void set_limits(double lower, double upper){
            if(lower == l[0] && upper == u[0]) return;
            for(int i=0; i<N; ++i) c[i] = std::numeric_limits<double>::quiet_NaN();
            l[0] = lower;
            u[0] = upper;
            for(int i=1; i<N+1; ++i){
                l[i] = l[i-1] * lower;
                u[i] = u[i-1] * upper;
            }
        }

        void set(const double* c){
            bool changed=false;
            for(int i=0; i<N; ++i){
                if(this->c[i] != c[i]) changed=true;
            }
            if(!changed) return;
            std::copy(c,c+N,this->c);

            norm = u[0] - l[0] +
                ((N<1) ? 0 : c[0]*0.5*(u[1] - l[1]) ) +
                ((N<2) ? 0 : c[1]*(((2.0/3.0*u[2]) - u[0]) - ((2.0/3.0*l[2]) - l[0])) ) +
                ((N<3) ? 0 : c[2]*((u[3] - (3.0/2.0*u[1])) - (l[3] - (3.0/2.0*l[1]))) ) +
                ((N<4) ? 0 : c[3]*(((8.0/5.0*u[4]) - (8.0/3.0*u[2]) + u[0]) - ((8.0/5.0*l[4]) - (8.0/3.0*l[2]) + l[0])) );
        }

        /** Calculate the fcn for a given value */
        double operator()(double x) const {
            if(x<l[0] || x>u[0]) return 0;

            double xp[N];
            xp[0]=x;
            for(int i=1; i<N; ++i) xp[i]=xp[i-1]*x;

            double value = 1.0 +
                ((N<1) ? 0 : c[0]*xp[0]) +
                ((N<2) ? 0 : c[1]*((2.0*xp[1])-1.0)) +
                ((N<3) ? 0 : c[2]*((4.0*xp[2]) - (3.0*xp[0]))) +
                ((N<4) ? 0 : c[3]*((8.0*xp[3]) - (8.0*xp[1]) + 1.0));
            return value/norm;
        }

        double getNorm() const { return norm; }

    protected:
        double c[N];
        double l[N+1];
        double u[N+1];
        double norm;
};

//FIXME Check function
class CrystalBall {
    public:
        CrystalBall(double lower, double upper):
            lower(lower),upper(upper), alpha(std::numeric_limits<double>::quiet_NaN()), n(0), mean(0), sigma(0) {}

        void set_limits(double lower, double upper){
            if(lower == this->lower && upper == this->upper) return;
            this->lower = lower;
            this->upper = upper;
            alpha = std::numeric_limits<double>::quiet_NaN();
        }

        void set(double alpha, double n, double mean, double sigma){
            if( this->alpha==alpha && this->n == n && this->mean==mean && this->sigma==sigma)
                return;

            this->alpha = alpha;
            this->n = n;
            this->mean = mean;
            this->sigma = sigma;
            A = std::pow(n/alpha,n) * std::exp(-alpha*alpha/2.);
            B = (n/alpha) - alpha;
            if(sigma==0) {
                norm = 0;
            } else {
                norm = 1./(integral(upper) - integral(lower));
            }
        }

        /** Calculate the fcn for a given value */
        double operator()(double x) const {
            if(x<lower || x>upper) return 0.0;
            if(sigma==0.0) return (x==mean)?std::numeric_limits<double>::max():0;
            if(!finite(sigma) || norm==0) return 0.0;
            const double dev = (x-mean)/sigma;
            if(dev>-alpha) return norm*std::exp(-dev*dev/2.);
            return norm*A*std::pow(B-dev,n);
        }

        double getNorm() const { return norm; }

    protected:
        double integral(double x){
            const double flip = mean-alpha*sigma;
            if(x<=flip) return integral_exp(x);
            return integral_exp(flip) + integral_gauss(x) - integral_gauss(flip);
        }

        double integral_exp(double x){
            return A*sigma*std::pow((mean+B*sigma-x)/sigma,1-n)/(n-1);
        }

        double integral_gauss(double x){
            return Gauss::sqrt_pi_2*(1+std::erf((x-mean)/(M_SQRT2*sigma)))*sigma;
        }

        double lower;
        double upper;
        double alpha;
        double n;
        double mean;
        double sigma;
        double A;
        double B;
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
            if(ratio==1) return fcn1(x);
            if(ratio==0) return fcn2(x);
            return ratio*fcn1(x) + (1.0-ratio)*fcn2(x);
        }

        void set_limits(double lower, double upper){
            fcn1.set_limits(lower, upper);
            fcn2.set_limits(lower, upper);
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
            if(ratio==1) return fcn1(x,y);
            if(ratio==0) return fcn2(x,y);
            return ratio*fcn1(x,y) + (1.0-ratio)*fcn2(x,y);
        }

        void set_limits(double lowerX, double upperX, double lowerY, double upperY){
            fcn1.set_limits(lowerX, upperX, lowerY, upperY);
            fcn2.set_limits(lowerX, upperX, lowerY, upperY);
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

        void set_limits(double lowerX, double upperX, double lowerY, double upperY){
            fcnx.set_limits(lowerX, upperX);
            fcny.set_limits(lowerY, upperY);
        }

        FCNX fcnx;
        FCNY fcny;
};


class DoubleGauss: public Add1DFcn<Gauss, Gauss> {
    public:
        DoubleGauss(double lower, double upper): Add1DFcn<Gauss, Gauss>(lower,upper) {}

        void set(double ratio, double mean, double meanshift, double sigma, double sigmascale, double sigma2=-1, double sigma2scale=-1){
            Add1DFcn<Gauss, Gauss>::set(ratio);
            double s2 = sigma2;
            if(sigma2!=-1) s2 = sigma*sigma2;
            if(ratio!=0) fcn1.set(mean, sigma, s2);
            if(sigma2 != -1){
                if(sigma2scale == -1) sigma2scale = sigmascale;
                s2 *= sigma2scale;
            }
            if(ratio!=1) fcn2.set(mean+meanshift, sigma*sigmascale, s2);
        }
};

template<int N> class MultiGauss {
    public:
        MultiGauss(double lower, double upper) {
            fcns = new Gauss[N];
            set_limits(lower, upper);
            norms[0]=1;
        }

        ~MultiGauss(){
            delete[] fcns;
        }

        void set(const double *par1){
            double mean(0);
            double sigma(1.0);
            for(int i=0; i<N; ++i){
                if(i>0) norms[i] = par1[3*i-1];
                mean += par1[3*i];
                sigma *= par1[3*i+1];
                if(norms[i]==0) continue;
                fcns[i].set(mean,sigma);
            }
        }

        void set_limits(double lower, double upper){
            for(int i=0; i<N; ++i){
                fcns[i].set_limits(lower,upper);
            }
        }

        void set(double mean, double sigma, const double *par1=0){
            fcns[0].set(mean,sigma);
            for(int i=1; i<N; ++i){
                norms[i] = par1[3*i-3];
                mean += par1[3*i-2];
                sigma *= par1[3*i-1];
                if(norms[i]==0) continue;
                fcns[i].set(mean,sigma);
            }
        }

        void set(double mean, double sigma, double norm1, double meanshift1, double sigmascale1, const double *par1=0){
            fcns[0].set(mean,sigma);
            mean += meanshift1;
            sigma *= sigmascale1;
            norms[1] = norm1;
            fcns[1].set(mean,sigma);
            for(int i=2; i<N; ++i){
                norms[i] = par1[3*i-6];
                mean += par1[3*i-5];
                sigma *= par1[3*i-4];
                if(norms[i]==0) continue;
                fcns[i].set(mean,sigma);
            }
        }

        void set(double mean, double sigma, double norm1, double meanshift1, double sigmascale1, double norm2, double meanshift2, double sigmascale2, const double *par1=0){
            fcns[0].set(mean,sigma);
            mean += meanshift1;
            sigma *= sigmascale1;
            norms[1] = norm1;
            fcns[1].set(mean,sigma);

            mean += meanshift2;
            sigma *= sigmascale2;
            norms[2] = norm2;
            fcns[2].set(mean,sigma);

            for(int i=3; i<N; ++i){
                norms[i] = par1[3*i-9];
                mean += par1[3*i-8];
                sigma *= par1[3*i-7];
                if(norms[i]==0) continue;
                fcns[i].set(mean,sigma);
            }
        }

        double operator()(double x) const {
            double result(0);
            double norm(0);
            for(int i=0; i<N; ++i){
                if(norms[i]==0) continue;
                result += norms[i]*fcns[i](x);
                norm += norms[i];
            }
            return result/norm;
        }

    private:
        Gauss* fcns;
        double norms[N];
};

#endif //OpenMPIFitter_Functions_h

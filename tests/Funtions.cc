#include "gtest/gtest.h"
#include "Functions.h"
#include "func.h"
#include <iostream>

double getStep(int i, int N, double min, double max){
    return 1.0*i/(N-1) * (max-min) + min;
}

#define step(var,min,max,N) for(double i##var=0, var##step=1.0*(max-min)/(N-1), var=min; i##var<N; i##var+=1, var=i##var*var##step+min)


TEST(getStepTest,Works){
    EXPECT_DOUBLE_EQ(-5.,getStep(0,101,-5,5));
    EXPECT_DOUBLE_EQ( 5.,getStep(100,101,-5,5));
    EXPECT_DOUBLE_EQ( 0.,getStep(50,101,-5,5));
    int i=0;
    step(x,-5,5,101){
        EXPECT_FLOAT_EQ(x, getStep(i++,101,-5,5));
    }
}

TEST(GaussTest,Mean) {
    Gauss g1(-5,5);
    step(mean,-5,5,1001){
        step(x,-5,5,1001){
            g1.set(mean,1);
            double norm = Belle::norm_gaussian(-5,5,mean,1);
            EXPECT_FLOAT_EQ(g1(x), Belle::gaussian(x, mean, 1)/norm)
                << " µ=" << mean << " s=1 x=" << x << std::endl;
        }
    }
}

TEST(GaussTest,Sigma) {
    Gauss g1(-5,5);
    step(sigma,0,5,1001){
        step(x,-5,5,1001){
            g1.set(0,sigma);
            double norm = Belle::norm_gaussian(-5,5,0,sigma);
            EXPECT_FLOAT_EQ(g1(x), Belle::gaussian(x, 0, sigma)/norm)
                << " µ=0 s=" << sigma << " x=" << x << std::endl;
        }
    }
}

TEST(GaussTest,Limits) {
    for(int ill=0; ill<99; ++ill){
        for(int iul=ill+1; iul<100; ++iul){
            double ul = getStep(iul,100,-10,10);
            double ll = getStep(ill,100,-10,10);
            Gauss g1(ll,ul);
            step(x,ll,ul,1001){
                g1.set(0,3);
                double norm = Belle::norm_gaussian(ll,ul,0,3);
                EXPECT_FLOAT_EQ(g1(x), Belle::gaussian(x,0,3)/norm);
            }
        }
    }
}

TEST(GaussTest,Bifur) {
    for(int ill=0; ill<10; ++ill){
        for(int iul=ill+1; iul<11; ++iul){
            double ul = getStep(iul,10,-10,10);
            double ll = getStep(ill,10,-10,10);
            Gauss g1(ll,ul);
            step(mean,-5,5,11){
                step(sigma1,1e-5,3,10){
                    step(sigma2,1e-5,3,10){
                        step(x,ul,ll,101){
                            g1.set(mean,sigma1,sigma2);
                            double norm = Belle::norm_bigauss(ll,ul,mean,sigma1,sigma2);
                            if(norm == 0) continue;
                            EXPECT_FLOAT_EQ(g1(x), Belle::bigauss(x, mean, sigma1, sigma2)/norm)
                                 << "ll: " << ll << " ul: " << ul << " µ: " << mean << " s1: " << sigma1 << " s2: " << sigma2 << " x: " << x;
                        }
                    }
                }
            }
        }
    }
}

TEST(FuncTest,Argus) {
    Argus a1(5.2,5.3);
    step(benergy,5.15,5.35,101){
        step(a,-5,5,101){
            a1.set(benergy,a);
            double norm = Belle::norm_argus(5.2,5.3,benergy,a);
            EXPECT_EQ(norm!=norm, a1.getNorm()!=a1.getNorm());
            if(norm!=norm) continue;
            EXPECT_FLOAT_EQ(a1.getNorm(), norm);
            if(norm==0) continue;
            step(x,5.2,5.3,1001){
                EXPECT_FLOAT_EQ(a1(x), Belle::argus(x,benergy,a)/norm) << "benergy=" << benergy << ", a=" << a << ", x=" << x << " norm=" << norm;
            }
        }
    }
}

TEST(FuncTest,Cheb1_T){
    Chebychev<1> c1(-5,5);
    step(c,-5,5,1001){
        c1.set(&c);
        double norm = Belle::norm_cheb1(-5,5,c);
        EXPECT_EQ(norm!=norm, c1.getNorm()!=c1.getNorm());
        if(norm!=norm) continue;
        EXPECT_FLOAT_EQ(c1.getNorm(), norm);
        if(norm==0) continue;

        step(x,-5,5,1001){
            EXPECT_FLOAT_EQ(c1(x),Belle::cheb1(x,c)/norm) << "c=" << c << " x=" << x;
        }
    }
}

TEST(FuncTest,Cheb2_T){
    Chebychev<2> cheb(-5,5);
    std::vector<double> cn;
    cn.resize(2);
    step(c1,-5,5,201){
        cn[0] = c1;
        step(c2,-5,5,201){
            cn[1] = c2;
            cheb.set(&cn[0]);
            double norm = Belle::norm_cheb2(-5,5,cn);
            EXPECT_EQ(norm!=norm, cheb.getNorm()!=cheb.getNorm());
            if(norm!=norm) continue;
            EXPECT_FLOAT_EQ(cheb.getNorm(), norm);
            if(norm==0) continue;

            step(x,-5,5,501){
                EXPECT_FLOAT_EQ(cheb(x),Belle::cheb2(x,cn)/norm) << "c=" << c1 << "," << c2 << " x=" << x;
            }
        }
    }
}

TEST(FuncTest,Cheb2_T_Bench){
    Chebychev<2> cheb(-5,5);
    double cn[2];
    step(c1,-5,5,201){
        cn[0] = c1;
        step(c2,-5,5,201){
            cn[1] = c2;
            cheb.set(cn);
            double foo(0);
            step(x,-5,5,5001){
                foo += cheb(x);
            }
            EXPECT_TRUE(fabs(1-(foo/500.))<0.002) << (foo/500.);
        }
    }
}


TEST(FuncTest,Cheb4_T){
    Chebychev<4> cheb(-5,5);
    std::vector<double> cn;
    cn.resize(4);
    step(c1,-5,5,11){
        cn[0] = c1;
        step(c2,-5,5,11){
            cn[1] = c2;
            step(c3,-5,5,11){
                cn[2] = c3;
                step(c4,-5,5,101){
                    cn[3] = c4;
                    cheb.set(&cn[0]);
                    double norm = Belle::norm_cheb4(-5,5,cn);
                    EXPECT_EQ(norm!=norm, cheb.getNorm()!=cheb.getNorm());
                    if(norm!=norm) continue;
                    EXPECT_FLOAT_EQ(cheb.getNorm(), norm);
                    if(norm==0) continue;

                    step(x,-5,5,101){
                        EXPECT_FLOAT_EQ(cheb(x),Belle::cheb4(x,cn)/norm) << "c=" << c1 << "," << c2 << " x=" << x;
                    }
                }
            }
        }
    }
}

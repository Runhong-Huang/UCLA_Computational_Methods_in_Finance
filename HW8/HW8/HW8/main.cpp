//
//  main.cpp
//  HW8
//
//  Created by HUANG runhong on 3/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Calculator {
public:
    
    Calculator(){
    };
    
    ~Calculator(){};
    
    template <class T>
    double average(std::vector<T> arr, int length){
        float sum =0;
        for (int i = 0 ; i  < length ; i++)
            sum = sum + arr[i];
        //    cout << sum << endl;
        return sum/length;
    };
    
    
    template <class T>
    double variance(std::vector<T> arr, int length) {
        //    calculate the variance of array with length of "length"
        double average = this->average(arr, length);
        double var = 0 ;
        for (int i = 0 ; i  < length ; i++)
            var = var + (arr[i] - average) * (arr[i] - average);
        return var/(length-1);
    };
    
    //    function used to calculate correlation
    template <class T,class U>
    double corr(std::vector<T> arr1,std::vector<U> arr2, int length) {
        //    calculate the variance of array with length of "length"
        double average1 = this->average(arr1, length);
        double var1 = this-> variance(arr1, length);
        double average2 = this->average(arr2, length);
        double var2 = this-> variance(arr2, length);
        double cumcov=0,cov = 0 ;
        
        for (int i = 0 ; i  < length ; i++){
            cumcov = cumcov + (arr1[i] - average1) * (arr2[i] - average2);
        }
        cov = cumcov/(length-1);
        double corr = cov/(sqrt(var1 * var2));
        return corr;
    }
};

class Normal{
private:
    int length = 200;
    double seed = 522313;
    double Amiu = 0, Asigma = 1, Bmiu = 0, Bsigma = 1, rou = 0;
    
    
public:
    Normal(){};
    Normal(int num):length(num) {
    }
    ~Normal(){};

    // define two vectors with varaible length.
    std::vector<float> a ;
    std::vector<float> b ;
    
    void setseed(long s){
        this->seed = s;
    };
    
    //    set the mean and varaince of A
    void setA(float miu, float sigma){
        this->Amiu = miu;
        this->Asigma = sigma;
    };
    
    //    set the mean and varaince of B
    void setB(float miu, float sigma){
        this->Bmiu = miu;
        this->Bsigma = sigma;
    }
    
    void setCorr(float correlation){
        this->rou = correlation;
    }
    
    void generator(){
        a.clear();
        b.clear();
        srand(seed);
        for (int i =0 ; i < length;){
            float v1 = rand()*1.0/(RAND_MAX*1.0) * 2 -1  ;
            float v2 = rand()*1.0/(RAND_MAX*1.0) * 2 -1 ;
            float w = v1* v1 + v2* v2;
            if (w <= 1 ) {
                i = i +1;
                double z1 = v1 * sqrt(-2 * log (w) /w);
                double z2 = v2 * sqrt(-2 * log (w) /w);
                a.push_back(Amiu + Asigma * z1);
                b.push_back(Bmiu + Bsigma * (rou* z1+ sqrt(1.0-rou*rou) * z2));
            }
        }
    };
};


class Vasicek {
    
protected:
    
    float r0, sigma, keppa, r_, FV, T;
    
public:
    
    Vasicek(float R0, float Sigma, float Keppa, float R_, float FV, float T): r0(R0), sigma(Sigma), keppa(Keppa), r_(R_), FV(FV), T(T) {
//        cout << R0 << endl;
    };
    ~Vasicek(){};
    
    void set_t(float time) {
        this-> T = time;
    }
    
    virtual float R(float r_start, float T, float t);
    float CIR(float r_start, float T, float t);
    virtual float R_MC(float r_start ,float T , float t );
    float semicoupon_bond(float coupon, float r_start, float T , float t );
    virtual float option_price1(float k , float mat_t);
    virtual float option_price2(float k , float mat_t);
};

class G2 {
    
    float phi,r0, rho, a,  b, sigma, ita;
    float T ;
    float k ;
    float mat_t;
    
public:
    
    G2(float Phi, float R0, float Rho, float A, float B, float Sigma, float Ita): phi(Phi),r0(R0), rho(Rho), a(A), b(B), sigma(Sigma), ita(Ita) {};
    ~G2(){};
    
public:
    
    void set_T (float time){
        this-> T = time;
    }
    void set_K (float k){
        this-> k = k;
    }
    void set_mat (float time){
        this-> mat_t = time;
    }
    
    float R_MC(float x , float y , float r_start ,float time_T, float time_t );
    float option_price();
    
};

float G2::option_price (){
    
    Normal norm;
    
    int n = 365 * mat_t ;
    int npath = 1000;
    float rsum = 1 ;
    float delta_t = mat_t / n ;
    float r[n];
    
    float x[n];
    float y[n];
    x[0] = 0 ;
    y[0] = 0 ;
    r[0] = r0;
    
    std::vector<float> r_values;
    std::vector<float> option;
    std::vector<float> test;
    
    for (int j = 0 ; j < npath ; j++){
        norm.setseed(j);
        norm.setCorr(rho);
        norm.generator();
        rsum = 1 ;
        for ( int i = 0 ; i < n; i ++){
            x[i+1] = x[i] -a * x[i] * delta_t + sigma * sqrt(delta_t) * norm.a[i];
            y[i+1] = y[i] -b * y[i] * delta_t + ita * sqrt(delta_t) * norm.b[i];
            r[i+1] = x[i+1] + y[i+1] + phi;
            rsum = rsum * exp (-r[i] * delta_t);
        }

        r_values.push_back(rsum);
//        cout << r[n] <<" " << rsum << endl;
//
        float bond = max((float) k - R_MC(x[n], y[n], r[n], 0.5 ,0) * 1000 ,(float) 0.0 );
    
//        cout   << R_MC(x[n], y[n], r[n], 0.5 ,0)  << endl;
        
//        cout << bond << " " << rsum << endl;
        
        option.push_back( bond * rsum );
    }

    
    Calculator cal;
    float option_price =  cal.average(option, npath) ;

    
    return option_price;

}

float G2::R_MC (float x_end, float y_end, float r_start ,float T = 0, float t = 0  ) {
    
    int n = 365 * (T-t) ;
    int npath = 2000;
    float rsum = 1 ;
    float delta_t = (T -t) / n ;
    Normal norm (n);

    
    float x[n];
    float y[n];
    float r[n];
    x[0] = x_end ;
    y[0] = y_end ;
    r[0] = r_start ;
    
    std::vector<float> r_values;
    
    for (int j = 0 ; j < npath ; j++){

        norm.setseed(j);
        norm.setCorr(rho);
        norm.generator();
        rsum = 1 ;
        for ( int i = 0 ; i < n; i ++){
            x[i+1] = x[i] -a * x[i] * delta_t + sigma * sqrt(delta_t) * norm.a[i];
            y[i+1] = y[i] -b * y[i] * delta_t + ita * sqrt(delta_t) * norm.b[i];
            r[i+1] = x[i+1] + y[i+1] + phi;
            rsum = rsum * exp (-r[i] * delta_t);
        }
        r_values.push_back(rsum);
    }
    
    Calculator cal;
    float ave = cal.average(r_values, npath);
    return ave;
    
}


class CIR: public Vasicek {
    
public:
    
    CIR(float R0, float Sigma, float Keppa, float R_, float FV, float T)  :  Vasicek(R0,  Sigma,  Keppa,  R_,  FV,  T) {
    };
    
    ~CIR(){};
    
public:
    
    float R_MC(float r_start ,float , float t );
    float option_price2(float k, float mat_t);
};

float CIR::R_MC (float r_start ,float T = 0, float t = 0) {
    int n = 365 * (T-t) ;
    int npath = 1000;
    float rsum = 1 ;
    float delta_t = (T -t) / n ;
    Normal norm (n);
    float r[n];
    r[0] = r_start ;
    
    std::vector<float> r_values;
    
    for (int j = 0 ; j < npath ; j++){
        norm.setseed(j);
        norm.generator();
        rsum = 1 ;
        for ( int i = 0 ; i < n; i ++){
            r[i+1] = r[i] + keppa * ( r_ - r[i]) * delta_t + sigma * sqrt(delta_t * r[i] ) * norm.b[i];
            cout << r[i+1] << endl;
            if (r[i+1] < 0){
                cout << "test" << endl;
                r[i+1] = 0;
            }
            //            cout << exp (-r[i] * delta_t) ;
            rsum = rsum * exp (-r[i] * delta_t);
        }
        r_values.push_back(rsum);
    }
    Calculator cal;
    float ave = cal.average(r_values, npath);
    return ave;
}

float CIR::option_price2(float k, float mat_t) {
    
    int n = 365 * mat_t ;
    int npath = 1000;
    float rsum = 1 ;
    float delta_t = mat_t / n ;
    Normal norm (n);
    float r[n];
    r[0] = r0;
    
    std::vector<float> r_values;
    std::vector<float> option;
    
    for (int j = 0 ; j < npath ; j++){
        //        cout << j << endl;
        norm.setseed(j);
        norm.generator();
        rsum = 1 ;
        for ( int i = 0 ; i < n; i ++){
            r[i+1] = r[i] + keppa * ( r_ - r[i]) * delta_t + sigma * sqrt(delta_t * r[i]) * norm.a[i];
//            cout << r[i] << endl;
            rsum = rsum *exp (-r[i] * delta_t);
        }
        r_values.push_back(rsum);

        float bond = max((float) R_MC(r[n], T- mat_t ,0) * FV - k ,(float) 0.0 );
//                cout << bond << endl;
        option.push_back( bond * rsum );
    }
    
    Calculator cal;
    float option_price =  cal.average(option, npath) ;
    return option_price;
}



float Vasicek::option_price1(float k, float mat_t) {
    
    int n = 365 * mat_t ;
    int npath = 1000;
    float rsum = 1 ;
    float delta_t = mat_t / n ;
    Normal norm (n);
    float r[n];
    r[0] = r0;
    
    std::vector<float> r_values;
    std::vector<float> option;
 
    for (int j = 0 ; j < npath ; j++){
        norm.setseed(j);
        norm.generator();
        rsum = 1 ;
        for ( int i = 0 ; i < n; i ++){
            r[i+1] = r[i] + keppa * ( r_ - r[i]) * delta_t + sigma * sqrt(delta_t) * norm.a[i];
            rsum = rsum *exp (-r[i] * delta_t);
        }
        r_values.push_back(rsum);
        
        float bond = max((float) R(r[n], T - mat_t ,0) * FV - k ,(float) 0.0 );
//        cout << bond << endl;
        option.push_back( bond * rsum );
    }
    

    Calculator cal;
    float option_price =  cal.average(option, npath) ;
    
    return option_price;
}

float Vasicek::option_price2(float k, float mat_t) {
    
    int n = 365 * mat_t ;
    int npath = 1000;
    float rsum = 1 ;
    float delta_t = mat_t / n ;
    Normal norm (n);
    float r[n];
//    float rterm[n];
    r[0] = r0;
    
    std::vector<float> r_values;
    std::vector<float> option;

    
    for (int j = 0 ; j < npath ; j++){
//        cout << j << endl;
        norm.setseed(j);
        norm.generator();
        rsum = 1 ;
        for ( int i = 0 ; i < n; i ++){
            r[i+1] = r[i] + keppa * ( r_ - r[i]) * delta_t + sigma * sqrt(delta_t) * norm.a[i];
            rsum = rsum *exp (-r[i] * delta_t);
        }

        float bond = max((float)semicoupon_bond(0.03 , r[n-1], 4, 0.25 ) * FV - 985 ,(float) 0.0 );
//        cout << bond << endl;
        option.push_back(bond);

    }
    
    Calculator cal;
    float option_price =  cal.average(option, npath) ;
//    cout << option_price *rsum << endl;
    return option_price*rsum ;

}



float Vasicek::R(float r_start ,float T = 0, float t = 0  ){
   
    float B = 1/ keppa * ( 1 - exp (- keppa * ( T-t )) );
    float A = exp( (r_ - sigma * sigma / (2 * keppa * keppa)) * (B - (T-t)) - sigma * sigma /(4 * keppa) * B * B );
    
    float R =  A * exp ( - B * r_start) ;
    return R;
}

float Vasicek::R_MC(float r_start ,float T = 0, float t = 0 ){
    
    int n = 365 * (T-t) ;
    int npath = 100;
    float rsum = 1 ;
    float delta_t = (T -t) / n ;
    Normal norm (n);
    float r[n];
    r[0] = r_start ;
    
    std::vector<float> r_values;
    
    for (int j = 0 ; j < npath ; j++){
        norm.setseed(j);
        norm.generator();
        rsum = 1 ;
        for ( int i = 0 ; i < n; i ++){
            r[i+1] = r[i] + keppa * ( r_ - r[i]) * delta_t + sigma * sqrt(delta_t) * norm.b[i];
//            cout << exp (-r[i] * delta_t) ;
             rsum = rsum *exp (-r[i] * delta_t);
        }
        r_values.push_back(rsum);
    }
    Calculator cal;
    float ave = cal.average(r_values, npath);
    return ave;
}



float Vasicek::semicoupon_bond(float coupon_rate, float r_start,float T = 0 , float t = 0){
    
    int n = ceil((T-t) * 2 );
    std::vector<float> price;
    float sum = 0 ;
    
    float rterm[n];
    
    for ( int i = 0 ; i < n ; i++) {
        float time = 0.5 * i + 0.5;
        rterm[i] = R_MC(r_start,time,0);
//        cout << rterm[i] <<" " << R(r0,time,0) << endl;
    }
    for ( int i = 0 ; i < n ; i++) {
        sum = sum + rterm[i] * coupon_rate;
    }
    
    sum = sum + R_MC(r_start,T - t,0) * 1;
    return sum;
}



void Q1 () {
    float r0 = 0.05;
    float sigma = 0.18;
    float keppa = 0.82;
    float r_ = 0.05;
    float FV = 1000.0;
    float T = 0.5;
    
//    (a)

    Vasicek Q1(r0, sigma, keppa, r_, FV , T);
    cout << "Q1.a The pure discouting bond price is : " << endl;
    cout << FV * Q1.R_MC(r0 , 0.5, 0) << endl;

////    (b)
    Q1.set_t (4);
    cout << "Q1.b The pure discouting bond price is : " << endl;
    cout << Q1.semicoupon_bond(0.03,r0, 4, 0) * FV  << endl;
    
//    (c)
    Q1.set_t (0.5);
    cout << Q1.option_price1(980, 0.25) << endl;
    
//    (d)
    Q1.set_t (4);
    cout << Q1.option_price2(980, 0.25) << endl;
}

void Q2 (){
    
    float r0 = 0.05;
    float sigma = 0.18;
    float keppa = 0.92;
    float r_ = 0.055;
    float FV = 1000.0;
    float T = 1;
    
    CIR Q2(r0, sigma, keppa, r_, FV , T);
//    cout << Q2.R_MC(0.04, 0.5, 0) << endl;
    
    cout << Q2.option_price2(980, 0.5) << endl;

}

void Q3(){
    
    float phi = 0.03;
    float r0 = 0.03;
    float rho = 0.7 ;
    float a = 0.1;
    float b = 0.3;
    float sigma = 0.03;
    float ita = 0.08;
    
    G2 Q3(phi,r0, rho, a, b, sigma, ita);
    Q3.set_K(985.0);
    Q3.set_mat(0.5);
    Q3.set_T(1.0);
    
//    cout << Q3.R_MC(0.03,1) << endl;
    cout << Q3.option_price() << endl;
    
}


int main() {
//    Q1();
    Q2();
//    Q3();
    return 0;
}

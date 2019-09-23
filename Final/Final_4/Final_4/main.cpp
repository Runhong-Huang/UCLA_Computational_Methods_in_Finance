
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
    
//    float CIR(float r_start, float T, float t);
    
};


class CIR {
    
public:
    
    float r0, alpha, beta, sigma, gamma,FV, T;
    
    
    CIR(float R0, float Alpha, float Beta, float Sigma, float Gamma, float FV, float T): r0(R0), alpha(Alpha), beta(Beta), sigma(Sigma), gamma(Gamma), FV(FV), T(T) {};
    
    
    ~CIR(){};
    
public:
    
    float R_MC(float r_start ,float , float t );
    float option_price2(float k, float mat_t);
};

float CIR::R_MC (float r_start ,float T = 0, float t = 0) {
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
            r[i+1] = r[i] + (alpha + beta * r[i]) * delta_t + sigma * pow(r[i], gamma) * sqrt(delta_t) * norm.b[i];
//            r[i+1] = r[i] + keppa * ( r_ - r[i]) * delta_t + sigma * sqrt(delta_t * r[i] ) * norm.b[i];
            if (r[i+1] < 0){
                cout << "test" << endl;
                r[i+1] = 0;
            }
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
    int npath = 10000;
    float rsum = 1 ;
    float delta_t = mat_t / n ;
    Normal norm (n);
    float r[n];
    r[0] = r0;
    
    std::vector<float> r_values;
    std::vector<float> option;
    
    for (int j = 0 ; j < npath ; j++){
//                cout << j << endl;
        norm.setseed(j);
        norm.generator();
        rsum = 1 ;
        for ( int i = 0 ; i < n; i ++){
            r[i+1] = r[i] + (alpha + beta * r[i]) * delta_t + sigma * pow(r[i], gamma) * sqrt(delta_t) * norm.a[i];
//            r[i+1] = r[i] + keppa * ( r_ - r[i]) * delta_t + sigma * sqrt(delta_t * r[i]) * norm.a[i];
            rsum = rsum *exp (-r[i] * delta_t);
        }
        r_values.push_back(rsum);
        
        float bond = max((float)k - R_MC(r[n], 0.5 ,0) * FV  ,(float) 0.0 );
        //                cout << bond << endl;
        option.push_back( bond * rsum );
    }
    
    Calculator cal;
    float option_price =  cal.average(option, npath) ;
    return option_price;
}



void Q2 (){
    
    float r0 = 0.05;
    float alpha = 0.36;
    float beta = -5.86;
    float sigma = 0.36;
    float gamma = 2.0;
    float FV = 10000.0;
    float T = 1;
    
    
    CIR Q2(r0, alpha, beta, sigma, gamma,FV, T);
    //    cout << Q2.R_MC(0.04, 0.5, 0) << endl;
    
    cout << Q2.option_price2(9800, 0.5) << endl;
    
}




int main() {
    //    Q1();
    Q2();
    //    Q3();
    return 0;
}


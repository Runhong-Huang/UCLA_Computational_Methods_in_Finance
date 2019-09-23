//
//  main.cpp
//  HW3
//
//  Created by HUANG runhong on 1/27/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//


#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;


class Calculator {
    
private:
    
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


// The class used to generate normal random variables
//The class takes three input: seeeds and length and correlation
//The class have two vectors a and b. each contains random variables
class Normal{
private:
    int length = 1000;
    double seed = 52233;
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



class Heston : public Normal{
    
private:
    float ro, r, s0,v0, sigma,alpha, beta, k,T;
    int nsimulation = 10000;
    int length = 1000;
    std::vector<float> price;
    
public:
    Heston (float Rou, float R, float S0, float V0, float Sigma, float Alpha, float Beta,float K, float T ): ro(Rou), r(R), s0(S0), v0(V0),sigma(Sigma), alpha(Alpha), beta(Beta),k(K), T(T) {
    };
    
    ~Heston(){};
    
    float reflection () {
        float dt = T/length;
        float dv, ds;
        std::vector<float> result1;
        
        Normal norm;
        norm.setCorr(ro);
        
        for ( int i = 0; i < nsimulation ; i ++){
            norm.setseed(i);
            norm.generator();
            float v = v0;
            float s = s0;
            for (int j = 0 ; j < length ; j++){
                dv = alpha * (beta - abs(v)) * dt + sigma* sqrt(v)* sqrt(dt) * norm.a[j];
                ds = r * s * dt + sqrt(abs(v)) * s * sqrt(dt) * norm.b[j];
                s = s + ds;
                v = v + dv;
            }
            Calculator cal;
            if (s > k){
                result1.push_back (exp(-r * T) * (s-k));
            }
            else {
                result1.push_back (0);
            }
        }
        Calculator cal;
        return (cal.average(result1, nsimulation));
    }
    
    float partial_truncation () {
        float dt = T/length;
        float dv, ds;
        std::vector<float> result1;
        
        Normal norm;
        norm.setCorr(ro);
        
        for ( int i = 0; i < nsimulation ; i ++){
            norm.setseed(i+21);
            norm.generator();
            float v = v0;
            float s = s0;
            float vplus;
            for (int j = 0 ; j < length ; j++){
                if (v < 0 ){
                    vplus = 0;
                }
                else {
                    vplus = v;
                }
                dv = alpha * (beta - v) * dt + sigma* sqrt(vplus)* sqrt(dt) * norm.a[j];
                ds = r * s * dt + sqrt(vplus) * s * sqrt(dt) * norm.b[j];
                s = s + ds;
                v = v + dv;
            }
            Calculator cal;
            if (s > k){
                result1.push_back (exp(-r * T) * (s-k));
            }
            else {
                result1.push_back (0);
            }
        }
        Calculator cal;
        return (cal.average(result1, nsimulation));
    }
    
    float full_truncation () {
        float dt = T/length;
        float dv, ds;
        std::vector<float> result1;
        
        Normal norm;
        norm.setCorr(ro);
        
        for ( int i = 0; i < nsimulation ; i ++){
            norm.setseed(i+1);
            norm.generator();
            float v = v0;
            float s = s0;
            float s_sum= 0 ;
            float vplus;
            for (int j = 0 ; j < length ; j++){
                if (v < 0 ){
                    vplus = 0;
                    cout << "Haha" << endl;
                }
                else {
                    vplus = v;
                }
//                dv = alpha * (beta - vplus) * dt + sigma* sqrt(vplus)* sqrt(dt) * norm.a[j];
                dv = (alpha + beta * vplus) * dt + sigma* sqrt(vplus)* sqrt(dt) * norm.a[j];
                ds = r * s * dt + sqrt(vplus) * s * sqrt(dt) * norm.b[j];
                s = s + ds;
                v = v + dv;
                s_sum = s_sum + s ;
            }
            
            k = s_sum /length;
//            cout << s_sum / length << endl;
            
            
            Calculator cal;
            if (s > k){
                result1.push_back (exp(-r * T) * (s-k));
            }
            else {
                result1.push_back (0);
            }
        }
        Calculator cal;
        return (cal.average(result1, nsimulation));
    }
    
    
};



int main(){
    

//    Heston (float Rou, float R, float S0, float V0, float Sigma, float Alpha, float Beta,float K, float T ): ro(Rou), r(R), s0(S0), v0(V0),sigma(Sigma), alpha(Alpha), beta(Beta),k(K), T(T) {
//    };
    
    Heston Q4(-0.75, 0.05, 20.0, 0.06, 0.25, 0.45 , -5.105, 50, 2);
//    cout << Q4.reflection() << endl;
//    cout << Q4.partial_truncation() << endl;
    cout << Q4.full_truncation() << endl;
    

    
    return 0;
};


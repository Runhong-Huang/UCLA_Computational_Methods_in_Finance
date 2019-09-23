//
//  main.cpp
//
//  Created by HUANG runhong on 2/20/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
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


class Stockprice : public Normal {
    
private:
    
protected:
    float s0, r , k , sigma , t;
    int n = 10000;
    int exe_step = 1000;
    int **index ;
    float **stock_price = new float *[n];
    
    
public:
    Stockprice (float S0, float R, float K, float Sigma, float T ) : s0{S0}, r(R), k(K), sigma(Sigma),t(T){
    };
    ~Stockprice(){};
    
    void stock_generator();
    void set_s0 (float s);
    void set_time(float time);
    void set_vol (float vol);
    float option();
    
};

void Stockprice::set_vol(float v){
    this-> sigma = v;
}

void Stockprice::set_s0 (float s){
    this -> s0 = s;
}

void Stockprice::set_time (float time){
    this -> t = time;
}

void Stockprice::stock_generator(){
    
    float dt = t/ exe_step;
    //    cout << dt << endl;
    Normal norm(exe_step);
    //    generate the stock price
    for(int i = 0; i < n ; i = i+1){
        stock_price[i] = new float[exe_step+1];
        stock_price[i][0] = s0;
        norm.setseed(i);
        norm.generator();
        for (int j =0 ; j < exe_step ; j++){
            stock_price[i][j+1] = stock_price[i][j] + r * stock_price[i][j] * dt + sigma * stock_price[i][j] * sqrt(dt) * norm.a[j];
        }
    }
    
    std::vector<float> option_price;
    float count = 0 ;
    float countl = 0 ;
    
    for(int i = 0; i < n ; i = i+1){
        for (int j =0 ; j < exe_step ; j++){
            float l = 50 * exp (0.138629 * j * dt);
            float u = 200 - 50 * exp (0.138629 * j * dt);
            
            if (stock_price[i][j] < l ){
                float tmp = k - stock_price[i][j];
                tmp = exp(-r * j * dt) * tmp;
                option_price.push_back(tmp);
//                cout << "l" << j << endl;
                count++;
                countl++;
                break;
            }
            
            if (stock_price[i][j] > u){
                float tmp = stock_price[i][j] - k;
                tmp = exp(-r * j * dt) * tmp;
                option_price.push_back(tmp);
//                cout << "u" << j << endl;
                count++;
                break;
            }
            
        }
    }
    
    Calculator cal;
    cout << cal.average(option_price, count) << endl;
    cout << countl/count   << endl;
        
    
}

float Stockprice::option(){

    return 0;
}



int main() {
    
    Stockprice Q1(100.0, 0.05, 100.0, 0.35, 5.0);
    
    Q1.stock_generator();
    
    return 0;
}


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


class BlackSchole_Numeric:public Normal {
    
private:
    int length = 10000;
    int seed = 9002350;
//    int T = 100;
    float s0,r, sigma,k , T;
    
    
public:
    BlackSchole_Numeric(float S0, float R , float Sigma, float K, float T): s0(S0),r(R),sigma(Sigma), k(K), T(T) {};
    
    ~BlackSchole_Numeric(){};
    
    std::vector<float> vector;
    std::vector<float> stock_price;
    std::vector<float> vector_red1;
    std::vector<float> vector_red2;
    std::vector<float> vector_red;
    
//    void sett(float T){
//        this->T = T;
//    };
    
    void sets0(float S0){
        this->s0 = S0;
    };
    
    void setseed(int s) {
        this-> seed = s;
    }
    
    void setlength(int l){
        this->length = l;
    };
    
    void stock_generator(){
        
        Normal norm(length);
        norm.generator();
        stock_price.clear();
        for (int i = 0 ; i < length ; i++){
            stock_price.push_back(s0 * exp((r- sigma*sigma/2) * T + sigma * sqrt(T) * norm.a[i]));
        }
    };
    
    void generator(){
        float tmp = 0;
        float stock;
        vector.clear();
        Normal norm(length);
        norm.a.clear();
        norm.generator();
//        Calculator cal;
//        cout << cal.average(norm.a, length) << endl;
//        cout << cal.variance(norm.a, length) << endl;
        for (int i = 0 ; i < length ; i++){
            stock = s0 * exp( (r- sigma*sigma/2.0) * T + sigma * sqrt(T) * norm.a[i]);
            if (stock > k){
                tmp = exp(-r * T) * ( stock - k);
            }
            else {
                tmp = 0;
            }
            vector.push_back(tmp);
        }
    };
    
    void generator_red(){
        float tmp1, tmp2, stock1,stock2;
        vector.clear();
        vector_red1.clear();
        vector_red2.clear();
        vector_red.clear();
        Normal norm(length);
        norm.setseed(seed);
        norm.setCorr(-1);
        norm.generator();
        
        for (int i = 0 ; i < length ;  i++){
            stock1 = s0 * exp( (r- sigma*sigma/2.0) * T + sigma * sqrt(T) * norm.a[i]);
            if (stock1 > k){
                tmp1 = exp(-r *T) * (stock1 - k);
            }
            else {
                tmp1 = 0 ;
            }
            
            vector_red1.push_back(tmp1);
        }
        
        for (int i = 0 ; i < length ;  i++){
            stock2 = s0 * exp( (r- sigma*sigma/2.0) * T + sigma * sqrt(T) * norm.b[i]);
            if (stock2 > k){
                tmp2 = exp(-r *T) * (stock2 - k);
            }
            else {
                tmp2 = 0 ;
            }
            vector_red2.push_back(tmp2);
        }
        
        for (int i = 0 ; i < length ;  i++){
            vector_red.push_back(vector_red1[i]/2.0 + vector_red2[i]/2.0);
        
        }
        

        
    };
};


class BlackSchole_Analytic: public Normal{
    
private:
    int length = 1000000;
    float s0,r, sigma,k , T;
    
    
public:
    BlackSchole_Analytic(float S0, float R , float Sigma, float K, float T): s0(S0),r(R),sigma(Sigma), k(K), T(T) {
    };
    
    ~BlackSchole_Analytic(){};
    
    void sets0(float s){
        this-> s0 = s;
    }
    
    template <class T>
    float normal_table_generator(T n){
        float N ;
        float d1 = 0.0498673470;
        float d2 = 0.0211410061;
        float d3 = 0.0032776263;
        float d4 = 0.0000380036;
        float d5 = 0.0000488906;
        float d6 = 0.0000053830;
        
        
        if (n >= 0.0){
            N = 1.0 - 1.0/2.0 * pow((1.0 + d1*n + d2*n*n + d3 * n*n*n + d4*n*n*n*n + d5*n*n*n*n*n + d6*n*n*n*n*n*n), -16);
        }
        
        else if (n < 0.0) {
            n = -n;
            N =  1.0/2.0 * pow((1.0 + d1*n + d2*n*n + d3 * n*n*n + d4*n*n*n*n + d5*n*n*n*n*n + d6*n*n*n*n*n*n), -16) ;
        }
        else {
            N = NULL;
        }
        return N;
    };
    
    
    float call_option(){
        float call_price;
//        float d1 = 1/(sigma*sqrt(T)) * (log(s0/k) + (r + 1.0/2.0 * sigma * sigma)*T) ;
        float d1 = 1/(sigma*sqrt(T)) * (log(s0/k) + (r + 1.0/2.0 * sigma * sigma)*T) ;
        float d2 = 1/(sigma*sqrt(T)) * (log(s0/k) + (r - 1.0/2.0 * sigma * sigma)*T);
        
//        cout << normal_table_generator(d1) << endl;
//        cout << normal_table_generator(d2) << endl;
        
        call_price = s0 * normal_table_generator(d1) - exp(-r * T ) * k * normal_table_generator(d2);
        
        return call_price;

    }
    
    };


float normal_table_generator(float n ){
    float N ;
    float d1 = 0.0498673470;
    float d2 = 0.0211410061;
    float d3 = 0.0032776263;
    float d4 = 0.0000380036;
    float d5 = 0.0000488906;
    float d6 = 0.0000053830;
    
    if (n >= 0.0){
        N = 1.0 - 1.0/2.0 * pow((1.0 + d1*n + d2*n*n + d3 * n*n*n + d4*n*n*n*n + d5*n*n*n*n*n + d6*n*n*n*n*n*n), -16);
    }
    
    else if (n < 0.0) {
        n = -n;
        N =  1.0/2.0 * pow((1.0 + d1*n + d2*n*n + d3 * n*n*n + d4*n*n*n*n + d5*n*n*n*n*n + d6*n*n*n*n*n*n), -16) ;
    }
    else {
        N = NULL;
    }
    return N;
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
            float vplus;
            for (int j = 0 ; j < length ; j++){
                if (v < 0 ){
                    vplus = 0;
                    cout << "Haha" << endl;
                }
                else {
                    vplus = v;
                }
                dv = alpha * (beta - vplus) * dt + sigma* sqrt(vplus)* sqrt(dt) * norm.a[j];
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
    
    
};





float BS_analytic (float s0, float r , float sigma, float k, float T){
    
    float call_price;
    float d1 = 1/(sigma*sqrt(T)) * (log(s0/k) + (r + 1.0/2.0 * sigma * sigma)*T) ;
    float d2 = 1/(sigma*sqrt(T)) * (log(s0/k) + (r - 1.0/2.0 * sigma * sigma)*T);
    
    //        cout << normal_table_generator(d1) << endl;
    //        cout << normal_table_generator(d2) << endl;
    
    call_price = s0 * normal_table_generator(d1) - exp(-r * T ) * k * normal_table_generator(d2);
    
    return call_price;
}

    


void Q1a(int seed, float span){
    int Q1_length = 1000;
    int Q1_nsimulation = 10000;
    
    float y0 = 3.0/4.0;
    std::vector<float> Y2;
    
    float t = 0.0;
    float dt = span/Q1_length;
    float dy = 0;
    int count = 0;

    for (int i = 0; i < Q1_nsimulation ; i++){
        float y = y0;
        t = 0;
        Normal Q1(Q1_length);
        Q1.setseed(seed+i);
        Q1.generator();
        for (int j = 0; j < Q1_length ; j ++){
            t = t + dt;
            dy = ((2.0/(1.0+t))*y + ((1.0+t*t*t)/3.0))* dt + (1+t*t*t)/3.0 * Q1.a[j]*sqrt(dt) ;
            y = y + dy;
        }
        Y2.push_back(y);
        if (y > 5) {
            count++;
        }
    }
    cout << "Q1 The probability of Y>5 is " << (float)count / Q1_nsimulation << endl;
}

void Q1b(int seed){
    //    Q1 b calculate the expextation of x
    int timespan = 2;
    int length = 1000;
    int nsimulation = 10000;
    float dt = (float)timespan/length;
    float dx = 0;
    float x,t;
    std::vector<float> result;
    
    for (int i = 0; i < nsimulation ; i ++){
        x = 1.0;
        t = 0;
        Normal Q1(length);
        Q1.setseed(seed+i);
        Q1.generator();
        for (int j = 0; j < length ; j ++){
            t = t + dt;
            dx = (1.0/5.0 - 1.0/2.0 * x) * dt + 2.0/3.0 * Q1.b[j] * sqrt(dt);
            x = x + dx;
        }
        result.push_back (cbrt(x));
//        cout << x << endl;
    }

    Calculator cal;
    cout << "Q1. E1 is " << cal.average(result, nsimulation) << endl;
}

void Q1c(int seed, float span){
    int Q1_length = 1000;
    int Q1_nsimulation = 10000;
    
    float y;
    std::vector<float> Y;
    
    float timespan = span;
    float t = 0.0;
    float dt = timespan/Q1_length;
    float dy = 0;
    
    for (int i = 0; i < Q1_nsimulation ; i ++){
        y = 3.0/4.0;
        t = 0;
        Normal Q1(Q1_length);
        Q1.setseed(seed+i);
        Q1.generator();
        for (int j = 0; j < Q1_length ; j ++){
            t = t + dt;
            dy = ((2/(1+t))*y + ((1+t*t*t)/3))* dt + (1+t*t*t)/3 * Q1.a[j]*sqrt(dt) ;
            y = y + dy;
        }
        Y.push_back(y);
    }
    //     calculate the probability of Y > 5;
    Calculator cal;
    cout << "Q1. E2 is " << cal.average(Y, Q1_nsimulation) << endl;
}

void Q1d(int seed){
    //    Q1 b calculate the expextation of x
    int timespan = 2;
    int length = 1000;
    int nsimulation = 10000;
    float dt = (float)timespan/length;
    float dx =0; float dy =0;
    float x,y,t;
    std::vector<float> result;
    
    for (int i = 0; i < nsimulation;i++){
        x = 1.0;
        y = 3.0/4.0;
        t = 0;
        Normal Q1(length);
        Q1.setseed(seed+i);
        Q1.setCorr(0);
        Q1.generator();
        for (int j = 0; j < length ; j ++){
            t = t + dt;
            dx = (1.0/5.0 - 1.0/2.0 * x) * dt + 2.0/3.0 * Q1.b[j] * sqrt(dt);
            x = x + dx;
            dy = ((2.0/(1.0+t))*y + ((1.0+t*t*t)/3.0))* dt + (1.0+t*t*t)/3.0 * Q1.a[j]*sqrt(dt);
            y = y + dy;
        }
        if (x > 1.0 ){
//            cout << x * y << endl;
            result.push_back(x * y);
        }
    }
    
//    cout << (int)result.size()/sizeof(result) << endl;
    
    Calculator cal;
    cout << "Q1. E4 is " << cal.average(result, (int)result.size()/sizeof(result)) << endl;
}


void Q2(int seed){
    //    Q1 b calculate the expextation of x
    int timespan = 3.0;
    int length = 1000;
    int nsimulation = 10000;
    float dt = (float)timespan/length;
    float dx =0;
    float x,y,t;
    std::vector<float> result1;
     std::vector<float> result2;
    
    for (int i = 0; i < nsimulation ; i ++){
        x = 1.0;
        y = 3.0/4.0;
        t = 0;
        Normal Q1(length);
        Q1.setseed(seed+i);
        Q1.setCorr(0);
        Q1.generator();
        for (int j = 0; j < length ; j ++){
            t = t + dt;
            dx = 1.0/4.0 * x * dt + 1.0/3.0 * x * sqrt(dt) * Q1.a[j] - 3.0/4.0 * x * sqrt(dt) * Q1.b[j] ;
            x = x + dx;
        }
        result1.push_back(cbrt(1 + x));
        
        y = exp(-0.08 * timespan + 1.0/3.0 * sqrt(timespan) * Q1.a[i] + 3.0/4.0 * sqrt(timespan) * Q1.b[i]);
        result2.push_back(cbrt(1 + y));
    }
    
    Calculator cal;
    cout << "Q2. E1 is " << cal.average(result1, length) << endl;
    cout << "Q2. E2 is " << cal.average(result2, length) << endl;
}


void Q3(int seed, float s0, float r, float sigma, float k, float time){
    
    for (int i = 15 ; i <26; i++ ){
        BlackSchole_Numeric Q3(i,  r, sigma, k, time );
        int length = 10000;
        Q3.setlength(length);
        Q3.generator();
        
        Calculator cal;
        cout << "Q3. numeric method, s0 :" << i << "call option" << cal.average(Q3.vector, length) << endl;
//        cout <<cal.average(Q3.vector, length) << endl;
    }
}

void Q3() {
    //    BlackSchole_Numeric Q4a(15,0.04, 0.25, 20, 0.5);
    std::vector<float> Q3_result1;
    std::vector<float> Q3_result2;
    BlackSchole_Numeric Q3a(15,0.04, 0.25, 20, 0.5);
    Calculator cal;
    Q3a.setseed(124154);
    Q3a.setlength(10000);
    
    for (float s = 15; s< 26 ; s++){
        Q3a.sets0(s);
        Q3a.generator_red();
        Q3_result1.push_back(cal.average(Q3a.vector_red, 10000));
//        cout << cal.average(Q3a.vector_red, 10000) << endl;
    }
    
//    BS method to calculate the call option
    for (float s = 15; s< 26 ; s++){
        cout << "Analytical methods, S0 = " << s << ",  " << BS_analytic(s,0.04, 0.25, 20, 0.5) << endl;
        Q3_result2.push_back(BS_analytic(s,0.04, 0.25, 20, 0.5));
    }
    
//    calculate the greeks
    double delta = 0.01;
    

    
    std::vector<float> Q3c_D;
    std::vector<float> Q3c_T;
    std::vector<float> Q3c_V;
    std::vector<float> Q3c_G;
    std::vector<float> Q3c_R;
    
    for (float s = 15; s< 26 ; s++){
        
        double D = (BS_analytic(s + delta,0.04, 0.25, 20, 0.5) - BS_analytic(s - delta,0.04, 0.25, 20, 0.5)) / (2 * delta) ;
        Q3c_D.push_back(D);
        
        double Dminus = (BS_analytic(s,0.04, 0.25, 20, 0.5) - BS_analytic(s - 2 * delta,0.04, 0.25, 20, 0.5)) / (2 * delta) ;
        
        double Dplus = (BS_analytic(s + 2* delta,0.04, 0.25, 20, 0.5) - BS_analytic(s,0.04, 0.25, 20, 0.5)) / (2 * delta) ;
        
        float G = (Dplus - Dminus) / (2 * delta);
        Q3c_G.push_back(G);
//        cout << G << endl;
        
        float T =(BS_analytic(s ,0.04, 0.25, 20, 0.5 + delta) - BS_analytic(s ,0.04, 0.25, 20, 0.5 - delta)) / (2 * delta) ;
        Q3c_T.push_back(T);
        
        
        float V =(BS_analytic(s ,0.04, 0.25 + delta, 20, 0.5 ) - BS_analytic(s ,0.04, 0.25 - delta, 20, 0.5 )) / (2 * delta) ;
        Q3c_V.push_back(V);
        
        float R =(BS_analytic(s ,0.04 + delta, 0.25, 20, 0.5 ) - BS_analytic(s ,0.04 - delta, 0.25, 20, 0.5 )) / (2 * delta) ;
        Q3c_V.push_back(R);
        
    }
    
}







float *Hatlon_generator (int base, long length) {
    
    
    int n = 1 + ceil(log(length)/log(base));
    
    float vet_base[n];
    float work_vec[n];
    float *result = new float[length];
    
    for (int i =0; i < n; i++){
        vet_base[i] = 0;
        work_vec[i] = 0;
    }
    for(int i = 0; i < n ; i ++ ){
        vet_base[i] = 1/(pow(base,i+1));
    }
    
    for (int i = 0; i < length ; i ++ ){
        int j = 0;
        int condition = 0;
        while (condition == 0){
            work_vec[j] = work_vec[j] + 1;
            if (work_vec[j] < base){
                break;
            }
            else{
                work_vec[j] = 0 ;
                j =j+1;
            }
        }
        
        result[i] = 0;
        for (int j = 0; j < n; j++){
//            cout << result[i] << endl;
            result[i] = result[i] + vet_base[j] * work_vec[j] ;
        }
    }
    return result;
}
    

void Q5(){
    
//    Q5.a
    
    ofstream myfile;
    myfile.open ("Q51.csv");
    float Q5a[100][2] ={};
    for (int i = 0 ; i < 100 ; i++){
        Q5a[i][0] = rand()*1.0/(RAND_MAX*1.0);
        Q5a[i][1] = rand()*1.0/(RAND_MAX*1.0);
        myfile << Q5a[i][0] << "," << Q5a[i][1] << endl;
    }
    myfile.close();
    
//    Q5.b
    float *Q5_base2;
    Q5_base2 = Hatlon_generator(2,100);
    
    float *Q5_base4;
    Q5_base4 = Hatlon_generator(4,100);
    float *Q5_base7;
    Q5_base7 = Hatlon_generator(7,100);
    
    myfile.open ("Q527.csv");
    for (int i = 0 ; i < 100 ; i++){
        myfile << Q5_base2[i] << "," << Q5_base7[i] << endl;
    }
    myfile.close();
    

//    Q5.c
    
    float *Q5c_base2;
    Q5c_base2 = Hatlon_generator(2,10000);
    
    float *Q5c_base4;
    Q5c_base4 = Hatlon_generator(4,10000);
    
    float *Q5c_base7;
    Q5c_base7 = Hatlon_generator(7,10000);
    
    float *Q5c_base5;
    Q5c_base5 = Hatlon_generator(5,10000);
    
    std::vector<float> Q5c1;
    float tmp = 0 ;
    for (int i = 0 ; i < 10000; i++){
        tmp = exp(-Q5c_base2[i] *  Q5c_base4[i])* (sin(6 * M_PI * Q5c_base2[i]) + cbrt((double)cos(2* M_PI * Q5c_base4[i])));
        Q5c1.push_back(tmp);
//        cout << tmp << endl;
    }
    Calculator cal;
    cout << "Q5.c1 if we choose base 2,4 " << cal.average(Q5c1, 10000)<< endl;
    
    std::vector<float> Q5c2;
    for (int i = 0 ; i < 10000; i++){
        tmp = exp(-Q5c_base2[i] *  Q5c_base7[i])* (sin(6 * M_PI * Q5c_base2[i]) + cbrt((double)cos(2* M_PI * Q5c_base7[i])));
        Q5c2.push_back(tmp);
    }
    cout << "Q5.c2 if we choose base 2,7 " << cal.average(Q5c2, 10000)<< endl;
    
    std::vector<float> Q5c3;
    for (int i = 0 ; i < 10000; i++){
        tmp = exp(-Q5c_base5[i] *  Q5c_base7[i])* (sin(6 * M_PI * Q5c_base5[i]) + cbrt((double)cos(2* M_PI * Q5c_base7[i])));
        Q5c3.push_back(tmp);
    }
    cout << "Q5.c3 if we choose base 5,7 " << cal.average(Q5c3, 10000)<< endl;
    
    
}




int main(){
    
// Question 1
//    seed in the bracket
      Q1a(1424,2.0);
    Q1b(1234);
    Q1c(23536,3);
    Q1d(25362);
    
// Question 2
    Q2(41356);
    
// Questin 3
    
    Q3(123215, 15,0.04, 0.25, 20, 0.5);
    Q3();

// Question 4
    Heston Q4(-0.6, 0.03, 48.0, 0.05, 0.42, 5.8 , 0.0625, 50, 0.5);
    cout << Q4.reflection() << endl;
    cout << Q4.partial_truncation() << endl;
    cout << Q4.full_truncation() << endl;
    
//    Question 5
    Q5();
    

    
    return 0;
};

//
//  main.cpp
//  HW2
//
//  Created by HUANG runhong on 1/21/19.
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




class GeometricBrownianMotion:public Normal {
    
private:
    int length = 10000;
    int T = 100;
    float s0,r, sigma,k ;
    
    
public:
    GeometricBrownianMotion(float S0, float R , float Sigma, float K): s0(S0),r(R),sigma(Sigma), k(K){
    }
    
    ~GeometricBrownianMotion(){};
    
    
    std::vector<float> vector;
    std::vector<float> stock_price;
    std::vector<float> vector_red1;
    std::vector<float> vector_red2;
    
    void sett(float T){
        this->T = T;
    };
    
    void sets0(float S0){
        this->s0 = S0;
    };

    
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
        float tmp1 = 0;
        vector.clear();
        Normal norm(length);
        norm.setA(0.0, 1.0);
        norm.a.clear();
        norm.generator();
        for (int i = 0 ; i < length ; i++){
            if (s0 * exp ( (r- sigma*sigma/2) * T + sigma * sqrt(T) * norm.a[i]) > k){
                tmp1 = exp(-r *T) * s0 * exp ( (r- sigma*sigma/2) * T + sigma * sqrt(T) * norm.a[i]) - k;
            }
            else {
                tmp1 = 0 ;
            }
            vector.push_back(tmp1);
        }
    };
    
    void generator_red(){
        float tmp1, tmp2;
        vector.clear();
        Normal norm(length);
        norm.setCorr(-1);
        norm.generator();
        
        for (int i = 0 ; i < length ;  i++){
            if (s0 * exp ( (r- sigma*sigma/2) * T + sigma * sqrt(T) * norm.a[i]) > k){
                tmp1 = exp(-r *T) * s0 * exp ( (r- sigma*sigma/2) * T + sigma * sqrt(T) * norm.a[i]) - k;
            }
            else {
                tmp1 = 0 ;
            }
        
        vector_red1.push_back(tmp1);
        }
        
        for (int i = 0 ; i < length ;  i++){
            if (s0 * exp ( (r- sigma*sigma/2) * T + sigma * sqrt(T) * norm.b[i]) > k){
                tmp2 = exp(-r *T)* s0 * exp ( (r- sigma*sigma/2) * T + sigma * sqrt(T) * norm.b[i]) - k;
            }
            else {
                tmp2 = 0 ;
            }
            
            vector_red2.push_back(tmp2);
        }
    };
    
    
    
};



int main(){
    
//Q1
//    create a normal object with seeds and length as input
    int Q1_length = 1000;
    Normal Q1_normal(Q1_length);
    
    
//    double Q1_seed = 2144;
//    cout << "Please type in seed for Q1: " << endl;
//    cin >> Q1_seed;
//    Q1_normal.setseed(Q1_seed);
    
    
    Q1_normal.setCorr(-0.7);
    Q1_normal.generator();

    
    Calculator cal;
    cout << "Q1. The correlation is: " << cal.corr(Q1_normal.a, Q1_normal.b, Q1_length)<< endl;
    
//  Q2
    int Q2_length = 10000;
    Normal Q2_normal(Q2_length);
    Q2_normal.setCorr(0.6);
    Q2_normal.generator();
    Q2_normal.setseed(32958);
    
    float Q2_sum = 0 ;
    float tmp;
    for (int i=0 ; i < Q2_length ; i++)
    {
        tmp = Q2_normal.a[i] * Q2_normal.a[i] * Q2_normal.a[i] + sin(Q2_normal.b[i]) + Q2_normal.a[i] * Q2_normal.a[i] * Q2_normal.b[i];
        
        Q2_sum = Q2_sum + (tmp>0.0? tmp : 0 );
    }
    
    cout <<"Q2. The expected value is about: "  << Q2_sum  / Q2_length << endl;
    
  
//    Q3 a1
    
    int Q3_length = 5000;
//    create an object Q3.
    Normal Q3(Q3_length);
    Q3.setseed(523);
    Q3.setCorr(0);
    Q3.generator();
    
    std::vector<float> Q3_a1;
    for (int i =0 ; i <Q3_length; i++) {
//        Q3_sum = Q3_sum + Q3.endpoint[i] * Q3.endpoint[i] + sin(Q3.endpoint[i]);
        Q3_a1.push_back(Q3.a[i] * sqrt(5) * Q3.a[i] * sqrt(5) + sin(Q3.a[i] * sqrt(5)));
    }
    

    cout << "Q3.a The expected value for a1 is " << cal.average(Q3_a1,Q3_length) << endl;
    cout << "Q3.a The variance for a1 is " << cal.variance(Q3_a1, Q3_length) << endl;
    
//    Q3 a calculated the expected value
    double Q3_time2 = 0.5;
    double Q3_time3 = 3.2;
    double Q3_time4 = 6.5;
    std::vector<float> Q3_a2;
    std::vector<float> Q3_a3;
    std::vector<float> Q3_a4;
    for (int i =0 ; i < Q3_length; i++) {
        Q3_a2.push_back(exp(Q3_time2/2) * cos(Q3.a[i] * sqrt(Q3_time2)));
        Q3_a3.push_back(exp(Q3_time3/2) * cos(Q3.a[i] * sqrt(Q3_time3)));
        Q3_a4.push_back(exp(Q3_time4/2) * cos(Q3.a[i] * sqrt(Q3_time4)));
    }

    cout << "Q3.a The expected value for a2 is " << cal.average(Q3_a2, Q3_length) << endl;
    cout << "Q3.a The expected value for a3 is " << cal.average(Q3_a3, Q3_length) << endl;
    cout << "Q3.a The expected value for a4 is " << cal.average(Q3_a4, Q3_length) << endl;
    

//    Q3 c1 variance reduction method
    
//    generating exponential distribution vector
    std::vector<float> Q3_c1, Q3_c1_T;
    for (int i=0; i< Q3_length ; i++)
    {
        Q3_c1.push_back(Q3.a[i]* Q3.a[i]*5);
        Q3_c1_T.push_back(Q3_a1[i] - 1 * (Q3_c1[i] - 5));
    }
//    cout << cal.corr(Q3_c1, Q3_a1, Q3_length) * sqrt(cal.variance(Q3_a1, Q3_length)/cal.variance(Q3_c1, Q3_length)) << endl;
    
    
    cout << "Q3.c c1 The VAR Before Reduction is " << cal.variance(Q3_a1, Q3_length) << endl;
    cout << "Q3.c c1 The VAR for Q3 c1 now is " << cal.variance(Q3_c1_T, Q3_length) << endl;


// Q3 c2
    
    std::vector<float> Q3_c2, Q3_c3, Q3_c4;
    
    for (int i=0; i< Q3_length ; i++)
    {
        Q3_c2.push_back(sin(Q3.a[i]* sqrt(Q3_time2)));
        Q3_c3.push_back(sin(Q3.a[i]* sqrt(Q3_time3)));
        Q3_c4.push_back(sin(Q3.a[i]* sqrt(Q3_time4)));
    }
    
//    cout << cal.corr(Q3_c4, Q3_a4, Q3_length) * sqrt(cal.variance(Q3_a4, Q3_length)/cal.variance(Q3_c4, Q3_length)) << endl;
    
    std::vector<float> Q3_c2_T, Q3_c3_T, Q3_c4_T;
    for (int i=0; i< Q3_length ; i++)
    {
        Q3_c2_T.push_back(Q3_a2[i] - 0.01 * (Q3_c2[i]));
        Q3_c3_T.push_back(Q3_a2[i] - 0.1 * (Q3_c2[i]));
        Q3_c4_T.push_back(Q3_a2[i] + 0.18 * (Q3_c2[i]));
    }

    cout << "Q3.c c2 The VAR Before Reduction is " << cal.variance(Q3_a2, Q3_length) << endl;
    cout << "Q3.c c2 The VAR now is " << cal.variance(Q3_c2_T, Q3_length) << endl;
    cout << "Q3.c c3 The VAR Before Reduction is " << cal.variance(Q3_a3, Q3_length) << endl;
    cout << "Q3.c c3 The VAR now is " << cal.variance(Q3_c3_T, Q3_length) << endl;
    cout << "Q3.c c4 The VAR Before Reduction is " << cal.variance(Q3_a4, Q3_length) << endl;
    cout << "Q3.c c4 The VAR now is " << cal.variance(Q3_c4_T, Q3_length) << endl;
    
//    Q4
    int Q4_length  = 10000;
//    float Q4_k = 100;
    GeometricBrownianMotion Q4_c(88, 0.04, 0.2, 100);
    Q4_c.sett(5);
    Q4_c.setlength(Q4_length);
    Q4_c.generator();
    
 
    cout << "Q4. The option price is " <<   cal.average(Q4_c.vector, Q4_length) << endl;
    cout << "Q4. The std of option price is " << cal.variance(Q4_c.vector, Q4_length) << endl;

    
// Q4 use the variance reduction methods atithetic methods
    
    Q4_c.generator_red();
    
    
    std::vector<float> Q4_red ;
    for (int i = 0 ; i < Q4_length ; i++){
        Q4_red.push_back( (Q4_c.vector_red1[i] + Q4_c.vector_red2[i])/2.0) ;
    }
    
    cout << "Q4. We used Anthithetic method to reduce varaice" << endl;
    cout << "Q4. The option price is " <<  cal.average( Q4_red, Q4_length) << endl;
    cout << "Q4. The VAR of random variables is " << cal.variance(Q4_red, Q4_length) << endl;
    
// Q5
    
//    Q5 (a)
    
    int Q5_length = 1000;
    GeometricBrownianMotion Q5_a1(88, 0.04, 0.18, 100);
    Q5_a1.setlength(Q5_length);
//    float Q5_a1_array[11][1000] = {};
    
//    ofstream myfile;
//    myfile.open ("Q5a.csv");
    for (int i = 0 ; i <= 10 ; i++) {
        Q5_a1.sett(i);
        Q5_a1.stock_generator();
//        cout << cal.average(Q5_a1.stock_price, Q4_length) << endl;
    }
    
//    Q5 b
    
    float Q5b[6][1000] = {};
    Normal Q5_norm(1000);
    Q5_norm.setCorr(0);
    
    
    ofstream myfile;
    myfile.open ("Q5b.csv");
    
    for (int i = 0; i < 6 ; i++){
        Q5b[i][0] = 88;
        Q5_norm.setseed(rand());
//        cout << rand() << endl;
        Q5_norm.generator();
        for (int j=0; j <1000 ; j++){
            Q5b[i][j+1] = Q5b[i][j] + Q5b[i][j] * 0.04 * 10.0/1000 + Q5b[i][j] * 0.18 * Q5_norm.a[j] * sqrt(10.0/1000);
            myfile << Q5b[i][j+1]  << "," ;
        }
        myfile << endl;
    }
    
    myfile.close();
    
    
//    Q6
    
//    Q6 a simple dicretiation to calculate the integral
    float Q6_a_sum = 0;
    for (float i = 0; i < 1 ; i = i + 0.0001){
        Q6_a_sum = Q6_a_sum + (sqrt(1 - i*i ) * 0.0001) * 4;
    }
    
    cout << "Q6.a The integral for Q6a is : " << Q6_a_sum << endl;
    
//    Q6 b esitmation from Monte Carlo Method
    std::vector<float> Q6_b_sum;
    float Q6_b_random = 0 ;
    for (int i = 0 ; i < 10000 ; i++){
        Q6_b_random = rand() *1.0/(RAND_MAX*1.0);
        Q6_b_sum.push_back (4 * sqrt(1-Q6_b_random * Q6_b_random));
    }
    
    cout << "Q6.b The integral for Q6b is : " << cal.average(Q6_b_sum, 10000)<< endl;
    cout << "Q6.b The variance is : " << cal.variance(Q6_b_sum, 10000)<< endl;
    
    
//   Q6 c importance sampling
    std::vector<float> Q6_c;
    
    float Q6_length = 10000 ;
    float Q6_c_random1, Q6_c_random2 ;
    float Q6_alpha = 0.74;
    
//    generate h(x) using acceptance& reject methods
    std::vector<float> Q6_tx;
    Q6_tx.clear();
    

    for (int i = 0 ; i < Q6_length ; ){
        Q6_c_random1 = rand() *1.0/(RAND_MAX*1.0);
        Q6_c_random2 = rand() *1.0/(RAND_MAX*1.0);
        if ( Q6_c_random1 < ((1.0- Q6_alpha * Q6_c_random2*Q6_c_random2)/ (1.0 - Q6_alpha/3.0))/2 ) {
            Q6_tx.push_back(Q6_c_random2);
            i++;
        };
        Q6_tx.push_back(Q6_c_random1);
        i++;
    }
    
    for (int i = 0 ; i < Q6_length ; i++){
        Q6_c.push_back (4.0 * sqrt(1.0 - Q6_tx[i]* Q6_tx[i]) / (( 1.0- Q6_alpha * Q6_tx[i]*Q6_tx[i])/ (1.0 - Q6_alpha/3.0)));

    }
    
    cout << "Q6.c The integral for Q6b is : " << cal.average(Q6_c, Q6_length)<< endl;
    cout << "Q6.c The variance is : " << cal.variance(Q6_c, Q6_length)<< endl;
    
  
    return 0;
};












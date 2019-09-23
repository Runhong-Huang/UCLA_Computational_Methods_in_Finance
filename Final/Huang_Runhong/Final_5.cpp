//
//  main.cpp
//  HW6
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



float exponetial( float lamda) {
    float arr = (-1/lamda * log ( rand()*1.0/(RAND_MAX*1.0)));
    return arr;
}


float final_5(){
    
    float s0 = 6000.0;
    float e0 = 0.0096;
    float r = 0.05;
    float q = 0.0;
    float rf = 0.04;
    float sigma1 = 0.10;
    float sigma2 = 0.15;
    float gamma = -0.04;
    float lambda = 1.5;
    float k = 60.0;
    float T = 1.0 ;
    
    
    int npath = 10000;
    int ntime = 600;
    
    float **s = new float *[npath];
    float **e = new float *[npath];
    float **jump_matrix = new float *[npath];
    
    
    //    generate two matrix, one is the value of the jump diffusion, the other is the jump matrix
    for (int i = 0 ; i < npath; i++){
        s[i] = new float [ntime+1];
        e[i] = new float [ntime+1];
        jump_matrix[i] = new float[ntime+1];
        s[i][0] = s0;
        e[i][0] = e0;
    }
    
//    initialize the jump matrix
    for (int i = 0 ; i < npath; i++){
        for (int j = 0 ; j < ntime+1; j++){
            jump_matrix[i][j] = 0;
        }
    }
    
    // generate the jump matrix, where each cell represents a jump incident
    for (int i = 0; i < npath; i++){
        float tmp = exponetial(lambda);
        while ( tmp < T){
            int colnumber =  (tmp * ntime/T) ;
            jump_matrix[i][colnumber]++;
            tmp = tmp + exponetial(lambda);
        }
    }
    
    // generate the matrix of collaterl, using name of v
    Normal norm(ntime+1);
    norm.setCorr(-0.25);
    float dt = T/ntime ;
    
    std::vector<float> sum;
    
    for (int i =0 ; i < npath ; i++){
        norm.setseed(i);
        norm.generator();
        for (int j = 0; j < ntime; j ++){
            
            s[i][j+1] = s[i][j] + s[i][j] * (r-q) * dt + s[i][j] * sigma1 * sqrt(dt) * norm.a[j];
            while (jump_matrix[i][j+1] > 0) {
                s[i][j+1] = s[i][j+1] * (1 + gamma);
                jump_matrix[i][j+1]--;
            }
            
            e[i][j+1] = e[i][j] + e[i][j] * (r-rf) * dt + sigma2 * e[i][j] * norm.b[i];
        }
        
//        float tmp = max((float)s[i][ntime] * e[i][ntime] - k, (float) 0.0 );
        
        float tmp = 0;
        float w = 0 ;
        tmp = s[i][ntime-1] * e[i][ntime-1] ;
        cout << tmp -k << endl ;
        if (tmp > 0){
            sum.push_back(tmp);
            w = w + tmp;
        }
        else {sum.push_back(0.0);};
        
    }
    
    Calculator cal;
    cout << cal.average(sum, npath) << endl;
//
    cout << w /npath << endl;
    sum.clear();
    return 0;
}
    



int main() {
    
    
    final_5();
    
    return 0;
}


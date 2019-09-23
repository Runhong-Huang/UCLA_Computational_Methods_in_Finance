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


class Stockprice : public Normal {
    
private:
    
protected:
    float s0, r , k , sigma , t;
    int n = 100000;
    int exe_step = 600;
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
    norm.setCorr(-1.0);
    
    //    float ** stock_price = new float *[n];
    
    //    initializing the matrix : stock matrix, index matrix and ECV ( the expected continuation value)
    
    //    generate the stock price
    for(int i = 0; i < n ; i = i+2){
        stock_price[i] = new float[exe_step+1];
        stock_price[i+1] = new float[exe_step+1];
        stock_price[i][0] = s0;
        stock_price[i+1][0] = s0;
        norm.setseed(i);
        norm.generator();
        for (int j =0 ; j < exe_step ; j++){
            stock_price[i][j+1] = stock_price[i][j] + r * stock_price[i][j] * dt + sigma * stock_price[i][j] * sqrt(dt) * norm.a[j];
            stock_price[i+1][j+1] = stock_price[i+1][j] + r * stock_price[i+1][j] * dt + sigma * stock_price[i+1][j]* sqrt(dt) * norm.b[j];
        }
    }
}


class Look_back_option : public Stockprice {
    
public:

    Look_back_option(float S0, float R, float K, float Sigma, float T )  :  Stockprice (S0, R,  K, Sigma, T ){};
    
    float call_option_calculator();
    float put_option_calculator();
};

float Look_back_option::put_option_calculator(){
    
    stock_generator();
    
    std::vector<float> payoff;
    
    for (int i =0; i < n ; i ++){
        float pmin = stock_price[i][0];
        for (int j = 0 ; j < exe_step ; j ++){
            if (stock_price[i][j] < pmin) pmin = stock_price[i][j];
        }
        
        float value = max((float) (k - pmin), (float)0.0);
        payoff.push_back(value);
        
    }
    
    Calculator cal;
    cout << cal.average(payoff, n) << endl;
    return 0;
    
}


float Look_back_option::call_option_calculator(){
    
    stock_generator();
    
    std::vector<float> payoff;
    
    
    for (int i =0; i < n ; i ++){
        float pmax = stock_price[i][0];
        for (int j = 0 ; j < exe_step ; j ++){
            if (stock_price[i][j] > pmax) pmax = stock_price[i][j];
        }
        
        float value = max((float) (pmax-k), (float)0.0);
        payoff.push_back(value);
       
    }
    
    Calculator cal;
    cout << cal.average(payoff, n) << endl;
    return 0;
    
}

void Q1 () {
    
    Look_back_option Q1(98.0 , 0.03 , 100, 0.06, 1);
    Q1.setseed(1345);
    
    cout << "Calculate the look back call option" << endl;
    for (int i = 0 ; i < 10 ; i++){
        float tmpvol =  0.12 + 0.04 * i;
        Q1.set_vol(tmpvol);
        Q1.call_option_calculator();
    }
//
    cout << "Calculate the look back put option" << endl;
    for (int i = 0 ; i < 10 ; i++){
        float tmpvol =  0.12 + 0.04 * i;
        Q1.set_vol(tmpvol);
        Q1.put_option_calculator();
    }
}


float exponetial( float lamda) {
    float arr = (-1/lamda * log ( rand()*1.0/(RAND_MAX*1.0)));
    return arr;
}

//float *exponetial( float lamda, int length = 100000) {
//    srand((int)time(0));
//    float *arr = new float [length];
//    for (int i=0; i< length ; i++)
//    {
//        arr[i] = (-1/lamda * log ( rand()*1.0/(RAND_MAX*1.0)));
//    }
//    return arr;
//}

float *Proj6_2function(float lambda1 = 0.2 , float lambda2 = 0.4 , float T = 5.0){
    
    float v0 = 20000;
    float L0 = 22000;
    float miu = -0.1;
    float sigma = 0.2;
    float gamma = -0.4;
    float r0 = 0.02;
    float delta = 0.25;
    float alpha = 0.7;
    float epsilon = 0.95;

    
    int npath = 10000;
    int ntime = 600;
    
    float **v = new float *[npath];
    float **jump_matrix = new float *[npath];
    
    
//    generate two matrix, one is the value of the jump diffusion, the other is the jump matrix
    for (int i = 0 ; i < npath; i++){
        v[i] = new float [ntime+1];
        jump_matrix[i] = new float[ntime+1];
        v[i][0] = v0;
    }
    
    for (int i = 0 ; i < npath; i++){
        for (int j = 0 ; j < ntime+1; j++){
            jump_matrix[i][j] = 0;
        }
    }
    
// generate the jump matrix, where each cell represents a jump incident
    for (int i = 0; i < npath; i++){
        float tmp = exponetial(lambda1);
        while ( tmp < T){
            int colnumber =  (tmp * ntime/T) ;
            jump_matrix[i][colnumber]++;
            tmp = tmp + exponetial(lambda1);
        }
    }
    
// generate the matrix of collaterl, using name of v
    Normal norm(ntime+1);
    float dt = T/ntime ;
    for (int i =0 ; i < npath ; i++){
        norm.setseed(i);
        norm.generator();
        for (int j = 0; j < ntime; j ++){
            v[i][j+1] = v[i][j] + v[i][j] * miu * dt + v[i][j]* sigma * sqrt(dt) * norm.a[j];
            while (jump_matrix[i][j+1] > 0) {
//                cout << i << " " << j  << " " << jump_matrix[i][j+1] << endl;
                v[i][j+1] = v[i][j+1] * (1 + gamma);
                jump_matrix[i][j+1]--;
            }
        }
    }
    
    //    calculate the necessary parameters
    float R  = r0 + delta * lambda2;
//    cout << R << endl;
    float r = R /12.0 ;
    float n = T * 12.0 ;
    float pmt = L0 * r /(1 - 1/pow((1+r),n));
   
    float a = pmt /r ;
    float b = pmt /( r * (pow((1+r), n)));
    float c = ( 1 + r );

    
    //    generate the matrix of L
    float lt[ntime+1];
    for (int j = 0 ; j < ntime +1 ; j++){
        lt[j] = a - b * pow(c, 12.0 * j/ntime * T );
    }
    
//    generate Qt
    float qt[ntime+1];
    float beta = (epsilon - alpha) / T ;
    for (int j = 0 ; j < ntime + 1 ; j++){
        qt[j] = alpha + beta * j/ntime * T ;
    }
    
    
//    generate matrix Q of the stopping time Vt < qt * Lt
//    generate matrix S of the stopping time with poisson process lambda2.
    float **Q = new float *[npath];
    float **S = new float *[npath];
    for (int i = 0 ; i < npath ; i++) {
        Q[i] = new float [ntime + 1];
        S[i] = new float [ntime + 1];
    }
    
    //    initialize the Q matrix and S matrix
    for (int i = 0; i < npath ; i++ ){
        for (int j = 0 ; j < ntime + 1 ; j++){
            Q[i][j] = 0 ;
            S[i][j] = 0 ;
        }
    }
    
    
//    generate the S matrix
    for (int i = 0; i < npath; i++){
        float tmp = exponetial(lambda2);
        while ( tmp < T){
            int colnumber =  (tmp * ntime/T);
            S[i][colnumber]++;
            tmp = tmp + exponetial(lambda2);
        }
    }
    
//    generate the Q matrix
    for (int i = 0; i < npath ; i++ ){
        for (int j = 0 ; j < ntime + 1 ; j++){
            if (v[i][j] < qt[j] * lt[j] ){
                Q[i][j] = 1 ;
            }
        }
    }

    
//    generate the payoff of the simulation
    std::vector<float> payoff;
    std::vector<float> exe;

    
    int count = 0;
    for (int i = 0; i < npath ; i++ ){
        for (int j = 0 ; j < ntime + 1 ; j++){
            if ((Q[i][j] == 1)  && (S[i][j] == 1) ){
                float tmp1 = lt[j] - epsilon * v[i][j];
                tmp1 = max ((float) tmp1 , (float) 0.0);
                float tmp2 = lt[j] - epsilon * v[i][j];
                tmp2 = abs(tmp2);
                float tmp = max((float) tmp1 , (float) tmp2);
                tmp = tmp / pow((1+r0/12), (int) ( ceil(j /10) ));
                
                payoff.push_back(tmp);
                exe.push_back(j);
                count ++;
                break;
            }
            else if (Q[i][j] == 1){
                float tmp1 = lt[j] - epsilon * v[i][j];
                tmp1 = max ((float) tmp1 , (float) 0.0);
                tmp1 = tmp1 / pow((1+r0/12), (int) ( ceil(j /10) ));
                payoff.push_back(tmp1);
                exe.push_back(j);
                count ++;
                break;
            }
            else if (S[i][j] == 1){
                float tmp2 = lt[j] - epsilon * v[i][j];
                tmp2 = abs(tmp2);
                tmp2 = tmp2 / pow((1+r0/12), (int) ( floor(j /10) ));
                payoff.push_back(tmp2);
                exe.push_back(j);
                count++;
                break;
            }
        }
    }
    
    float Prob = (float) count / npath;


    float sum_payoff = 0 ;
    float exe_time = 0 ;
    for (int i = 0 ; i < count ; i ++){
        sum_payoff = sum_payoff + payoff[i];
        exe_time = exe_time + exe[i];
    }
    
    exe_time = exe_time /10;
    
    float *result = new float[3];
    
    result[0] = sum_payoff/npath;
    result[1] = Prob;
    
    result[2] = exe_time /count;
    result[2] = exe_time /npath;
    
    exe.clear();
    
    return result;
}



void Q2 (){
    float *ptr = Proj6_2function();
    cout << ptr[0] << endl;
    cout << ptr[1] << endl;
    cout << ptr[2] << endl;
    
//    generate data for graph a
    float l1 = 0.2 ;
    for (int i =0 ; i < 9 ; i ++){
        for (int t = 3 ; t < 9 ; t++){
            float l2 = i * 0.1  ;
            ptr = Proj6_2function(l1, l2 , t);
            cout << ptr[2] << ",";
        }
    }
    cout << endl;
    cout << endl;
    
    float l2 = 0.4 ;
    for (int i =0 ; i < 9 ; i ++){
        for (int t = 3 ; t < 9 ; t++){
            float l1 = i * 0.05 ;
            ptr = Proj6_2function(l1, l2 , t);
            cout << ptr[2] << ",";
        }
    }
    
}

int main() {
    Q1();
    Q2();
    
    
    return 0;
}

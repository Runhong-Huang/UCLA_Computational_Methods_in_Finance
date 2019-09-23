//
//  main.cpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "Binomial.hpp"

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

class Binomial {
    
private:
    
    float s0, r ,sigma, k, T;

    
    float u, d , p;


public:
    
    int nlayers;

    Binomial(float S0, float R , float Sigma, float K, float T, int N): s0(S0),r(R),sigma(Sigma), k(K), T(T),nlayers(N){};
    
    ~Binomial(){};
    
    void setn(int N);
    void sets0(float s);
    void set_t(float t);
    void set_sigma(float sigma);
    void set_r (float r);


    void SetUDP1();
    void SetUDP2();
    void SetUDP3();
    void SetUDP4();

    float call_pricing();
    float put_pricing();
};

void Binomial::setn(int N){
    this -> nlayers = N;
}

void Binomial::sets0(float s){
    this -> s0 = s;
}

void Binomial::set_t(float t){
    this-> T = t;
}

void Binomial::set_sigma(float s){
    this-> sigma = s;
}

void Binomial::set_r(float s){
    this-> r = s;
}


void Binomial::SetUDP1(){
    float dt = T/nlayers;
    
    float c = 1.0/2.0 * ( exp (-r * dt) + exp((r + sigma* sigma) * dt));

    this -> d = c - sqrt(c*c -1) ;
    this -> u = 1/d ;
    this -> p = (exp(r * dt) -d )/(u -d);

//    cout << u << " " << d << " "<< p << endl;
}



void Binomial::SetUDP2(){
    
    float dt = T/nlayers;
    
    this -> u = exp(r * dt) *( 1 + sqrt(exp(sigma*sigma*dt)-1)) ;
    this -> d = exp(r * dt) *( 1 - sqrt(exp(sigma*sigma*dt)-1));
    this -> p = 0.5;
    
//    cout << u << " " << d << " "<< p << endl;
}

void Binomial::SetUDP3(){
    float dt = T/nlayers;
    
    this -> u = exp((r - 0.5 * sigma* sigma)* dt + sigma * sqrt(dt));
    this -> d = exp((r - 0.5 * sigma* sigma)* dt - sigma * sqrt(dt));
    this -> p = 0.5;
    
//    cout << u << " " << d << " "<< p << endl;
}

void Binomial::SetUDP4(){
    float dt = T/nlayers;
    
    this -> u = exp(sigma*sqrt(dt));
    this -> d = exp(-sigma*sqrt(dt));
    this -> p = 0.5 + 0.5 * (((r- 0.5 * sigma * sigma) * sqrt(dt))/(sigma));
    
//    cout << u << " " << d << " "<< p << endl;
}

float Binomial::call_pricing () {

    int n = nlayers;
 
    float dt = T/n;
    
//    float stock_price[n+1][n+1];
//    float call[n+1][n+1];
    
    float** stock_price = new float*[n+1];
    float** call = new float*[n+1];
    for(int i = 0; i < n+1; i++){
        stock_price[i] = new float[n+1];
        call[i] = new float[n+1];
    }
    
//    set the initial price as so
    stock_price[0][0] = s0;
    
//    generate the stock pricing
    for (int i = 0; i < n ; i++){
        for (int j = 0; j < i+1 ; j++){
            stock_price[j][i+1] = stock_price[j][i] * u;
        }
        stock_price[i+1][i+1] = stock_price[i][i] * d ;
    }
    
//    calculate the call pricing
    for (int j = 0; j < n + 1 ; j++){
        call[j][n] =  std::max((float)(stock_price[j][n] -k), float(0.0));
    }
    
    
    for (int i = n; i > 0 ; i--){
        for (int j = 0; j < i ; j++){
            call[j][i-1] = exp(-r * dt) * (call[j][i] * (p) + call[j+1][i] * (1.0 - p));
        }
    }

    return (call[0][0]);
}


float Binomial::put_pricing () {

    int n = nlayers;
    float dt = T/n;

    //    float stock_price[n+1][n+1];
    //    float call[n+1][n+1];

    float** stock_price = new float*[n+1];
    float** put = new float*[n+1];
    for(int i = 0; i < n+1; i++){
        stock_price[i] = new float[n+1];
        put[i] = new float[n+1];
    }

    //    set the initial price as so
    stock_price[0][0] = s0;

    //    generate the stock pricing
    for (int i = 0; i < n ; i++){
        for (int j = 0; j < i+1 ; j++){
            stock_price[j][i+1] = stock_price[j][i] * u;
        }
        stock_price[i+1][i+1] = stock_price[i][i] * d ;
    }

    //    calculate the call pricing
    for (int j = 0; j < n + 1 ; j++){
        put[j][n] =  std::max((float)(k - stock_price[j][n]), float(0.0));
//        cout << put[j][n] << endl;
    }

    for (int i = n; i > 0 ; i--){
        for (int j = 0; j < i ; j++){
            float tmp1 = exp(-r * dt) * (put[j][i] * (p) + put[j+1][i] * (1.0 - p));
            float tmp2 = k - stock_price[j][i-1];
//            float tmp = max(tmp1,tmp2);
            float tmp = tmp1;
            if (tmp > 0){
                put[j][i-1] = tmp;
            }
            else {
                put[j][i-1] = 0;
            }
                
        }
    }

    return (put[0][0]);
}



class Trinomial{
    
private:
    
    float s0, r ,sigma, k, T;
    
    
    float u, d , pu, pd, pm;
    
    float up_delta, down_delta;
    
    
public:
    
    int nlayers;
    
    Trinomial(float S0, float R , float Sigma, float K, float T, int N): s0(S0),r(R),sigma(Sigma), k(K), T(T),nlayers(N){};
    
    ~Trinomial(){};
    
    float tri_call_pricing();
    float tri_logcall_pricing();
    void set_UDP1();
    void set_UDP2();
    void setn(int N);
 
};

void Trinomial::setn(int N){
    this -> nlayers = N;
}

void Trinomial::set_UDP1(){
    float dt = T/nlayers;
    this -> u = exp ( 1.0 * sigma * sqrt(3.0 * dt) );
    this -> d = 1/u;
    this -> pd = ((r * dt * (1-u)) + (r * dt * r * dt) + (sigma* sigma * dt) )/((u-d)* (1.0-d));
    this -> pu = ((r * dt * (1-d)) + (r * dt * r * dt) + (sigma* sigma * dt))/((u-d)*(u-1.0));
    this -> pm = 1-pd -pu;
//    cout << u << " " << d << " "<< pu << endl;
}

void Trinomial::set_UDP2(){
    float dt = T/nlayers;
    
    
//    this -> u = exp ( 1.0 * sigma * sqrt(3.0 * dt) );
//    this -> d = 1/u;
    this -> up_delta = sigma * sqrt(3 * dt);
    this -> down_delta = -sigma * sqrt(3 * dt);
    this -> pd = 0.5 * ( (((sigma * sigma * dt) + (r - 0.5 * sigma * sigma)*(r - 0.5 * sigma * sigma ) * dt * dt )/(up_delta* up_delta)) + (((r - 0.5 * sigma*sigma)*dt)/(down_delta)) ) ;;
    this -> pu = 0.5 * ( (((sigma * sigma * dt) + (r - 0.5 * sigma * sigma)*(r - 0.5 * sigma * sigma ) * dt * dt )/(up_delta* up_delta)) - (((r - 0.5 * sigma*sigma)*dt)/(down_delta)) ) ;
    
    this -> pm = 1-pd -pu;
//    cout << up_delta << " " << pu << " "<< pd << endl;
}


float Trinomial::tri_call_pricing(){
    
    int n = nlayers;
    float dt = T/n;
    
    //    float stock_price[n+1][n+1];
    //    float call[n+1][n+1];
    
    float** stock_price = new float*[2*n+1];
    float** call = new float*[2*n+1];
    for(int i = 0; i < 2*n+1; i++){
        stock_price[i] = new float[n+1];
        call[i] = new float[n+1];
    }
    
    //    set the initial price as so
    stock_price[0][0] = s0;
    
    //    generate the stock pricing
    for (int i = 0; i < n ; i++){
        for (int j = 0; j < 2*i+1 ; j++){
            stock_price[j][i+1] = stock_price[j][i] * u ;
//            cout <<stock_price[j][i+1] << endl;
        }
        stock_price[2*i+1][i+1] = stock_price[2*i][i] ;
        stock_price[2*i+2][i+1] = stock_price[2*i][i] * d ;
    }
    
//    for (int i = 0; i < 2*n+1 ; i++){
//        cout << stock_price[i][n] << endl;
//    }
//
    
//    calculate the call pricing
    for (int j = 0; j < 2 * n + 1 ; j++){
        call[j][n] =  std::max((float)(stock_price[j][n] - k), float(0.0));
    }
    

//
//
    for (int i = n; i > 0 ; i--){
        for (int j = 0; j < 2*i -1 ; j++){
            call[j][i-1] = exp(-r * dt) * (call[j][i] * (pu) + call[j+1][i] * (pm) + call[j+2][i] * (pd));
        }
    }
    
    return (call[0][0]);
}


float Trinomial::tri_logcall_pricing(){

    int n = nlayers;
    float dt = T/n;

    //    float stock_price[n+1][n+1];
    //    float call[n+1][n+1];

    float** stock_price = new float*[2*n+1];
    float** call = new float*[2*n+1];
    for(int i = 0; i < 2*n+1; i++){
        stock_price[i] = new float[n+1];
        call[i] = new float[n+1];
    }

    //    set the initial price as so
    stock_price[0][0] = s0;

    //    generate the stock pricing
    for (int i = 0; i < n ; i++){
        for (int j = 0; j < 2*i+1 ; j++){
            stock_price[j][i+1] = stock_price[j][i] * exp( up_delta) ;
            //            cout <<stock_price[j][i+1] << endl;
        }
        stock_price[2*i+1][i+1] = stock_price[2*i][i]  ;
        stock_price[2*i+2][i+1] = stock_price[2*i][i] *exp(down_delta) ;
    }



    //    calculate the call pricing
    for (int j = 0; j < 2*n + 1 ; j++){
        call[j][n] =  std::max((float)(stock_price[j][n] -k), float(0.0));
    }


    for (int i = n; i > 0 ; i--){
        for (int j = 0; j < 2*i -1 ; j++){
            call[j][i-1] = exp(-r * dt) * (call[j][i] * (pu) + call[j+1][i] * (pm) + call[j+2][i] * (pd));
        }
    }

    return (call[0][0]);
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


void Q1(){
    
    Binomial Q1(32, 0.05, 0.24, 30, 0.5, 1); //S0, R ,sigma, K , T ,N
 
    
    vector<int> nrange = { 10, 20, 40, 80, 100, 200,500 };
    
//    Method1
    cout << "Q1, Method 1" << endl;
    for (auto itr:nrange){
        Q1.setn(itr);
        Q1.SetUDP1();
        cout << Q1.call_pricing() << "," ;
    }
    cout << endl;
    
//    Method2
    cout << "Q1, Method 2" << endl;
    for (auto itr:nrange){
        Q1.setn(itr);
        Q1.SetUDP2();
        cout << Q1.call_pricing() << ",";
    }
    cout << endl;
    
//    Method3
    cout << "Q1, Method 3" << endl;
    for (auto itr:nrange){
        Q1.setn(itr);
        Q1.SetUDP3();
        cout << Q1.call_pricing() << ",";
    }
    cout << endl;
    
//    Method4
    cout << "Q1, Method 4" << endl;
    for (auto itr:nrange){
        Q1.setn(itr);
        Q1.SetUDP4();
        cout << Q1.call_pricing() << ",";
    }
    cout << endl;
}


void Q2() {
    cout << "Question 2 " << endl;
    
    float S0 = 1110.75;
    float r = 0.02;
    float k = 1220;
    float vol = 0.2331;
    
    
    Binomial Q2(S0, r, vol, k, 1, 1000); //S0, R ,sigma, K , T ,N
    Q2.SetUDP1();
    cout << Q2.call_pricing() << endl;
    
    
//    73.80 is the midium price of bid and ask
    float newvol =0.24260;
    
    Binomial Q2b(S0, r, newvol, k, 1, 1000); //S0, R ,sigma, K , T ,N
    Q2b.setn(1000);
    Q2b.SetUDP1();
    cout << Q2b.call_pricing() << endl;
}


void Q3(){
    
    cout << "Question 3 " << endl;
    
    float S0 = 49.0;
    float r = 0.03;
    float k = 50.0;
    float sigma = 0.2;
    float T = 0.3846;
    float miu = 0.14;
    float delta = 0.1;
    
    
    Binomial Q3(S0, r, sigma, k, T, 1000); //S0, R ,sigma, K , T ,N
    
//    (i)
//    for (float s = 20; s <= 80 ; s= s+2){
//        Q3.sets0(s + delta);
//        Q3.SetUDP1();
//        float plus = Q3.call_pricing();
//        Q3.sets0(s - delta);
//        Q3.SetUDP1();
//        float minus = Q3.call_pricing();
//        cout <<(plus-minus)/(2.0 * delta) << ",";
//    }
    
    
//    (ii)
    Q3.sets0(49.0);
    for (float t = 0; t <= T ; t = t + 0.01){
        Q3.set_t(t);
        Q3.sets0(S0 + delta);
        Q3.SetUDP1();
        float plus = Q3.call_pricing();
        Q3.sets0(S0 - delta);
        Q3.SetUDP1();
        float minus = Q3.call_pricing();
        cout <<(plus-minus)/(2.0 * delta) << ",";
    }
    
    cout << endl;
    
//    iii theta
//    Q3.set_t(T);
//    for (float s = 20; s <= 80 ; s= s+2){
//        Q3.sets0(s);
//        Q3.set_t(T + delta);
//        Q3.SetUDP1();
//        float plus = Q3.call_pricing();
//        Q3.set_t(T - delta);
//        Q3.SetUDP1();
//        float minus = Q3.call_pricing();
//        cout <<(plus-minus)/(2.0 * delta) << ",";
//    }
//
//    iv Gamma
//    Q3.set_t(T);
//    Q3.sets0(S0);
//    delta = 0.2;
//
//    for (float s = 20; s <= 80 ; s= s+2){
//
//        Q3.sets0(s);
//        Q3.SetUDP1();
//        float zero = Q3.call_pricing();
//
//        Q3.sets0(s + 2.0 * delta);
//        Q3.SetUDP1();
//        float plus2 = Q3.call_pricing();
//
//        Q3.sets0(s - 2.0 * delta);
//        Q3.SetUDP1();
//        float minus2 = Q3.call_pricing();
//
//        float Detalminus = (plus2 - zero) / (2 * delta);
//        float Detalplus = (zero - minus2) / (2 * delta);
//        float Gamma = (Detalplus - Detalminus) /2.0;
//        float Gamma = (plus2 - minus2) / (4 * delta );
//
//        cout <<Gamma << ",";
//    }
    
//    Vega
    
//    Q3.set_t(T);
//    Q3.sets0(S0);
//
//    for (float s = 20; s <= 80 ; s= s+2){
//        Q3.sets0(s);
//        Q3.set_sigma(sigma + delta);
//        Q3.SetUDP1();
//        float plus = Q3.call_pricing();
//        Q3.set_sigma(sigma - delta);
//        Q3.SetUDP1();
//        float minus = Q3.call_pricing();
//        float vega = (plus - minus) / (2 * delta);
//        cout <<vega << ",";
//    }
    
    
    
//    Rho
//    Q3.set_sigma(sigma);
//    Q3.sets0(S0);
//
//    for (float s = 20; s <= 80 ; s= s+2){
//        Q3.sets0(s);
//        Q3.set_r(r + delta);
//        Q3.SetUDP1();
//        float plus = Q3.call_pricing();
//        Q3.set_r(r - delta);
//        Q3.SetUDP1();
//        float minus = Q3.call_pricing();
//        float rho = (plus - minus) / (2 * delta);
//        cout <<rho << ",";
//    }

   
}

void Q4(){
    float S0 = 100.0;
    float r = 0.05;
    float k = 100.0;
    float sigma = 0.3;
    float T = 1;
    Binomial Q3(S0, r, sigma, k, T, 100);
    
    for (int i =80; i <= 120 ; i = i+4){
        Q3.sets0(i);
        Q3.SetUDP1();
        cout << Q3.put_pricing() << "," ;
    }
}

void Q5(){
    float S0 = 32.0;
    float r = 0.05;
    float k = 30.0;
    float sigma = 0.24;
    float T = 0.5;
    
    Trinomial Q5(S0, r, sigma, k, T, 500);
    
    vector<int> nrange = { 10, 15, 20, 40,70, 80, 100, 200,500 };
    
    for (auto itr:nrange){
        Q5.setn(itr);
        Q5.set_UDP1();
        cout << Q5.tri_call_pricing() << "," ;
    }
    
    for (auto itr:nrange){
        Q5.setn(itr);
        Q5.set_UDP2();
        cout << Q5.tri_logcall_pricing() << "," ;
    }
    

   
    
};

void Q6(){
//    Get input
    float S0, k, T, r, sigma;
    int N, base1 ,base2;
    
    cout << endl;

    cout << "Please type in initial price s0 :" << endl;
    cin >> S0 ;
    cout << "Please type in strike price K :" << endl ;
    cin >> k ;
    cout << "Please type in Maturity time :" << endl ;
    cin >> T ;
    cout << "Please type in interests rate :" << endl ;
    cin >> r ;
    cout << "Please type in voltility :" << endl ;
    cin >> sigma;
    cout << "Please type in Numer of points :" << endl ;
    cin >> N;
    cout << "Please type in Base1 for Halton number :" << endl ;
    cin >> base1;
    cout << "Please type in Base2 for Halton number :" << endl ;
    cin >> base2;
    
//    float S0 = 32.0;
//    float r = 0.05;
//    float k = 30.0;
//    float sigma = 0.24;
//    float T = 0.5;
//    int N = 100000;
//
//    int base1 = 3;
//    int base2 = 7;
    
    float *H1 , *H2;
    H1 = Hatlon_generator(base1, N);
    H2 = Hatlon_generator(base2, N);
    
    std::vector<float> call;
    std::vector<float> Z1;
    std::vector<float> Z2;
    
    for (int i = 0 ; i < N ;i++ ){
        float z1 = sqrt(-2.0 * log (H1[i])) * cos(2* M_PI * H2[i]);
        float z2 = sqrt(-2.0 * log (H1[i])) * sin(2 * M_PI * H2[i]);
        Z1.push_back(z1);
        Z2.push_back(z2);
        float tmp1 = S0 * exp ((r - 0.5 * sigma * sigma)* T + sigma * sqrt(T) * z2) - k;
        if (tmp1 > 0 ){
            call.push_back(tmp1 * exp(-r * T));
        }
        if (tmp1 < 0){
            call.push_back(0);
        }
        float tmp2 = S0 * exp ((r - 0.5 * sigma * sigma)* T + sigma * sqrt(T) * z2) - k;
        if (tmp1 > 0 ){
            call.push_back(tmp2* exp(-r * T));
        }
        if (tmp1 < 0){
            call.push_back(0);
        }
    }

    Calculator cal;
    cout << "The call option price is " << cal.average(call, 2*N) << endl;
    
}

int main() {
    
    Q1();
    Q2();
    Q3();
    Q4();
    Q5();
    Q6();
    
    return 0;
};

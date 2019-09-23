//
//  main.cpp
//  HW9
//
//  Created by HUANG runhong on 3/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
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

float sg(int t){
    
    float result = t/30.0 ;
    if (result > 1.0){
        result = 1.0;
    }
    
    return result;
}


float sy(int t){
    std::map <int, float> SY;
    
    int month = t % 12  + 1;
    SY[1] = 0.94;
    SY[2] = 0.76;
    SY[3] = 0.74;
    SY[4] = 0.95;
    SY[5] = 0.98;
    SY[6] = 0.92;
    SY[7] = 0.98;
    SY[8] = 1.10;
    SY[9] = 1.18;
    SY[10] = 1.22;
    SY[11] = 1.23;
    SY[12] = 0.98;
    return SY[month];
}

float bu (float pv_t1, float pv0 ) {
    float result = 0.3 + 0.7 * pv_t1 / pv0;
    return result;
}

float ri (float R, float r){
    float result = 0.28 + 0.14 *  atan ( -8.57 + 430 * (R - r)) ;
    return result;
}

class CIR {
    
public:
    float r0, sigma, keppa, r_;
    
    CIR(float R0, float Sigma, float Keppa, float R_): r0(R0), sigma(Sigma), keppa(Keppa), r_(R_) {
    };
    
    ~CIR(){};
    
public:
    
    float R_MC(float r_start ,float T, float t );
    float R (float r_start, float T , float t );
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
//            cout << r[i+1] << endl;
            if (r[i+1] < 0){

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

float CIR::R (float r_start ,float T = 0, float t = 0) {
    
//    float result = A * exp( -1.0 * B * r_start) ;
    float h = sqrt (keppa * keppa + 2 * sigma * sigma);
    float B = (2.0 * (exp ((T-t)*h ) -1.0) )/ (2 * h + (keppa + h) *(exp ((T-t) *h) -1 ) );
    float A =  (2.0 * h * exp ( (keppa + h ) * (T -t ) / 2.0  )) / (2.0 * h + (keppa + h ) * ( exp ((T- t) * h)  -1));
    A = pow(A, (2.0 * keppa * r_ / (sigma * sigma)) );
    float result = A * exp ( - B * r_start);
    
    return result;
}


class MBS {
    
    float WAC ;
    float r0 ;
    float keppa ;
    float r_ ;
    float sigma ;
    float T ;
    float loan ;
    float oas_x = 0;
    
public:
    
    MBS(float WAC, float R0, float Keppa, float R_, float Sigma, float T , float Loan): WAC(WAC), r0(R0), keppa(Keppa), r_(R_), sigma(Sigma), T(T), loan(Loan) {
    };
    ~MBS(){};
    
    float CPR(float, float , int);
    float MBS_pricing();
    void set_keppa(float k);
    void set_r_(float r);
    void set_sigma(float sigma);
    void set_oas(float oas);
    
    
};
void MBS::set_sigma(float s ){
    this -> sigma = s;
}

void MBS::set_keppa(float k){
    this-> keppa = k;
}

void MBS::set_r_(float r){
    this-> r_ = r;
}

void MBS::set_oas(float x){
    this-> oas_x = x;
}

float MBS::MBS_pricing(){
    
    int nmonth = T * 12 ;
    int npath = 10000;
    float pv[nmonth+1];
    float c[nmonth+1];
    float cpr[nmonth+1];
    float rate[nmonth+1];
    float ip[nmonth+1];
    float tpp[nmonth+1];
    float discount[nmonth+1];
    
    float r = WAC /12.0 ;
    rate[0] = r0;
    pv[0] = loan;
    cpr[0] = 0.0;
    c[0] = 0.0 ;
    ip[0] = 0.0;
    c[0] = 0.0;
    discount[0] = 1.0;
    
    Normal norm(nmonth);
    
    float delta_t = 1.0 /12.0;
    
    float sum = 0;
    
    std::vector<float> values;
    
    for (int j = 0 ; j < npath; j ++){
        
        norm.setseed(j*j);
        norm.generator();
        
        for (int i = 0 ; i < nmonth+1 ; i++){
            rate[i+1] = rate[i] + keppa * ( r_ - rate[i]) * delta_t + sigma * sqrt(delta_t * rate[i] ) * norm.b[i];
                    if (rate[i+1] < 0){
                        rate[i+1] = 0;
                    }
            cpr[i+1] = CPR(rate[i], pv[i], i);
            c[i+1] = pv[i] * r / ( 1.0 - pow( (1+r), -nmonth + i)) + (pv[i] - pv[i] * r * ( 1.0 / ( 1.0 - pow( (1+r), -nmonth + i)  ) -1.0 )) * (1.0 - pow ((1.0 - cpr[i+1]), 1.0/12.0 )) ;
            
            tpp[i+1] = pv[i] * r  * (1.0 / ( 1.0 - pow( (1+r), -nmonth + i)) -1.0 ) + (pv[i] - pv[i] * r * ( 1.0 / ( 1.0 - pow( (1+r), -nmonth + i)  ) -1.0 )) * (1.0 - pow ((1.0 - cpr[i+1]), 1.0/12.0 )) ;
            ip[i+1] = pv[i] * r ;
            pv[i+1] = pv[i] - (tpp[i+1]);
//            cout << i+1<< " " << c[i] << endl;
        }
        
        float value = 0.0;

        for (int i = 0 ; i < nmonth +1 ; i++){
//            discount[i+1] = discount[i] * (exp(-rate[i+1] * delta_t));
            discount[i+1] = discount[i] * (exp((-(rate[i+1] + oas_x))  * delta_t));
            value = value + discount[i] * c[i];
        }
        
        values.push_back(value);
        sum = sum + value;
    }
    
    Calculator cal;
//    cout << cal.average(values, npath) << endl;
    
    return cal.average(values, npath) ;
}


float MBS::CPR(float rate, float pv ,int time ){
    
    CIR interest(r0, sigma, keppa, r_);
    
    float r_10 = -1.0 / 10 * log(interest.R(rate,10,0));

    float RI = ri( WAC , r_10 );
    float BU = bu(pv, loan );
    float SG = sg(time);
    float SY = sy(time);
    
    float cpr = RI * BU * SG * SY ;

    return cpr;
}




int main () {
    
    float WAC = 0.08;
    float r0 = 0.078;
    float keppa = 0.6;
    float r_ = 0.08;
    float sigma = 0.12 ;
    float T = 30;
    float loan = 100000;
    
    MBS Q1(WAC, r0, keppa, r_, sigma, T , loan);
//  Q1
    
    cout << Q1.MBS_pricing() << endl;
    
//    Qb
    for ( int i  = 3 ; i < 10; i++){
        float ke = i /10.0;
        cout << ke << endl;
        Q1.set_keppa(ke);
        cout << Q1.MBS_pricing() << endl;
    }
    
///    Qc
    for ( int i  = 3 ; i < 10; i++){
        float rate = i /100.0;
        cout << rate << endl;
        Q1.set_r_(rate);
        cout << Q1.MBS_pricing()<< endl;
    }
    
///    Qd
    for ( int i  = 10 ; i < 21; i++){
        float sigma = i /100.0;
        cout << sigma << endl;
        Q1.set_sigma(sigma);
        cout << Q1.MBS_pricing() << endl;
    }
//
    
//    Q2
    Q1.MBS_pricing();

    for (int i = 0 ; i < 100 ; i++){
        cout << -i * 0.0001 - 0.012 << endl;
        Q1.set_oas(-i * 0.0001 - 0.012);
        cout << Q1.MBS_pricing() << endl;
    }
    
    //    the OAS is -0.0128;
    
//    Q3
    float y = 0.0005;
    Q1.set_oas(-0.0125);
    float p0 = Q1.MBS_pricing();
    Q1.set_oas(-0.0125 + y);
    float p_plus = Q1.MBS_pricing();
    Q1.set_oas(-0.0125 - y);
    float p_minus = Q1.MBS_pricing();
    float duration = (p_minus - p_plus) / ( 2 * y * p0);
    cout << duration << endl;

    cout << p_plus << " " << p_minus << " " << p0 << endl;
    float convexity = (p_plus + p_minus - 2.0 * p0) / ( 2 * p0 * y * y);
    cout << convexity << endl;
    
    return 0;
}


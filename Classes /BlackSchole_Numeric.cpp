//
//  BlackSchole_Numeric.cpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include "BlackSchole_Numeric.hpp"
#include <cmath>
#include <vector>


void BlackSchole_Numeric::sets0(float S0){
        this->s0 = S0;
};

void BlackSchole_Numeric::setseed(int s) {
    this-> seed = s;
}

void BlackSchole_Numeric::setlength(int l){
    this->length = l;
};

void BlackSchole_Numeric::stock_generator(){
    Normal norm(length);
    norm.generator();
    stock_price.clear();
    for (int i = 0 ; i < length ; i++){
        stock_price.push_back(s0 * exp((r- sigma*sigma/2) * T + sigma * sqrt(T) * norm.a[i]));
    }
};
    
void BlackSchole_Numeric::generator(){
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
    
void BlackSchole_Numeric::generator_red(){
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

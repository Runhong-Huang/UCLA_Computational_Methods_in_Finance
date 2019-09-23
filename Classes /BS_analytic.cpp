//
//  BS_analytic.cpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include "BS_analytic.hpp"
#include <cmath>
#include <vector>

void BlackSchole_Analytic::sets0(float s){
    this-> s0 = s;
}


float BlackSchole_Analytic::normal_table_generator(float n){
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


float BlackSchole_Analytic::call_option(){
    float call_price;
    //        float d1 = 1/(sigma*sqrt(T)) * (log(s0/k) + (r + 1.0/2.0 * sigma * sigma)*T) ;
    float d1 = 1/(sigma*sqrt(T)) * (log(s0/k) + (r + 1.0/2.0 * sigma * sigma)*T) ;
    float d2 = 1/(sigma*sqrt(T)) * (log(s0/k) + (r - 1.0/2.0 * sigma * sigma)*T);
    
    //        cout << normal_table_generator(d1) << endl;
    //        cout << normal_table_generator(d2) << endl;
    
    call_price = s0 * normal_table_generator(d1) - exp(-r * T ) * k * normal_table_generator(d2);
    
    return call_price;
};






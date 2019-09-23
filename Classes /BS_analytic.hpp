//
//  BS_analytic.hpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#ifndef BS_analytic_hpp
#define BS_analytic_hpp

#include <stdio.h>
#include <cmath>
#include <vector>
#include "Normal.hpp"


#endif /* BS_analytic_hpp */

class BlackSchole_Analytic: public Normal{
    
private:
    int length = 1000000;
    float s0,r, sigma,k , T;
    
public:
    BlackSchole_Analytic(float S0, float R , float Sigma, float K, float T): s0(S0),r(R),sigma(Sigma), k(K), T(T) {
    };
    
    ~BlackSchole_Analytic(){};
    
    void sets0(float s);
    
    float normal_table_generator(float n);
   
    float call_option();
};

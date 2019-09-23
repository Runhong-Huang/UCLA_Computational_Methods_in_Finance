//
//  BlackSchole_Numeric.hpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#ifndef BlackSchole_Numeric_hpp
#define BlackSchole_Numeric_hpp

#include <stdio.h>

#endif /* BlackSchole_Numeric_hpp */

#include <cmath>
#include <vector>
#include "Normal.hpp"

class BlackSchole_Numeric:Normal{
    
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
    
    void sets0(float S0);
    void setseed(int s);
    void setlength(int l);
    void stock_generator();
    void generator();
    void generator_red();
    
};


//
//  Normal.hpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#ifndef Normal_hpp
#define Normal_hpp

#include <stdio.h>
#include <cmath>
#include <vector>

#endif /* Normal_hpp */

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
    
    void setseed(long s);
    
    //    set the mean and varaince of A
    void setA(float miu, float sigma);
    
    //    set the mean and varaince of B
    void setB(float miu, float sigma);
    
    void setCorr(float correlation);
    
    void generator();
};

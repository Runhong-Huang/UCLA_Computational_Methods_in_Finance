//
//  Normal.cpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include "Normal.hpp"
#include <cmath>
#include <vector>

    

void Normal::setseed(long s){
        this->seed = s;
    };
    
    //    set the mean and varaince of A
void Normal::setA(float miu, float sigma){
        this->Amiu = miu;
        this->Asigma = sigma;
    };
    
    //    set the mean and varaince of B
void Normal::setB(float miu, float sigma){
        this->Bmiu = miu;
        this->Bsigma = sigma;
    }
    
void Normal::setCorr(float correlation){
        this->rou = correlation;
    }
    
void Normal::generator(){
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
}

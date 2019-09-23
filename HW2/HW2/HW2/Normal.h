//
//  Header.h
//  HW2
//
//  Created by HUANG runhong on 1/21/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#ifndef Header_h
#define Header_h


#endif /* Header_h */

class Calculator {
    
private:
    
    
public:
    
    Calculator(){
    };
    
    ~Calculator(){};
    
    
};



// The class used to generate normal random variables
//The class takes three input: seeeds and length and correlation
//The class have two vectors a and b. each contains random variables
class Normal{
private:
    int length;
    double seed;
    double Amiu, Asigma, Bmiu, Bsigma, rou;
    Calculator c;
    
public:
    Normal(double s, int num):seed(s),length(num) {
        Amiu = 0; Asigma = 1;
        Bmiu = 0; Bsigma = 1;
        rou = 0;
    }
    ~Normal(){};
    
    // define two vectors with varaible length.
    std::vector<float> a;
    std::vector<float> b;
    
    //    set the mean and varaince of A
    void setA(float miu, float sigma){
        this->Amiu = miu;
        this->Asigma = sigma;
    }
    
    //    set the mean and varaince of B
    void setB(float miu, float sigma){
        this->Bmiu = miu;
        this->Bsigma = sigma;
    }
    
    void setCorr(float correlation){
        this->rou = correlation;
    }
    
    void generator(){
        srand(seed);
        for (int i =0 ; i < length;){
            float v1 = rand()*1.0/(RAND_MAX*1.0) * 2 -1  ;
            float v2 = rand()*1.0/(RAND_MAX*1.0) * 2 -1 ;
            float w = v1* v1 + v2* v2;
            if (w <= 1 ) {
                float z1 = v1 * sqrt(-2 * log (w) /w);
                float z2 = v2 * sqrt(-2 * log (w) /w);
                a.push_back(Amiu + Asigma * z1);
                b.push_back(Bmiu + Bsigma * (rou * z1 + sqrt(1-rou*rou)) * z2);
                i = i +1;
            }
        }
    }
};

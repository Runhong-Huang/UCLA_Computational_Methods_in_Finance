//
//  Calculator.hpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#ifndef Calculator_hpp
#define Calculator_hpp

#include <iostream>
#include <cmath>
#include <vector>


using namespace std;

class Calculator {
    
private:
    
public:
    template <class T>
    double average(std::vector<T> arr, int length);
    template <class T>
    double variance(std::vector<T> arr, int length);
    template <class T,class U>
    double corr(std::vector<T> arr1,std::vector<U> arr2, int length);
};

#endif /* Calculator_hpp */

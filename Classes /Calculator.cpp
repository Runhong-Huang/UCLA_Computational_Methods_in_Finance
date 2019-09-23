//
//  Calculator.cpp
//  HW4
//
//  Created by HUANG runhong on 2/2/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//


#include "Calculator.hpp"
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;


template <class T>
double Calculator::average(std::vector<T> arr, int length){
    float sum =0;
    for (int i = 0 ; i  < length ; i++)
        sum = sum + arr[i];
    //    cout << sum << endl;
    return sum/length;
};

template <class T>
double Calculator::variance(std::vector<T> arr, int length) {
    //    calculate the variance of array with length of "length"
    double average = this->average(arr, length);
    double var = 0 ;
    for (int i = 0 ; i  < length ; i++)
        var = var + (arr[i] - average) * (arr[i] - average);
    return var/(length-1);
};


template <class T,class U>
double Calculator::corr(std::vector<T> arr1,std::vector<U> arr2, int length) {
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

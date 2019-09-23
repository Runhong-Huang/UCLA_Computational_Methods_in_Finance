//
//  main.cpp
//  Final_1
//
//  Created by HUANG runhong on 3/20/19.
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



int s (float n ){

    float sum = 0;
    int count = 0 ;
    while (sum < n){
        sum = sum + rand()*1.0/(RAND_MAX*1.0) ;
//        cout << sum << endl;
        count ++;
    }
    return count;
}


int main (){
    
//    cout << s(1.1);
    
    std::vector<float> sum;
    float tmp;
    
    for (int i = 0 ; i < 10000; i++){
        srand(i*i);
        tmp = max((float)4.54 - s(1.1), float (0));
//        cout << s(1.1) << endl;
        sum.push_back(tmp);
    }
    
    Calculator cal;
    cout << cal.average(sum, 10000) << endl;
    
    return 0;
}

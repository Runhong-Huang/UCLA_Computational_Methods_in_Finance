//
//  main.cpp
//  HW1
//
//  Created by HUANG runhong on 1/10/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <ctime>
#include <fstream>
using namespace std;

//create a function to generate the random varaiables
double  random_generator ( double arr[] , int n = 100 , int m  = pow(2,31) - 1 , int a = pow(7,5)) {
    //    n is the numer of the random varaible
    // x is the number from 0 to m, v is uniform from [0,1]
    int seed = 19900523;
    long x[n];
    x[0] = seed;
    for (int i = 0 ; i < n ; i++ ){
        x[i+1] = (a * x[i]) % m ;
        arr[i] = ((float) x[i] / m);
        //        cout << arr[i] << endl;
    }
    return *arr;
}

template <class T>
double takeAverage(T arr[], int length) {
    float sum =0;
    for (int i = 0 ; i  < length ; i++)
        sum = sum + arr[i];
//    cout << sum << endl;
    return sum/length;
}

template <class T>
//function to return variance
double takeVariance(T arr[], int length) {
    //    calculate the variance of array with length of "length"
    T average = takeAverage(arr, length);
    double var = 0 ;
    for (int i = 0 ; i  < length ; i++)
        var = var + (arr[i] - average) * (arr[i] - average);
    return var/length;
}


int Q2_random_generator (int arr[], int n) {
    //    generate rv for Q2. -1 wp 0.3 , 0 wp 0.35, 1 wp 0.2, 2 wp 0.15
    srand((int)time(0));
    for (int i = 0 ; i < n ; i++ ){
        float tmp = rand()*1.0/(RAND_MAX*1.0);
//        cout << tmp << endl;
        if (tmp < 0.3)
            arr[i] = -1;
        else if (tmp < 0.65)
            arr[i] = -0;
        else if (tmp < 0.85)
            arr[i] = 1;
        else
            arr[i] = 2;
    }
    return 0;
}


double uniform (double arr[], int n){
    for (int i=0; i< n ; i++)
    {
    arr[i] = rand()*1.0/(RAND_MAX*1.0);
    }
    return 0;
}

int Binomial(int arr[],int length, int n , float p ) {
    srand((int)time(0));
    float tmp = 0;
    for (int i = 0; i < length ; i++){
        for (int j =0; j < n ; j++){
            tmp = rand()*1.0/(RAND_MAX*1.0);
            if (tmp < p)
                arr[i]++;
        }
//        cout << arr[i] << endl;
    }
    return 0;
}


float exponetial(float arr[], float lamda, int length) {
    srand((int)time(0));
    for (int i=0; i< length ; i++)
    {
        arr[i] = -1/lamda * log ( rand()*1.0/(RAND_MAX*1.0));
    }
    return 0;
}

double normal_BoxMuller (double arr[],int length){
    for (int i = 0 ; i < length;i = i +2){
        float u1 = rand()*1.0/(RAND_MAX*1.0);
        float u2 = rand()*1.0/(RAND_MAX*1.0);
        arr[i] = sqrt((-2 * log (u1))) * cos( 2 * acos(-1) * u2);
        arr[i+1] = sqrt((-2 * log (u1))) * sin( 2 * acos(-1) * u2);


//        cout << arr[2*i] << endl;
    }
    return 0;
}

double normal_pm (double arr[], int length){
    for (int i = 0 ; i < length; ){
        float v1 = rand()*1.0/(RAND_MAX*1.0) * 2 -1  ;
        float v2 = rand()*1.0/(RAND_MAX*1.0) * 2 -1 ;
        float w = v1* v1 + v2* v2;
        if (w <= 1 ) {
            arr[i] = v1 * sqrt(-2 * log (w) /w) ;
            arr[i+1] = v2 * sqrt(-2 * log (w) /w) ;
            i = i +2;
        }
    }
    return 0;
}

int main() {
    
//    create array with length of 10000
    double Q1_random[10000] ={};
    random_generator(Q1_random, 10000);

    int Q1_length =  sizeof(Q1_random)/sizeof(*Q1_random);
    
//    calculate the average value of Q1
    float Q1_average = takeAverage(Q1_random, Q1_length);
    cout << "The average of the Q1  is :" << endl;
    cout << Q1_average << endl;
    
//    calculate the variance of the array
    cout << "The std of the Q1  is :" << endl;
    cout << sqrt(takeVariance(Q1_random, Q1_length)) << endl;

    
//    Q1b repeat the Q1a using random function
    double Q1b_random [10000] = {};
    uniform (Q1b_random , 10000);
    
    cout << "The average of the Q1(b)  is :" << endl;
    cout << takeAverage(Q1b_random, Q1_length) << endl;
    
    cout << "The std of the Q1(b)  is :" << endl;
    cout << sqrt(takeVariance(Q1b_random, Q1_length)) << endl;
    
// Q2 question 2
    
    int Q2[10000] = {};  // define a empty array
    Q2_random_generator(Q2, 10000);
    cout << "Q2: the average is " <<  takeAverage(Q2, 10000) << endl;
    cout << "Q2: the std is " <<  sqrt(takeVariance(Q2, 10000)) << endl;
    
    
//    output the data to Q2.csv
    ofstream myfile;
    myfile.open ("Q2.csv");
    for (int i =0 ; i < 10000 ; i ++){
        myfile << Q2[i]  << endl;
    }
    myfile.close();
    
    
//  Q3 Question 3

    //    Q3 a
    int Q3[1000] = {};
    Binomial(Q3, 1000,44, 0.64);
    cout << "Q3: the average is " <<  takeAverage(Q3, 1000) << endl;
    
    myfile.open ("Q3.csv");
    for (int i =0 ; i < 1000 ; i ++){
        myfile << Q3[i]  << endl;
    }
    myfile.close();
    
    //   Q3 b
    int Q3b[1000000] = {};
    Binomial(Q3b, 1000000,44, 0.64);
    int Q3_count = 0;
    for (long long int j = 0; j < 1000000 ; j++){
        if (Q3b[j] > 39)
                Q3_count++;
    }
    
    cout << "Q3 from 100000 test, the number of times of x>= 40 is " <<Q3_count << endl;
    
//    Q4 exponential distribution
    float Q4[10000] = {};
    
    exponetial(Q4, 1.5, 10000);
    
    cout << "Q4 average" << takeAverage(Q4,10000) << endl;
    cout << "Q4 standard deviation" << sqrt(takeVariance(Q4,10000)) << endl;
    
    int Q4_count1 = 0; int Q4_count2 = 0;
    for ( int i = 0 ; i < 10000 ; i++){
        if (Q4[i] >= 1){
            Q4_count1++;
            if (Q4[i] >= 4){
                Q4_count2++;
            }
        }
    }
    
    cout << "Q4: (2) P(X>1) is " << Q4_count1/10000.0  << endl;
    cout << "Q4: (2) P(X>4) is " << Q4_count2/10000.0  << endl;
    
    myfile.open ("Q4.csv");
    for (int i =0 ; i < 10000 ; i ++){
        myfile << Q4[i]  << endl;
    }
    myfile.close();
    
//    Q5 normal distribution
    
  
    //    Q5 (b) random numer generater from box muller methods
    double Q5[5000] = {};
    
    std::clock_t start;
    double duration1, duration2;
    
//   track the time used for pm method
    start = std::clock();
    normal_BoxMuller (Q5, 5000);
    duration1 = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    
    //    Q5 (b)
    cout << "Q5 (c): The average value of normal from box muller is " << takeAverage(Q5, 5000) << endl;
    cout << "Q5 (c): The variance value of normal from box muller is " <<takeVariance(Q5, 5000) << endl;

    //    Q5 (c) random numer generater from polar marsaglia
    double Q5c[5000] = {};
    
    
    start = std::clock();
    normal_pm (Q5c, 5000);
    duration2 = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
   
    
    //    Q5(d)
    cout << "Q5 (d): The average value of normal from PM is " << takeAverage(Q5c, 5000) << endl;
    cout << "Q5 (d): The variance value of normal from PM is " <<takeVariance(Q5c, 5000) << endl;
    
  //    Q5(f)
    cout << "The time used for bm method is " << duration1 << endl;
    cout << "The time used for pm method is " << duration2 << endl;

    return 0;
}



//
//  main.cpp
//  HW7
//
//  Created by HUANG runhong on 2/25/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include <iostream>
#include <cmath>
using namespace std;

float *regression(float **matrix_A , float *vec_B , int n ){
    //  find the parameter in x for the Ax = b, where A is k by k matrix, b is k by 1 matrix
    
    float *a = new float [n];
    
    for (int i = 0 ; i < n ; i++){
        a[i] =0;
    }
    
    //    convert
    //    perform Gaussian eliminiation on A and B
    for (int i = 0 ; i < n ; i ++){
        for (int j = i+1 ; j < n ; j++){
            float tmp = matrix_A[j][i] / matrix_A[i][i];
            vec_B[j] = vec_B[j] - tmp * vec_B[i];
            for (int k = 0 ; k < n ; k++){
                matrix_A[j][k] = matrix_A[j][k] - tmp * matrix_A[i][k];
            }
        }
    }
    
    ////    check the determinate the gaussian eliminated matrix a
    try{
        
        if (matrix_A[n-1][n-1] == 0){
            throw 1;
        }
        
        //    calculating A
        for (int i = n-1; i >= 0 ; i--){
            float tmp = 0;
            for (int j = i+1 ; j < n; j++ ){
                tmp = tmp + matrix_A[i][j] * a[j];
            }
            a[i] = (vec_B[i] - tmp) / matrix_A[i][i];
        }
    }
    
    catch (int x ){
        cout << x << "Error: The matrix A is singular and not invertible " << endl;
    }
    
    return (a);
}


class BS_PDE{

private:
    float s0 = log(10.0) ;
    float sigma = 0.20 ;
    float r = 0.04 ;
    float k =  10.0;
    float t = 0.5;
    float s0_min = log(4.0);
    float s0_max = log(16.0);

    float delta_t = 0.002;
    float delta_x = sigma * sqrt( delta_t) ;

public:

    float *explicit_FD();
    float *implicit_FD();
    float *cn_FD();
    float *am_explicit_call();
    float *am_explicit_put();
    float *am_implicit_call();
    float *am_implicit_put();
    float *am_cn_call();
    float *am_cn_put();

};

float *BS_PDE::explicit_FD(){

    int ntime =  t/ delta_t ;
    int nprice = ceil((s0_max - s0_min) / (delta_x));
    
    //    cout << "delta t :" << delta_t <<  endl;
    //    cout << "delta x :" << delta_x <<  endl;
    //    cout << s0_max << " " << s0_min << endl;
    cout << nprice << endl;
    
    float **option = new float *[nprice+1];
    float *stock = new float[nprice+1];
    
    for (int i = 0 ; i <= nprice +1 ; i++) {
        option[i] = new float [ntime+1];
    }
    
    for (int i = 0 ; i < nprice + 1 ; i++){
        stock[i] = log(10) + 52 * delta_x - i * delta_x;
        cout << i << " " << stock[i] << endl;
    }
    
    for (int i = 0 ; i < nprice +1; i++){
        option[i][ntime+1] = max ((float)(k -exp(stock[i])), (float) 0.0);
        option[i][ntime+1] = exp(-r * t) * option[i][ntime+1];
        //        cout << option[i][ntime] <<endl;
    }
    
    float pu = delta_t * ( (sigma* sigma)/(2.0 * delta_x * delta_x) + (r - 0.5 * sigma * sigma)/(2.0 * delta_x));
    float pm = 1.0 - delta_t * sigma * sigma / ( delta_x * delta_x) - r * delta_t;
    float pd = delta_t * ((sigma * sigma) / (2.0 * (delta_x * delta_x)) - (r - 0.5 * sigma * sigma)/(2.0 * delta_x));
    
    cout << "Pu : " << pu << endl;
    cout << "Pm : " << pm << endl;
    cout << "Pd : " << pd << endl;
    
    for (int j = ntime ; j >= 0 ; j--){
        for (int i = 1 ; i < nprice ; i++){
            option[i][j] = pu * option[i+1][j+1] + pm * option[i][j+1] + pd * option[i-1][j+1];
        }
        option[0][j] = option[1][j] ;
        option[nprice][j] = option[nprice-1][j] - ( exp(stock[nprice]) - exp(stock[nprice-1]));
    }
    
    for (int i = 0; i < nprice+1 ; i++){
        //        cout << option[i][0] << endl;
    }
    //
    //    for (int i =0 ; i < nprice ; i ++){
    //        cout << option[i][0] << " ,";
    //    }
    
    cout << option[50][0] << endl;
    
    return 0;
}



float *BS_PDE::implicit_FD(){

    int ntime =  t/ delta_t ;
    //int nprice = (s0_max - s0_min) / (delta_x);
    int nprice = 2*52+1;
//    delta_x =(s0_max - s0_min) / nprice;


    float *option = new float [nprice];
    float *stock = new float[nprice];

    for (int i = 0 ; i < nprice  ; i++){
        stock[i] = log(10) + 52 * delta_x  - i * delta_x;
    }

    for (int i = 0 ; i < nprice ; i++){
        option[i] = max ((float)(k - exp(stock[i])), (float) 0.0);
//        option[i][ntime] = exp(-r * t) * option[i][ntime] ;
    }
    
    float v = r - (0.5 * sigma * sigma);
    float pu = - 0.5 * delta_t * ( (sigma * sigma)/ (delta_x * delta_x)  + (v)/delta_x );
    float pm = 1 + delta_t * sigma *sigma /(delta_x * delta_x ) + r * delta_t;
    float pd = - 0.5 * delta_t * ( (sigma * sigma)/ (delta_x * delta_x)  - (v)/delta_x );

//    cout << delta_t *v/delta_x << endl;
    cout << "Pu : " << pu << endl;
    cout << "Pm : " << pm << endl;
    cout << "Pd : " << pd << endl;

    
//    define matrix A and B
    float **A = new float *[nprice];
    
    for (int i = 0 ; i < nprice +1 ; i++) {
//        cout << i << endl;
        A[i] = new float [nprice];
    }
    for(int i = 0; i < nprice;i++){
        for(int j = 0; j < nprice; j++){
            A[i][j]=0;
        }
    }
 
//    define B vector
    for (int j = ntime-1 ; j >= 0 ; j--){
        
        A[0][0] = 1.0;
        A[0][1] = -1.0;
        A[nprice-1][nprice-2] = 1.0;
        A[nprice-1][nprice-1] = -1.0;
        
        for (int i = 1 ; i < nprice ; i++){
            A[i][i-1] = pu;
            A[i][i] = pm;
            A[i][i+1] = pd;
        }
        option[0] = 0;
        option[nprice-1] = exp(stock[nprice-2]) - exp(stock[nprice-1])  ;
        option = regression(A, option, nprice);
//        cout<<option[52]<<endl;
        
//        cout << option[0][j] << " " << option[1][j] << endl;
    }
    
    
    cout << option[52] << endl;
    
    return 0;

}


float *BS_PDE::cn_FD(){
    
    int ntime =  t/ delta_t ;
    //int nprice = (s0_max - s0_min) / (delta_x);
    int nprice = 2*52+1;
    //    delta_x =(s0_max - s0_min) / nprice;
    
    
    float *option = new float [nprice];
    float *stock = new float[nprice];
    
    
    for (int i = 0 ; i < nprice  ; i++){
        stock[i] = log(10) + 52 * delta_x  - i * delta_x;
    }
    
    for (int i = 0 ; i < nprice ; i++){
        option[i] = max ((float)(k - exp(stock[i])), (float) 0.0);
    }
    
    float v = r - (0.5 * sigma * sigma);
    float pu = - 0.25 * delta_t * ((sigma* sigma)/(delta_x * delta_x) + (v)/(delta_x) );
    float pm = 1.0 + delta_t * sigma *sigma /(2.0 * delta_x * delta_x ) + r * delta_t / 2.0;
    float pd = -0.25 * delta_t *  ((sigma * sigma)/(delta_x * delta_x) - (v)/(delta_x));
    
    //    cout << delta_t *v/delta_x << endl;
    cout << "Pu : " << pu << endl;
    cout << "Pm : " << pm << endl;
    cout << "Pd : " << pd << endl;
    
    //    define matrix A and B
    float **A = new float *[nprice];
    
    for (int i = 0 ; i < nprice +1 ; i++) {
        //        cout << i << endl;
        A[i] = new float [nprice];
    }
    
// define Z
    float *Z = new float [nprice];
    
    //    define B vector
    for (int j = ntime-1 ; j >= 0 ; j--){
        
        for (int i = 1 ; i < nprice - 1 ; i++){
            Z[i] = - pu * option[i-1] - ( pm-2.0) * option[i] - pd * option[i+1];
        }
        
        Z[nprice-1] = 0;
        Z[0] = exp(stock[nprice-2] ) - exp(stock[nprice-1]);

        for(int i = 0; i < nprice;i++){
            for(int j = 0; j < nprice; j++){
                A[i][j]=0;
            }
        }
        
        A[0][0] = 1.0;
        A[0][1] = -1.0;
        A[nprice-1][nprice-2] = 1.0;
        A[nprice-1][nprice-1] = -1.0;
        
        for (int i = 1 ; i < nprice ; i++){
            A[i][i-1] = pu;
            A[i][i] = pm;
            A[i][i+1] = pd;
        }
        
        option = regression(A, Z, nprice-1);
        
    }
    
    
    cout << option[52] << endl;

    
    return 0;
}



float *BS_PDE::am_explicit_put(){
    
    int ntime =  t/ delta_t ;
    int nprice = (s0_max - s0_min) / (delta_x);
    
    //    cout << "delta t :" << delta_t <<  endl;
    //    cout << "delta x :" << delta_x <<  endl;
    //    cout << s0_max << " " << s0_min << endl;
    cout << nprice << endl;
    
    float **option = new float *[nprice+1];
    float *stock = new float[nprice+1];
    
    for (int i = 0 ; i <= nprice +1 ; i++) {
        option[i] = new float [ntime+1];
    }
    
    for (int i = 0 ; i < nprice + 1 ; i++){
        stock[i] = log(10) + 52 * delta_x - i * delta_x;
        
        //        cout << i << " " << stock[i] << endl;
    }
    
    for (int i = 0 ; i < nprice +1; i++){
        option[i][ntime+1] = max ((float)(k -exp(stock[i])), (float) 0.0);
        option[i][ntime+1] = exp(-r * t) * option[i][ntime+1];
        //        cout << option[i][ntime] <<endl;
    }
    
    float pu = delta_t * ( (sigma* sigma)/(2.0 * delta_x * delta_x) + (r - 0.5 * sigma * sigma)/(2.0 * delta_x));
    float pm = 1.0 - delta_t * sigma * sigma / ( delta_x * delta_x) - r * delta_t;
    float pd = delta_t * ((sigma * sigma) / (2.0 * (delta_x * delta_x)) - (r - 0.5 * sigma * sigma)/(2.0 * delta_x));
    
    //    cout << "Pu : " << pu << endl;
    //    cout << "Pm : " << pm << endl;
    //    cout << "Pd : " << pd << endl;
    
    for (int j = ntime ; j >= 0 ; j--){
        for (int i = 1 ; i < nprice ; i++){
            option[i][j] = pu * option[i+1][j+1] + pm * option[i][j+1] + pd * option[i-1][j+1];
        }
        option[0][j] = option[1][j] ;
        option[nprice][j] = option[nprice-1][j] - ( exp(stock[nprice]) - exp(stock[nprice-1]));
        
        for (int i = 0 ; i <= nprice; i++){
            float tmp = k - exp(stock[i]);
            if (option[i][j] < tmp){
                option[i][j] = tmp;
            }
        }
    }
    
    for (int i =nprice ; i >0 ; i--){
        cout << option[i][0] << "," ;
    }
    return 0;
}



float *BS_PDE::am_explicit_call(){
    
    int ntime =  t/ delta_t ;
    int nprice = (s0_max - s0_min) / (delta_x);
    
    //    cout << "delta t :" << delta_t <<  endl;
    //    cout << "delta x :" << delta_x <<  endl;
    //    cout << s0_max << " " << s0_min << endl;
    cout << nprice << endl;
    
    float **option = new float *[nprice+1];
    float *stock = new float[nprice+1];
    
    for (int i = 0 ; i <= nprice +1 ; i++) {
        option[i] = new float [ntime+1];
    }
    
    for (int i = 0 ; i < nprice + 1 ; i++){
        stock[i] = log(10) + 52 * delta_x - i * delta_x;
        
        //        cout << i << " " << stock[i] << endl;
    }
    
    for (int i = 0 ; i < nprice +1; i++){
        option[i][ntime+1] = max ((float)(exp(stock[i])) - k , (float) 0.0);
        option[i][ntime+1] = exp(-r * t) * option[i][ntime+1];
        //        cout << option[i][ntime] <<endl;
    }
    
    float pu = delta_t * ( (sigma* sigma)/(2.0 * delta_x * delta_x) + (r - 0.5 * sigma * sigma)/(2.0 * delta_x));
    float pm = 1.0 - delta_t * sigma * sigma / ( delta_x * delta_x) - r * delta_t;
    float pd = delta_t * ((sigma * sigma) / (2.0 * (delta_x * delta_x)) - (r - 0.5 * sigma * sigma)/(2.0 * delta_x));
    
    //    cout << "Pu : " << pu << endl;
    //    cout << "Pm : " << pm << endl;
    //    cout << "Pd : " << pd << endl;
    
    for (int j = ntime ; j >= 0 ; j--){
        for (int i = 1 ; i < nprice ; i++){
            option[i][j] = pu * option[i+1][j+1] + pm * option[i][j+1] + pd * option[i-1][j+1];
        }
        option[0][j] = option[1][j] ;
        option[nprice-1][j] = option[nprice-1][j] - ( exp(stock[nprice]) - exp(stock[nprice-1]));
        
        for (int i = 0 ; i <= nprice; i++){
            float tmp = exp(stock[i]) - k;
            if (option[i][j] < tmp){
                option[i][j] = tmp;
            }
        }
    }
    
    
    for (int i =0 ; i < nprice ; i ++){
        cout << option[i][0] << "," ;
    }
    return 0;
}



    

float *BS_PDE::am_cn_call(){
    
    int ntime =  t/ delta_t ;
    //int nprice = (s0_max - s0_min) / (delta_x);
    int nprice = 2*52+1;
    //    delta_x =(s0_max - s0_min) / nprice;
    
    float *option = new float [nprice];
    float *stock = new float[nprice];
    
    
    for (int i = 0 ; i < nprice  ; i++){
        stock[i] = log(10) + 52 * delta_x  - i * delta_x;
    }
    
    for (int i = 0 ; i < nprice ; i++){
        option[i] = max ((float)(exp(stock[i]) - k ), (float) 0.0);
    }
    
    float v = r - (0.5 * sigma * sigma);
    float pu = - 0.25 * delta_t * ((sigma* sigma)/(delta_x * delta_x) + (v)/(delta_x) );
    float pm = 1.0 + delta_t * sigma *sigma /(2.0 * delta_x * delta_x ) + r * delta_t / 2.0;
    float pd = -0.25 * delta_t *  ((sigma * sigma)/(delta_x * delta_x) - (v)/(delta_x));
    
    //    cout << delta_t *v/delta_x << endl;
    cout << "Pu : " << pu << endl;
    cout << "Pm : " << pm << endl;
    cout << "Pd : " << pd << endl;
    
    //    define matrix A and B
    float **A = new float *[nprice];
    
    for (int i = 0 ; i < nprice +1 ; i++) {
        //        cout << i << endl;
        A[i] = new float [nprice];
    }
    
    // define Z
    float *Z = new float [nprice];
    //    define B vector
    for (int j = ntime-1 ; j >= 0 ; j--){
        
        for (int i = 1 ; i < nprice - 1 ; i++){
            Z[i] = - pu * option[i-1] - ( pm-2.0) * option[i] - pd * option[i+1];
        }
        
        Z[0] = 0;
        Z[nprice-1] = exp(stock[nprice-2] ) - exp(stock[nprice-1]);
        
        for(int i = 0; i < nprice;i++){
            for(int j = 0; j < nprice; j++){
                A[i][j]=0;
            }
        }
        
        A[0][0] = 1.0;
        A[0][1] = -1.0;
        A[nprice-1][nprice-2] = 1.0;
        A[nprice-1][nprice-1] = -1.0;
        
        for (int i = 1 ; i < nprice ; i++){
            A[i][i-1] = pu;
            A[i][i] = pm;
            A[i][i+1] = pd;
        }
        
        option = regression(A, Z, nprice-1);
        
        for (int i = 0 ; i <= nprice; i++){
            float tmp = exp(stock[i] - k );
            if (option[i] < tmp){
                option[i] = tmp;
            }
        }
        
    }
    
    for (int i =0 ; i < nprice ; i ++){
        cout << option[i] << "," ;
    }
    
    return 0;
}

float *BS_PDE::am_cn_put(){
    
    int ntime =  t/ delta_t ;
    //int nprice = (s0_max - s0_min) / (delta_x);
    int nprice = 2*52+1;
    //    delta_x =(s0_max - s0_min) / nprice;
    
    float *option = new float [nprice];
    float *stock = new float[nprice];
    
    
    for (int i = 0 ; i < nprice  ; i++){
        stock[i] = log(10) + 52 * delta_x  - i * delta_x;
    }
    
    for (int i = 0 ; i < nprice ; i++){
        option[i] = max ((float)(k - exp(stock[i])), (float) 0.0);
    }
    
    float v = r - (0.5 * sigma * sigma);
    float pu = - 0.25 * delta_t * ((sigma* sigma)/(delta_x * delta_x) + (v)/(delta_x) );
    float pm = 1.0 + delta_t * sigma *sigma /(2.0 * delta_x * delta_x ) + r * delta_t / 2.0;
    float pd = -0.25 * delta_t *  ((sigma * sigma)/(delta_x * delta_x) - (v)/(delta_x));
    
    //    cout << delta_t *v/delta_x << endl;
    cout << "Pu : " << pu << endl;
    cout << "Pm : " << pm << endl;
    cout << "Pd : " << pd << endl;
    
    //    define matrix A and B
    float **A = new float *[nprice];
    
    for (int i = 0 ; i < nprice +1 ; i++) {
        //        cout << i << endl;
        A[i] = new float [nprice];
    }
    
    // define Z
    float *Z = new float [nprice];
    //    define B vector
    for (int j = ntime-1 ; j >= 0 ; j--){
        
        for (int i = 1 ; i < nprice - 1 ; i++){
            Z[i] = - pu * option[i-1] - ( pm-2.0) * option[i] - pd * option[i+1];
        }
        
        Z[0] = 0;
        Z[nprice -1] = exp(stock[nprice-2] ) - exp(stock[nprice-1]);
        
        for(int i = 0; i < nprice;i++){
            for(int j = 0; j < nprice; j++){
                A[i][j]=0;
            }
        }
        
        A[0][0] = 1.0;
        A[0][1] = -1.0;
        A[nprice-1][nprice-2] = 1.0;
        A[nprice-1][nprice-1] = -1.0;
        
        for (int i = 1 ; i < nprice ; i++){
            A[i][i-1] = pu;
            A[i][i] = pm;
            A[i][i+1] = pd;
        }
        
        option = regression(A, Z, nprice-1);
        
        for (int i = 0 ; i <= nprice; i++){
            float tmp = k - exp(stock[i]);
            if (option[i] < tmp){
                option[i] = tmp;
            }
        }
        
        
        
    }
    
    for (int i =0 ; i < nprice ; i ++){
        cout << option[i] << "," ;
    }
    
    return 0;
}


float *BS_PDE::am_implicit_put(){
    
    int ntime =  t/ delta_t ;
    //int nprice = (s0_max - s0_min) / (delta_x);
    int nprice = 2*52+1;
    //    delta_x =(s0_max - s0_min) / nprice;
    
    
    float *option = new float [nprice];
    float *stock = new float[nprice];
    
    
    for (int i = 0 ; i < nprice  ; i++){
        stock[i] = log(10) + 52 * delta_x  - i * delta_x;
    }
    
    for (int i = 0 ; i < nprice ; i++){
        option[i] = max ((float)(k - exp(stock[i])), (float) 0.0);
        //        option[i][ntime] = exp(-r * t) * option[i][ntime] ;
    }
    
    float v = r - (0.5 * sigma * sigma);
    float pu = - 0.5 * delta_t * ( (sigma * sigma)/ (delta_x * delta_x)  + (v)/delta_x );
    float pm = 1 + delta_t * sigma *sigma /(delta_x * delta_x ) + r * delta_t;
    float pd = - 0.5 * delta_t * ( (sigma * sigma)/ (delta_x * delta_x)  - (v)/delta_x );
    
    //    cout << delta_t *v/delta_x << endl;
    cout << "Pu : " << pu << endl;
    cout << "Pm : " << pm << endl;
    cout << "Pd : " << pd << endl;
    
    
    //    define matrix A and B
    float **A = new float *[nprice];
    
    for (int i = 0 ; i < nprice +1 ; i++) {
        //        cout << i << endl;
        A[i] = new float [nprice];
    }
    for(int i = 0; i < nprice;i++){
        for(int j = 0; j < nprice; j++){
            A[i][j]=0;
        }
    }
    
    //    define B vector
    for (int j = ntime-1 ; j >= 0 ; j--){
        
        A[0][0] = 1.0;
        A[0][1] = -1.0;
        A[nprice-1][nprice-2] = 1.0;
        A[nprice-1][nprice-1] = -1.0;
        
        for (int i = 1 ; i < nprice ; i++){
            A[i][i-1] = pu;
            A[i][i] = pm;
            A[i][i+1] = pd;
        }
        option[0] = 0;
        option[nprice-1] = exp(stock[nprice-2]) - exp(stock[nprice-1])  ;
        option = regression(A, option, nprice);

        
        for (int i = 0 ; i <= nprice; i++){
            float tmp = k - exp(stock[i]);
            if (option[i] < tmp){
                option[i] = tmp;
            }
        }
        
        
    }
    
    for (int i =0 ; i < nprice ; i ++){
        cout << option[i] << "," ;
    }
    
    
    return 0;
    
}


float *BS_PDE::am_implicit_call(){
    
    int ntime =  t/ delta_t ;
    //int nprice = (s0_max - s0_min) / (delta_x);
    int nprice = 2*52+1;
    //    delta_x =(s0_max - s0_min) / nprice;
    
    
    float *option = new float [nprice];
    float *stock = new float[nprice];
    
    
    for (int i = 0 ; i < nprice  ; i++){
        stock[i] = log(10) + 52 * delta_x  - i * delta_x;
    }
    
    for (int i = 0 ; i < nprice ; i++){
        option[i] = max ((float)( exp(stock[i])- k), (float) 0.0);
        //        option[i][ntime] = exp(-r * t) * option[i][ntime] ;
    }
    
    float v = r - (0.5 * sigma * sigma);
    float pu = - 0.5 * delta_t * ( (sigma * sigma)/ (delta_x * delta_x)  + (v)/delta_x );
    float pm = 1 + delta_t * sigma *sigma /(delta_x * delta_x ) + r * delta_t;
    float pd = - 0.5 * delta_t * ( (sigma * sigma)/ (delta_x * delta_x)  - (v)/delta_x );
    
    //    cout << delta_t *v/delta_x << endl;
    cout << "Pu : " << pu << endl;
    cout << "Pm : " << pm << endl;
    cout << "Pd : " << pd << endl;
    
    
    //    define matrix A and B
    float **A = new float *[nprice];
    
    for (int i = 0 ; i < nprice +1 ; i++) {
        //        cout << i << endl;
        A[i] = new float [nprice];
    }
    for(int i = 0; i < nprice;i++){
        for(int j = 0; j < nprice; j++){
            A[i][j]=0;
        }
    }
    
    //    define B vector
    for (int j = ntime-1 ; j >= 0 ; j--){
        
        A[0][0] = 1.0;
        A[0][1] = -1.0;
        A[nprice-1][nprice-2] = 1.0;
        A[nprice-1][nprice-1] = -1.0;
        
        for (int i = 1 ; i < nprice ; i++){
            A[i][i-1] = pu;
            A[i][i] = pm;
            A[i][i+1] = pd;
        }
        option[nprice-1] = 0;
        option[0] = exp(stock[nprice-2]) - exp(stock[nprice-1])  ;
        option = regression(A, option, nprice);
        cout<<option[52]<<endl;
        
        for (int i = 0 ; i <= nprice; i++){
            float tmp = exp(stock[i]) - k;
            if (option[i] < tmp){
                option[i] = tmp;
            }
        }
        
       
        
    }
    
    for (int i =0 ; i < nprice ; i ++){
        cout << option[i] << "," ;
    }
//    cout<<option[52]<<endl;
    
    return 0;
    
}

void Q1(){
    BS_PDE Q1;
    Q1.explicit_FD();
    
    Q1.implicit_FD();
    Q1.cn_FD();
    
}

void Q2(){
    BS_PDE Q2;
    
//    Q2.am_explicit_call();
//    Q2.am_explicit_put();
    
//    Q2.am_implicit_call();
//    Q2.am_implicit_put();
//    Q2.am_cn_call();
//    Q2.am_cn_put();
    
}

int main() {
    
    Q1();
    Q2();
 
    return 0;
}


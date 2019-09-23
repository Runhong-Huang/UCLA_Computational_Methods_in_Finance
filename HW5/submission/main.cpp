//
//  main.cpp
//  HW5
//
//  Created by HUANG runhong on 2/9/19.
//  Copyright © 2019 HUANG runhong. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Normal{
private:
    int length = 200;
    double seed = 522313;
    double Amiu = 0, Asigma = 1, Bmiu = 0, Bsigma = 1, rou = 0;
    
    
public:
    Normal(){};
    Normal(int num):length(num) {
    }
    
    ~Normal(){};
    
    
    // define two vectors with varaible length.
    std::vector<float> a ;
    std::vector<float> b ;
    
    void setseed(long s){
        this->seed = s;
    };
    
    //    set the mean and varaince of A
    void setA(float miu, float sigma){
        this->Amiu = miu;
        this->Asigma = sigma;
    };
    
    //    set the mean and varaince of B
    void setB(float miu, float sigma){
        this->Bmiu = miu;
        this->Bsigma = sigma;
    }
    
    void setCorr(float correlation){
        this->rou = correlation;
    }
    
    void generator(){
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
    };
};


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
    
//    for (int i = 0 ; i < n ; i ++){
//        for (int j =0  ; j < n ; j++){
//            cout << matrix_A[i][j] << " "   ;
//        }
//        cout << endl;
//    }
//
//    for (int i = 0 ; i < n ; i ++){
//        cout << vec_B[i] << endl;
//    }
    
    
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



class LSMC : public Normal {
    
private:
    

protected:
    float s0, r , k , sigma , t;
    int n = 100000;
    int exe_step = 200;
    int **index ;
    float **stock_price = new float *[n];
    
    
public:
    LSMC (float S0, float R, float K, float Sigma, float T ) : s0{S0}, r(R), k(K), sigma(Sigma),t(T){
    };
    ~LSMC(){};
    
    void stock_generator();
    float calculator(int n );
    virtual float lx(int term , float num);
    void set_s0 (float s);
    void set_time(float time);

};

void LSMC::set_s0 (float s){
    this -> s0 = s;
}

void LSMC::set_time (float time){
    this -> t = time;
}

void LSMC::stock_generator(){
    
    float dt = t/ exe_step;
//    cout << dt << endl;
    Normal norm(exe_step);
    norm.setCorr(-1.0);
    
//    float ** stock_price = new float *[n];
    
//    initializing the matrix : stock matrix, index matrix and ECV ( the expected continuation value)
    
//    generate the stock price
    for(int i = 0; i < n ; i = i+2){
        stock_price[i] = new float[exe_step+1];
        stock_price[i+1] = new float[exe_step+1];
        stock_price[i][0] = s0;
        stock_price[i+1][0] = s0;
        norm.setseed(i);
        norm.generator();
        for (int j =0 ; j < exe_step ; j++){
            stock_price[i][j+1] = stock_price[i][j] + r * stock_price[i][j] * dt + sigma * stock_price[i][j] * sqrt(dt) * norm.a[j];
            stock_price[i+1][j+1] = stock_price[i+1][j] + r * stock_price[i+1][j] * dt + sigma * stock_price[i+1][j]* sqrt(dt) * norm.b[j];
        }
    }
}


float LSMC::lx(int term , float num){
    
    try {
        if (term == 0 ){
            return (exp(-num/2.0));
//            return (1);
        }
        
        if (term == 1 ){
            return (exp(-num/2.0) * ( 1.0- num));
//            return (num);
        }
        
        if (term == 2 ){
//            return (num * num);
            return (exp(-num/2.0) * ( 1.0 - 2.0 * num + num*num /2.0 ));
        }
        
        if (term == 3 ){
//            return (num * num * num);
            return (exp(-num/2.0) * ( 1.0 - 3.0 * num + num*num * 1.5 - num * num * num /6.0));
        }
        else throw 1;
    }
    
    catch (int x ){
        cout << "Term number out of definition. " << endl;
    }
    
    return 0;
    
}



float LSMC::calculator(int order ){
    
    stock_generator();
    

    float dt = t/ exe_step;
    int** index = new int*[n];
    
    for(int i = 0; i < n; i++){
        index[i] = new int[exe_step+1];
    }
    
//    initialize the index matrix
    for (int i = 0 ; i < n ; i++){
        for (int j = 0 ; j < exe_step+1; j ++){
            index[i][j] = 0;
        }
    }
    
    
//    initialize the final step
    for (int i = 0 ; i < n; i++){
        float EV;
        EV = k - stock_price[i][exe_step];
        if (EV > 0.0){
            index[i][exe_step] = 1.0;
            
        }
    }
    
    vector<float> x;
    vector<float> y;
    
    
            cout << "We are using order of " << order << endl;
            
            for (int j = exe_step ; j > 0 ; j--){ //start from the second last coloumn
                
                for (int i = 0 ; i < n; i++){
                    if (k - stock_price[i][j]   > 0 ){
                        x.push_back(stock_price[i][j]);
                        float tmp = 0.0;
                        for ( int w = j+1; w <exe_step+1 ; w++){
                            if (index[i][w] == 1){
//                                cout << j << " " << i << " " << w<< endl;
//                                tmp = exp(-r * (w-j) *dt) * max((float) (k - stock_price[i][w]),(float) 0.0);
                                tmp = exp(-r * (w-j) *dt) * (k - stock_price[i][w]);
                            }
//                            cout << i << " "<< stock_price[i][j] <<" " << tmp << endl;
                        }
                        y.push_back(tmp);
                    }
                }

                
                float **A = new float *[order];
                float *B= new float [order];
                
                for (int i =0; i < order ; i ++){
                    B[i] = 0;
                    A[i] = new float [n];
                    for (int j =0; j < order ; j ++){
                        A[i][j] = 0;
                    }
                }
                
                //  generate matrix A nd B
                for (int i = 0; i < x.size(); i ++){
                    for (int row = 0; row < order; row++){
                        for (int col = 0; col < order; col++){
                            A[row][col] = A[row][col] + lx(row,x[i]) * lx(col, x[i]);
                        }
                        B[row] = B[row] + lx(row,x[i]) * y[i];
                    }
                }
            
//                regress to Ax = b, return the x (x is a vector)
                float *a = regression (A, B ,order);
                
//                cout << a[0] << " "<<  a[1] << endl;
                
//                compare the execise value with the continuation value
//                float count = 0;

                for (int i = 0 ; i < n; i++){
                    if ( k - stock_price[i][j] > 0  ){
                        float CV = 0;
                        for ( int w = 0 ; w < order ; w ++){
                            CV = CV + a[w] * lx(w, stock_price[i][j]);
                        }
//                        cout << CV << endl;
                        float EV = k - stock_price[i][j];

//                        remove the trailing zero
                        if (EV > CV) {
                            index[i][j] = 1 ;
                            for (int w = j+1 ; w < exe_step +1 ; w++){
                                index[i][w] = 0;
                            }
                        }
                    }
                }
                //            clean the index matrix
                x.clear();
                y.clear();
                //        cout << count << endl;
            }
    
//            now we have the matrix and we will calcultate the derivitie price
            float tmp = 0;
            for (int i = 0 ; i < n ; i++) {
                for (int j = 0 ; j < exe_step+1; j++){
                    if (index[i][j] == 1 ){
//                        cout << i << " " << j << endl;
                        tmp = tmp + exp(-r * j * dt) * (k - stock_price[i][j]) ;
//                        count++;
                    }
                }
            }
 
    return tmp/n;
}


class Monomials_LSMC : public LSMC {
    
private:
   
public:
    
    Monomials_LSMC(float S0, float R, float K, float Sigma, float T )  :  LSMC (S0, R,  K, Sigma, T ){};
    ~Monomials_LSMC(){};

    float lx(int term , float num);
    
};

float Monomials_LSMC::lx(int term , float num){
    
    try {
        if (term == 0 ){
            return (1);
           
        }
        
        if (term == 1 ){
            return (num);
        }
        
        if (term == 2 ){
            return (num * num);
        }
        
        if (term == 3 ){
            return (num * num * num);
        }
        else throw 1;
    }
    
    catch (int x ){
        cout << "Term number out of definition. " << endl;
    }
    return 0;
}

class Hermit_LSMC : public LSMC {
    
private:
    
public:
    
    Hermit_LSMC(float S0, float R, float K, float Sigma, float T )  :  LSMC (S0, R,  K, Sigma, T ){};
    ~Hermit_LSMC(){};
    
    float lx(int term , float num);
    
};

float Hermit_LSMC::lx(int term , float num){
    
    try {
        if (term == 0 ){
            return (1.0);
            
        }
        
        if (term == 1 ){
            return (2.0 * num);
        }
        
        if (term == 2 ){
            return (4.0 * num *num - 2.0);
        }
        
        if (term == 3 ){
            return (8.0 * num * num * num - 12.0 * num);
        }
        else throw 1;
    }
    
    catch (int x ){
        cout << "Term number out of definition. " << endl;
    }
    return 0;
}


class Forward_start : public Monomials_LSMC {
public:
    Forward_start (float S0, float R, float K, float Sigma, float T )  :  Monomials_LSMC (S0, R,  K, Sigma, T ){};
    ~Forward_start(){};
    float discret_calculator(float calltime);
    float continuous_calculator (float calltime, int order);
    
};

float Forward_start::discret_calculator(float calltime){
    
    stock_generator();
    
    int number = calltime * exe_step;
    float sum = 0;
    for (int i = 0 ; i < n ; i++){
        float tmp = stock_price[i][number] - stock_price[i][exe_step];
//        cout << tmp << endl;
        sum = sum + max (tmp , (float) 0.0 );
    }
    return (exp(-r * (t )) *sum/n);
}



float Forward_start::continuous_calculator(float calltime, int order){

    stock_generator();

    int num = calltime * exe_step ;
    float dt = t/ exe_step;
    int** index = new int*[n];

    for(int i = 0; i < n; i++){
        index[i] = new int[exe_step+1];
    }

    //    initialize the index matrix
    for (int i = 0 ; i < n ; i++){
        for (int j = 0 ; j < exe_step+1; j ++){
            index[i][j] = 0;
        }
    }


    //    initialize the final step
    for (int i = 0 ; i < n; i++){
        float EV;
        EV = stock_price[i][num] - stock_price[i][exe_step];
        if (EV > 0.0){
            index[i][exe_step] = 1;
        }
    }

    vector<float> x;
    vector<float> y;


    cout << "We are using order of " << order << endl;

    for (int j = exe_step-1 ; j > num ; j--){ //start from the second last coloumn
        for (int i = 0 ; i < n; i++){
            if (stock_price[i][num] - stock_price[i][j]   > 0 ){
                x.push_back(stock_price[i][num]-stock_price[i][j]);
                float tmp = 0.0;
                for ( int w = j+1; w <exe_step+1 ; w++){
                    if (index[i][w] == 1){
//                                                        cout << j << " " << i << " " << w<< endl;
                        //                                tmp = exp(-r * (w-j) *dt) * max((float) (k - stock_price[i][w]),(float) 0.0);
//                        tmp = exp(-r * (w-j) *dt) * (stock_price[i][num] - stock_price[i][w]);
                        tmp = exp(-r * (w-j) * dt) * (stock_price[i][num] - stock_price[i][w]);
//                        cout << i << " "<< w<< " "<< j << " "<< stock_price[i][j] <<" " << tmp << endl;
                    }

                }
                y.push_back(tmp);
            }
        }


        float **A = new float *[order];
        float *B= new float [order];

        for (int i =0; i < order ; i ++){
            B[i] = 0;
            A[i] = new float [n];
            for (int j =0; j < order ; j ++){
                A[i][j] = 0;
            }
        }

        //  generate matrix A nd B
        for (int i = 0; i < x.size(); i ++){
            for (int row = 0; row < order; row++){
                for (int col = 0; col < order; col++){
                    A[row][col] = A[row][col] + lx(row,x[i]) * lx(col, x[i]);
                }
                B[row] = B[row] + lx(row,x[i]) * y[i];
            }
        }

        //                regress to Ax = b, return the x (x is a vector)
        float *a = regression (A, B ,order);

//                        cout << a[0] << " "<<  a[1] << endl;
        //                compare the execise value with the continuation value

        for (int i = 0 ; i < n; i++){
             float CV = 0;
            if ( stock_price[i][num] - stock_price[i][j] > 0 ){
                for ( int w = 0 ; w < order ; w ++){
                    CV = CV + a[w] * lx(w, stock_price[i][num] - stock_price[i][j]);
                }
//                                        cout << CV << endl;
                float EV = stock_price[i][num] - stock_price[i][j];

                //                        remove the trailing zero
                if ((EV > CV ) ) {
                    index[i][j] = 1 ;
//                    cout << j << " "<< i << " "<< EV <<" " << CV << endl;
                    for (int w = j+1 ; w < exe_step +1 ; w++){
                        index[i][w] = 0;
                    }
                }
            }
        }
        //            clean the index matrix
        x.clear();
        y.clear();
        //        cout << count << endl;
    }

    //            now we have the matrix and we will calcultate the derivitie price
    float tmp = 0;
    for (int i = 0 ; i < n ; i++) {
        for (int j = 0 ; j < exe_step+1; j++){
            if (index[i][j] == 1 ){
                tmp = tmp + exp(-r * j * dt) * (stock_price[i][num] - stock_price[i][j]) ;
//                cout << i << " "<<j << endl;
//                cout << tmp << endl;
//                tmp = tmp + exp(-r * j * dt) * (stock_price[i][num] - stock_price[i][j]) ;
                //                        count++;
            }
        }
    }
    

    return tmp/n;


}


void Q1() {
    
    
    //    LSMC Q1(36.0, 0.06 ,40.0, 0.2, 0.5);
    ////    Q1.set_time(2);
    ////    The Basic Method use Laguerre method
    //    cout << "Using Laguerre Method, Time = 0.5, S0 = 40" << endl;
    //    cout << Q1.calculator(2)<< endl;
    //    cout << Q1.calculator(3)<< endl;
    //    cout << Q1.calculator(4)<< endl;
    
    
    
    //    LSMC lan(44.0, 0.06 ,40.0, 0.2, 0.5);
    //    Hermit_LSMC hemit(44.0, 0.06 ,40.0, 0.2, 0.5);
    Monomials_LSMC mono(40.0, 0.06 ,40.0, 0.2, 0.5);
    //    cout << mono.calculator(2)<< endl;
    //
    //
    cout << mono.calculator(2)<< endl;
    cout << mono.calculator(3)<< endl;
    cout << mono.calculator(4)<< endl;
    //
    //    mono.set_time(1);
    //    cout << "Using Monomials Method, Time = 0.5, S0 = 40" << endl;
    //    cout << mono.calculator(2)<< endl;
    //    cout << mono.calculator(3)<< endl;
    //    cout << mono.calculator(4)<< endl;
    //
    //    mono.set_time(2);
    //    cout << "Using Monomials Method, Time = 0.5, S0 = 40" << endl;
    //    cout << mono.calculator(2)<< endl;
    //    cout << mono.calculator(3)<< endl;
    //    cout << mono.calculator(4)<< endl;
    //
    //
    //
    
}



void Q2(){
    Forward_start fs(65.0, 0.06 ,60.0, 0.2, 1.0);
    cout << fs.discret_calculator(0.20) << endl;

    cout << fs.continuous_calculator(0.2,2) << endl;

    
}



int main() {
    
    Q1();
    Q2();

    
    cout << "end" << endl;
    return 0;
}


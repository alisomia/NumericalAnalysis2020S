#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <cstring>
#include <omp.h> 


#define M 32
#define n1 100000
#define tol 0.01


using namespace std;
int spin[2][M][M];
long long Hbar[2], Hsqbar[2],H[2];
double beta;
default_random_engine generator;
uniform_real_distribution<double> distribution(0.0, 1.0);
uniform_int_distribution<int> distribution2(0, M - 1);
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());


void init();
void warmingup();
void metropolis();

int main(){
    fstream f;
    f.open("result5.txt", ios::out);  // save the result
    init();  // initilization 
    for(int B = 1; B <99; B++){
        beta = 0.01*B;
        Hbar[0] = Hbar[1] = Hsqbar[0] = Hsqbar[1] = 0;  //simulate two paths
        int test = 0;
        warmingup();  // warming up 
        do{
            test+=1; 
            metropolis();
            //cout << Hbar[1] - Hbar[0] << endl;
            if(test==1000){Hbar[1] = Hbar[0] = Hsqbar[1] = Hsqbar[0] = 0; test = 0;}  // retry until the result converges
        }while(test==0 || abs(Hbar[1]-Hbar[0])>tol*M*M*10000*test || abs(Hsqbar[1]-Hsqbar[0])>tol*M*M*10000*test);
    double U = 0.5*(Hbar[1]+Hbar[0]), C = 0.5*(Hsqbar[0]+Hsqbar[1]);
    U/=10000*test; C = (C/(10000*test)-U*U)*beta*beta/M/M; U/=M*M;
    cout << beta << ' ' << U << ' '<< C << ' '<< test << endl;
    f << beta << ' ' << U << ' '<< C << ' '<< test << endl;
 }
    f.close();
}

void init(){
    for(int k = 0; k<2;k++){
    H[k] = 0;
    for(int i = 0; i<M;i++)
    for(int j = 0; j<M;j++)
    spin[k][i][j] = 1;
    for (int i = 0; i < M; i++)
    for (int j = 0; j < M; j++)
        H[k] += spin[k][i][j] * spin[k][(i + 1) % M][j]+ spin[k][i][j] * spin[k][i][(j + 1) % M];
    }
}


void warmingup(){
    for(int k = 0; k < 2; k++){
    for(int iter = 0; iter < n1; iter++){
        int i = rng()%M, j = rng()%M, delta = 2*spin[k][i][j] * (spin[k][(i + M - 1) % M][j] + spin[k][(i + 1) % M][j] + spin[k][i][(j + M - 1) % M] + spin[k][i][(j + 1) % M]);
        if (1.0 * rng() / (0xffffffff) < exp(-beta*delta)){
            spin[k][i][j] = -spin[k][i][j], H[k] -= delta;
        }
    }
    }
}
void metropolis(){
    for(int k = 0; k < 2; k++){
    for(int iter = 0; iter < 10000; iter++){
        int i = rng()%M, j = rng()%M, delta = 2*spin[k][i][j] * (spin[k][(i + M - 1) % M][j] + spin[k][(i + 1) % M][j] + spin[k][i][(j + M - 1) % M] + spin[k][i][(j + 1) % M]);
        if (1.0 * rng() / (0xffffffff) < exp(-beta*delta)){
            spin[k][i][j] = -spin[k][i][j], H[k] -= delta;
        }
        Hbar[k] -=  H[k], Hsqbar[k] +=  H[k] * H[k];
    }
    }
}

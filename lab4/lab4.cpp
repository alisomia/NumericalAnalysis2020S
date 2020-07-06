/*
This code implements FFT algorithm and FFT-based Filter.
Author： Ting Lin @PKU
*/
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstring>

using namespace std;
const double PI=3.1415926535897384626;
struct com{  // Complex 
    double re, im;
    com(double rr, double ii){this->re = rr; this->im = ii;}
    com(){this->re=this->im=0;}

};
com operator +(com a,com b){return com (a.re+b.re,a.im+b.im);}
com operator -(com a,com b){return com (a.re-b.re,a.im-b.im);}
com operator *(com a,com b){return com (a.re*b.re-a.im*b.im, a.im*b.re+a.re*b.im);}
com operator /(com a, double b){return com (a.re/b, a.im/b);}
com operator *(com a, double b){return com (a.re*b, a.im*b);}


inline com root(int n){return com(cos(2*PI/n),sin(2*PI/n));}  // n-th root
inline com inv_root(int n){return com(cos(2*PI/n),-sin(2*PI/n));} // n-th root inverse
inline double f(double t){return exp(-t*t/10)*(sin(2*t)+2*cos(4*t)+0.4*sin(t)*sin(50*t));} //function
const int maxn=256;
com a[maxn];
double b[maxn];

void change(com *a,int k){ // re-oder the array.
	for(int i=0,j=0;i<k;i++){
		if(i<j)swap(a[i],a[j]);
		for(int l=(k>>1);(j^=l)<l;l>>=1);
	}
}
void fft(com *a,int n){ // assembling 
	for(int m=2;m<=n;m<<=1){
		com e=com(cos(2*PI/m),sin(2*PI/m));
		for(int k=0;k<n;k+=m){
			com w=com(1,0);
			for(int j=0;j<m/2;j++){
				com t=w*a[k+j+(m/2)];
				com u=a[k+j];
				a[k+j]=u+t;
				a[k+j+(m/2)]=u-t;
				w=w*e;
			}
		}
	}
}
void inv_fft(com*a,int n){
	for(int m=2;m<=n;m<<=1){
		com e=inv_root(m);
		for(int k=0;k<n;k+=m){
			com w=com(1,0);
			for(int j=0;j<m/2;j++){
				com t=w*a[k+j+(m/2)];
				com u=a[k+j];
				a[k+j]=u+t;
				a[k+j+(m/2)]=u-t;
				w=w*e;
			}
		}
	}
	for(int i=0;i<n;i++) a[i]=a[i]/n;
}
void filter(com *a,int n, int m){
	change(a,n);
	fft(a,n);
	for(int i=m; i<=n-m;i++) a[i] = com(0,0);

	change(a,n);
	inv_fft(a,n);
}
int main(){
    
    //for(int i=0;i<256;i++) b[i] = f(PI*i/256);
    for(int m:{0,1,2,3,4,5,6,7,8,9,10,20,30,40,50}){
    for(int i=0;i<256;i++) a[i] = com(f(PI*i/256),0);
    filter(a,256,m);
    for(int i=0;i<256;i++) b[i] = a[i].re;
    ofstream f ("./out"+to_string(m), ios::out);
    for(int i=0;i<256;i++) f << b[i] << endl;
    f.close();
    }
	return 0;
}
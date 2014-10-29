#include<iostream>
#include<cmath>
#include<vector>
#include<climits>
#include<stdlib.h>
using namespace std;

/*
Author: Qing Li
E-Mail: qingli@cse.tamu.edu

Description: This file implements density evolution for BEC channel. Also density evolution 
for joint decoder of content-replicated coders. (refer to my paper).
Date: 10/29/2014  
*/

//get lambda_i
double getLambda_i(int dv, int dc){
  return 1;
  //return (double)(dc - dv)/(double)dc;
}

//get lambda_p
double getLambda_p(int dv, int dc){
  return 1;
  //return (double) dv/(double)dc;
}

//computer lambda_i(x)
double lambda_iX(int dv, int dc, double x){
  double coeff = getLambda_i(dv, dc);
  double power = pow(x,  2*dv - 1);
  return coeff*power;
}

//computer nchoosek
int nChoosek(int n, int k){
  if(k > n) return 0;
  if(k*2 > n) k = n - k;
  if(k==0) return 1;

  int result = n;
  for(int i = 2; i <=k; i++){
    result *= (n-i+1);
    result /=i;
  }
  return result;
}

//computer lambda_p(x)
double lambda_pX(int dv, int dc, double x){
  double coeff = getLambda_p(dv,  dc);
  double power = pow(x, dv - 1);
  return coeff*power;
}

double lambda(int dv, int dc, double x){
  return pow(x, dv - 1);
}

//get coefficients, rho^{(i)}_{i,j}
void coeffRhoI(int dv, int dc, vector<double> &coeff){
  double sum = 0;
  double con = pow((double)dv/dc, dc);
  for(int i = 1; i <=dc; i++){
    //cout<<"nchooseK ("<<dc<<", "<< i<<") is "<<nChoosek(dc, i)<<endl;
    //cout<<pow((double)(dc - dv)/dc, i)<<endl;
    //cout<<pow((double)dv/dc, dc - i )<<endl;
    double temp =(double)((nChoosek(dc, i))*pow((double)(dc - dv)/dc, i)*pow((double)dv/dc, dc - i))/(1-con);
    //cout<<"The " <<i<<" -th coefficient is "<<temp<<endl;
    sum += temp;
    //cout<<"Now sum is "<<sum<<endl; 
    coeff.push_back(temp);
  }
  //cout<<"==================="<<endl;
}

//compute Rho^{(i)} 
double RhoI(int dv, int dc, double x, double y){
  double result = 0.0;
  vector<double> coeff;
  //cout<<"x  is " <<x<<endl; 
  coeffRhoI(dv, dc, coeff);
  double sum = 0;
  for(int i = 0; i <dc; i++){
    // cout<<i<<"-th coefficient is "<<coeff[i]<<endl;
     sum += coeff[i];
     //cout<<"pow one is "<<pow(x, i)<<endl;
     //cout<<"pow two is "<<pow(y, dc - i - 1)<<endl;
     //cout<<"sum of coefficient is "<<sum<<endl;
    result += coeff[i]*pow(x, i)*pow(y,dc-i-1);
  }
  return result;
}
//get coefficients, rho^{(p)}_{i,j}
void coeffRhoP(int dv, int dc, vector<double> &coeff){
  double sum = 0;
  double con = pow((double)(dc - dv)/dc, dc);
  for(int i = 0; i <dc; i++){
    //cout<<"nchooseK ("<<dc<<", "<< i<<") is "<<nChoosek(dc, i)<<endl;
     //cout<<pow((double)(dc - dv)/dc, i)<<endl;
    //cout<<pow((double)dv/dc, dc - i)<<endl;
     double temp = (double)(nChoosek(dc, i))*pow((double)(dc - dv)/dc, i)*pow((double)dv/dc, dc - i)/(1-con);
    sum += temp;
    //cout<<"The " <<i<<"-th coefficient is "<<temp<<endl;
    //cout<<"Now sum is "<<sum<<endl;
    coeff.push_back(temp);
  }
}

//compute Rho^{(p)}
double RhoP(int dv, int dc, double x, double y){
  double result = 0.0;
  vector<double> coeff;
  coeffRhoP(dv, dc, coeff);
  for(int i = 0; i <dc; i++){
    //cout<<i<<"-th coefficient is "<<coeff[i]<<endl;
    //cout<<"pow one is "<<pow(x, i)<<endl;
    //cout<<"pow two is "<<pow(y, dc - i -1)<<endl;
    result += coeff[i]*pow(x, i)*pow(y, dc- 1 -i);
  }
  return result;
}


int f(int dv, int dc, double e){
  double x = e;
  for(int i = 0; i < 1000; i++){
    long double xnow = e*lambda(dv, dc, 1- RhoI(dv, dc, 1 - x, 1 - x));
    //cout<<"now the probability is "<<x<<endl;
    x = xnow;			      
  }
  if(x < (double)1/pow(10.0, 32)){
    return 1;
  }
  else return -1;
}

int fi(int dv, int dc, double e){
  long double x = e*e;
  long double y = e;
  for(int i = 0; i < 1000; i++){
    //cout<<"x(e) is "<<x<<"        ";
    //cout<< RhoI(dv, dc, 1 - x, 1- y)<<endl;
    //cout<< 1 - RhoI(dv, dc, 1 - x, 1- y)<<endl;
    //cout<<lambda_iX(dv, dc, 1 - RhoI(dv, dc, 1 - x, 1- y))<<endl;
    long double xnow = (e*e)*lambda_iX(dv, dc, 1 - RhoI(dv, dc, 1 - x, 1- y));
    //cout<<"y(e) is "<<y<<endl;
     //cout<< RhoP(dv, dc, 1 - x, 1- y)<<endl;
     //cout<< 1 - RhoP(dv, dc, 1 - x, 1- y)<<endl;
     //cout<<lambda_pX(dv, dc, 1 - RhoP(dv, dc, 1 - x, 1- y))<<endl;
    long double ynow = e*lambda_pX(dv, dc, 1- RhoP(dv, dc, 1 - x, 1 - y));
    if(ynow == y){
      cout<<"converge "<<endl;
      break;
    }
    else{
      x = xnow;
      y = ynow;
      cout<<"now f_i is "<<x<<"    ";
      cout<<"now f_p is "<<y<<endl;
    }
  }
  if(y < (double)1/pow(10.0, 60))
    {
      //cout<<x<<endl;
    //cout<<e<<" is good"<<endl; 
      return 1;
    }
  else return -1;
}

int fi2(int dv, int dc, int dv3, int dc3, double e){
  long double x = e*e;
  long double y = e;
  for(int i = 0; i < 10000; i++){
    long double xnow = (e*e)*lambda_iX(dv, dc, 1 - RhoI(dv, dc, 1 - x, 1 - y))*lambda(dv3, dc3, 1 - RhoI(dv3, dc3, 1 - x, 1 - x));
    long double ynow = e*lambda_pX(dv, dc, 1 - RhoP(dv, dc, 1 - x, 1- y));
    if(xnow == x){
      break;
    }
    else{
      x = xnow;
      y = ynow;
    }
  }
  if(x<(double)1/pow(10.0,60)){
      return 1;
    }
  else 
    return -1;
}

int fp(int dv, int dc, double e){
  long double x = e*e;
  long double y = e;
  for(int i = 0; i < 1000; i++){
    long double xnow = (e*e)*lambda_iX(dv, dc, 1 - RhoI(dv, dc, 1 - x, 1- y));
    long double ynow = e*lambda_pX(dv, dc, 1- RhoP(dv, dc, 1 - x, 1 - y));
    if(ynow == y)
      break;
    else{
      x = xnow;
      y = ynow;
    }
    //cout<<"now f_i is "<<x<<"     ";
    //cout<<"now f_p is "<<y<<endl;
  }
  if(y < (double)1/pow(10.0, 60))
  {
    //cout<<y<<endl;
    //cout<<e<<" is good"<<endl; 
      return 1;
    }
  else return -1;
}


double bp(int dv, int dc){
  for(double e = 0; e < 1; e = e + 0.0001){
    if(f(dv, dc, e) == 1) continue;
    else{
      //cout<<"Threshold for conventional code ("<<dv<<"," <<dc<<") is "<<e<<endl;
      return e;
    }
  }
}
//bp for i
double bpI(int dv, int dc){
  for(double e = 0; e < 1; e = e + 0.0001){
    if(fi(dv, dc, e) == 1) continue;
    else{
      //cout<<"threshold for information bit is for ( "<<dv<<" "<<dc<<" ) is "<<e<<endl;
      return e;
    }
  }
}

//bp for p
double bpP(int dv, int dc){
  for(double e = 0; e < 1; e = e + 0.0001){
    if(fp(dv, dc, e) == 1) continue;
    else{
      //cout<<"threshold for check node is for ( "<<dv<<" "<<dc<<" ) is "<<e<<endl;
      return e;
    }
  }
}


double bpAnother(int dv, int dc, int dv3, int dc3){
  for(double e = 0; e < 1; e = e + 0.0001){
    if(fi2(dv, dc, dv3, dc3,e)==1) continue;
    else
      return e;
  }
}

int main(int argc, char *argv[]){
  int dv  = atoi(argv[1]);
  int dc  = atoi(argv[2]);
  //int dv3 = atoi(argv[3]);
  //int dc3 = atoi(argv[4]);

  /*cout<<"test for coeffRhoI"<<endl;
  vector<double>coeff;
  double temp = 0;
  coeffRhoI(dv, dc, coeff);
  for(int i = 0; i < coeff.size(); i++){
    cout<<coeff[i]<<" ";
    temp += coeff[i];
  }
  cout<<endl<<"sum is "<<temp<<endl;
  cout<<"==========================="<<endl;*/

  //cout<<"test for RhoI"<<endl;
  //cout<<lambda_iX(dv, dc, 1- RhoI(dv, dc, 0.1, 0.1))<<endl; 
  //cout<<"=================="<<endl;
   //cout<<"test for RhoP"<<endl;
   //cout<<RhoP(dv, dc, 0.9, 0.9)<<endl;
  //cout<<"test for lambda"<<endl;
  //cout<<lambda_iX(dv, dc, 0.5)<<endl;
  //cout<<lambda_pX(dv, dc, 0.5)<<endl;

  //cout<<"traditional BP threshold for ("<<dv<<" , "<<dc<<") is "<<bp(dv, dc)<<endl;
  //cout<<"BP threshold for parity check bit is "<<bpP(dv, dc)<<endl;
  //cout<<"BP threshold for information bit is "<<bpAnother(dv, dc, dv3, dc3)<<endl;
  //cout<<"==================="<<endl;
  //f(dv, dc, 0.8);
  fi(dv, dc, 0.62);
  //bp(dv, dc);
  //cout<<"BP threshold for information bit is "<<bpI(dv, dc)<<endl;
  //bpP(dv, dc);
}

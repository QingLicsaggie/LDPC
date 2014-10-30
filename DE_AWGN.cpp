//Gaussian Approximation Density Evolution: version 5 Based on Analysis of Sum-Product Decoding of Low-Density Parity-Check Codes Using a Gaussian Approximation.
//Author: Zhu, Shenli
//        Qing Li (qingli@cse.tamu.edu)
//Change Log: 
//	7/24	Version 2 change linear search to binary search (bsearch)
//      8/7     Version 3 bug fixed and add usage()
//      8/8     Version 4 add input file
//      10/30/2014 Version 5 made it more readable.
//Output:
//(1) threshold
//(2) average column weight
extern "C"
{
    #include <stdlib.h>
    #include <assert.h>
    #include <string.h>
}
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
//constant used in Phi()
//refer to the paper page 661 for the meanings of those parameters.
static const double ALFA=-0.4527;
static const double BETA= 0.0218;
static const double GAMA= 0.86;
static const double PI= 3.1415926535897932384626433832795;

//iteration parameter
static const int ITERATION= 1000;
static const double SIGMA_BEGIN= 0.2;
static const double SIGMA_PRECISION= 0.0001;

//check irregular: 1;
//check regular: 0;
bool chk_irre;

// parameter used in phi_inv
// phi_inv_sample
static const int PHI_SAMPLE_LENGTH = 100000;
static const double PHI_SAMPLE_RANGE = 40;
static const double PHI_SAMPLE_STEP = PHI_SAMPLE_RANGE/PHI_SAMPLE_LENGTH;

// phi_inv function in and out vector
double phi_inv_in[PHI_SAMPLE_LENGTH]={0};
double phi_inv_out[PHI_SAMPLE_LENGTH]={0};

//---------------------------------
//regular LDPC degree
//---------------------------------
int dv,dc;

double code_rate;

//----------------------------------
//irregular LDPC degree distribution
//----------------------------------
int vdd_cnt=0;
int cdd_cnt=0;

//edge perspective
double vdd_e[20][2]={0};
double cdd_e[20][2]={0};
//node perspective
double vdd_n[20][2]={0};
double cdd_n[20][2]={0};

//-------------------------
//function declearation
//-------------------------
void usage();
int check_threshold(double);
double re_mu_update (double,double);
double irre_mu_update (double,double);
double Phi(double);
double Phi_inv(double);
int compare (const void *, const void *);
//some transfer function
void vdd_n2vdd_e();
void cdd_n2cdd_e();
void vdd_e2vdd_n();
void cdd_e2cdd_n();
void vdd_n2cdd_n();
void vdd_e2cdd_e();
void print_degree_dis();
void check_degree_dis();

void usage()
{
    cout<<"usage: GADE_ALL_NEW"<<endl;
    cout<<"the file named argv.dat must be contain in the same fold of the program"<<endl;
    cout<<"regular code example:"<<endl;
    cout<<"---------------"<<endl;
    cout<<"-r"<<endl;
    cout<<"3"<<endl;
    cout<<"6"<<endl;
    cout<<"---------------"<<endl<<endl;
    cout<<"irregular code example:"<<endl;
    cout<<"---------------"<<endl;
    cout<<"-i-e or -i-n"<<endl;
    cout<<"3(vdd counter)"<<endl;
    cout<<"2 0.3"<<endl;
    cout<<"3 0.4"<<endl;
    cout<<"4 0.3"<<endl;
    cout<<"3(cdd counter)"<<endl;
    cout<<"6 0.8"<<endl;
    cout<<"7 0.2"<<endl;
    cout<<"---------------"<<endl;
    cout<<"-i-e-r -i-n-r"<<endl;
    cout<<"0.5(code rate)"<<endl;
    cout<<"3(vdd counter)"<<endl;
    cout<<"2 0.3"<<endl;
    cout<<"3 0.4"<<endl;
    cout<<"4 0.3"<<endl;
    cout<<"---------------"<<endl;
    exit(1);
}

//----------------
//this main is DE
//----------------
int main(int argc,char** argv){
    ifstream is1;
    int i;
    
    char method[10]; 
    
    //check the inputs
    char* argvfile="argv.dat";
    is1.open(argvfile);
    
    is1>>method;
    
    if(strcmp(method,"-r")==0) {
        is1>>dv>>dc;
        chk_irre=0;
    }
    else if (strcmp(method,"-i-n")==0) {
        is1>>vdd_cnt;
        for(i=0;i<vdd_cnt;i++) {
            is1>>vdd_n[i][0]>>vdd_n[i][1];
        }
        is1>>cdd_cnt;
        for(i=0;i<cdd_cnt;i++) {
            is1>>cdd_n[i][0]>>cdd_n[i][1];
        }
        chk_irre=1;
        vdd_n2vdd_e();
        cdd_n2cdd_e();
    }
    else if (strcmp(method,"-i-e")==0) {
        is1>>vdd_cnt;
        for(i=0;i<vdd_cnt;i++) {
            is1>>vdd_e[i][0]>>vdd_e[i][1];
        }
        is1>>cdd_cnt;
        for(i=0;i<cdd_cnt;i++) {
            is1>>cdd_e[i][0]>>cdd_e[i][1];
        }
        chk_irre=1;
		vdd_e2vdd_n();
		cdd_e2cdd_n();
    }
    else if (strcmp(method,"-i-n-r")==0) {
        is1>>code_rate;
		is1>>vdd_cnt;
        for(i=0;i<vdd_cnt;i++) {
            is1>>vdd_n[i][0]>>vdd_n[i][1];
        }
        cdd_cnt=2;
		vdd_n2cdd_n();
        chk_irre=1;
        vdd_n2vdd_e();
		cdd_n2cdd_e();
    }
    else if (strcmp(method,"-i-e-r")==0) {
		is1>>code_rate;
        is1>>vdd_cnt;
        for(i=0;i<vdd_cnt;i++) {
            is1>>vdd_e[i][0]>>vdd_e[i][1];
        }
        cdd_cnt=2;
		vdd_e2cdd_e();
        chk_irre=1;
        vdd_e2vdd_n();
		cdd_e2cdd_n();
    }
    else {
        usage();
    }
            
    // calc the phi_inv data
    // the data will be used in the program's lifetime
    for (i=0;(PHI_SAMPLE_RANGE-i*PHI_SAMPLE_STEP)>0;i++) { 
    	phi_inv_out[i]=PHI_SAMPLE_RANGE-i*PHI_SAMPLE_STEP;
		phi_inv_in[i]=Phi(phi_inv_out[i]);
    }		
    
	double sigma;
	for(sigma=SIGMA_BEGIN;;sigma=sigma+SIGMA_PRECISION) {
		if (check_threshold(sigma)==0)
			break;
	}
	
	print_degree_dis();
	check_degree_dis();
	cout<<((chk_irre==0)?"regular":"irregular")<<" code threshold"<<endl;
	cout<<(sigma-SIGMA_PRECISION)<<endl<<endl;
	return 1;
}

//check if the threshold is satisfied
int check_threshold(double thre) {
	assert(thre>=0);
	double mu0 = 2/(thre*thre);
	double mu;
	double mu_next=0;
	for (int j=0;j<ITERATION;j++) {
		mu = mu_next;
		mu_next = (chk_irre == 0)?
		                            re_mu_update(mu0,mu)//check regular LDPC
		                         :
		                            irre_mu_update(mu0,mu);//check irregular LDPC
		if (mu_next >= (PHI_SAMPLE_RANGE-1)) {
			return 1;
		}
	}
	return 0;
}

//regular mu update: with constant dv and dc
double re_mu_update (double mu0,double mu) {
	assert(mu0>=0);
	assert(mu>=0);
	double temp1=0;
	double temp2=0;
	temp1 = 1-Phi(mu0+(dv-1)*mu);
	temp2 = 1-pow(temp1,(dc-1));
	return Phi_inv(temp2);
}

//irregular mu update: with different dv and different dc
//the output is wrong
double irre_mu_update (double mu0,double mu){
	assert(mu0>=0);
	assert(mu>=0);
	
	double temp1=0;
	double temp2=0;
	
	for (int v=0;v<vdd_cnt;v++) {
		temp1 += vdd_e[v][1]*Phi(mu0+(vdd_e[v][0]-1)*mu);
	}
	temp1 = 1 - temp1;
	
	for (int c=0;c<cdd_cnt;c++) {
		temp2 += cdd_e[c][1]*Phi_inv(1-pow(temp1,(cdd_e[c][0]-1)));
	}
	return temp2;
}

//Phi function
double Phi(double x) {
	//assert x>0
	assert(x>=0);

    if(x<=10) {
		return exp(ALFA*pow(x,GAMA)+BETA);
	}
	else if(x>10) {
		return sqrt(PI/x)*exp(-x/4)*(1+1/(14*x)-3/(2*x));
	}
}

//inverse Phi function v2: using bsearch
double Phi_inv(double y) {
    assert(y>=0);
    
	double * pItem;
	int index;
	pItem = (double*) bsearch (&y, phi_inv_in, PHI_SAMPLE_LENGTH, sizeof (double), compare);
    index = pItem-phi_inv_in;
    return  (pItem!=NULL) ? (phi_inv_out[index]): PHI_SAMPLE_RANGE;
}
    
// used in Phi_inv()'s bsearch()
int compare (const void * kp, const void * vp)
{
    //key is the item I want to search
    double key = *(const double*)kp;
    const double *vec = (const double*)vp;
    
	//////////////////////////////////////////////////
	//	vec[-1]            vec[0]               vec[1] 
	//     |                |                       |
	//              |                     |
	//  (vec[-1]+vec[0])/2.0        (vec[0]+vec[1])/2.0
	//////////////////////////////////////////////////
	if (key > (vec[0]+vec[1])/2.0)
		return 1;
	else if (key < (vec[-1]+vec[0])/2.0)
		return -1;
	else
		return 0;
}

//transfer node perspective degree distribution to edge perspective
//vdd_n->vdd-e
void vdd_n2vdd_e(){
	int i;
    double vdd_n_sum=0;
    for(i=0;i<vdd_cnt;i++) {
        vdd_n_sum += vdd_n[i][0]*vdd_n[i][1];
    }
    for(i=0;i<vdd_cnt;i++) {
        vdd_e[i][0]=vdd_n[i][0];
        vdd_e[i][1] = (vdd_n[i][0]*vdd_n[i][1])/vdd_n_sum;
    } 
}
//cdd_n->cdd_e	
void cdd_n2cdd_e (){
    int i;
    double cdd_n_sum=0;
    for(i=0;i<cdd_cnt;i++) {
        cdd_n_sum += cdd_n[i][0]*cdd_n[i][1];
    }
    for(i=0;i<vdd_cnt;i++) {
        cdd_e[i][0]=cdd_n[i][0];
        cdd_e[i][1] = (cdd_n[i][0]*cdd_n[i][1])/cdd_n_sum;
    }     
}

//transfer edge perspective degree distribution to node perspective
//vdd_e->vdd_n
void vdd_e2vdd_n() {
	int i;
	double vdd_e_sum=0;
    for(i=0;i<vdd_cnt;i++) {
        vdd_e_sum += (vdd_e[i][1]/vdd_e[i][0]);
    }
    for(i=0;i<vdd_cnt;i++) {
        vdd_n[i][0]=vdd_e[i][0];
        vdd_n[i][1] = (vdd_e[i][1]/vdd_e[i][0])/vdd_e_sum;
    } 
}

//cdd_e->cdd_n
void cdd_e2cdd_n (){
    int i;
    double cdd_e_sum=0;
    for(i=0;i<cdd_cnt;i++) {
        cdd_e_sum += cdd_e[i][1]/cdd_e[i][0];
    }
    for(i=0;i<vdd_cnt;i++) {
        cdd_n[i][0]=cdd_e[i][0];
        cdd_n[i][1] = (cdd_e[i][1]/cdd_e[i][0])/cdd_e_sum;
    }     
}

//vdd_n -> cdd_n
void vdd_n2cdd_n() {
	double vdd_sum=0;
	double cdd_sum=0;
	int i;
	for(i=0;i<vdd_cnt;i++) {
		vdd_sum += vdd_n[i][0]*vdd_n[i][1];
	}
	cdd_sum = vdd_sum/(1-code_rate);
	cdd_n[0][0]=(int) cdd_sum;
	cdd_n[1][0]=cdd_n[0][0]+1;
	cdd_n[1][1]=cdd_sum-(double)cdd_n[0][0];
	cdd_n[0][1]=1-cdd_n[1][1];
}

//vdd_e->vdd_n->cdd_n->cdd_e
void vdd_e2cdd_e() {
	
	vdd_e2vdd_n();
	vdd_n2cdd_n();
	cdd_n2cdd_e();
	
}

void print_degree_dis() {
	int i;
	double vdd_sum=0;
	double cdd_sum=0;
	cout<<"degree distribution(edge perspective):"<<endl;
	cout<<"Varible Node:"<<endl;
	for(i=0;i<vdd_cnt;i++) {
		cout<<vdd_e[i][0]<<":\t"<<vdd_e[i][1]<<endl;
	}
	cout<<"Check Node:"<<endl;
	for(i=0;i<cdd_cnt;i++) {
		cout<<cdd_e[i][0]<<":\t"<<cdd_e[i][1]<<endl;
	}
	cout<<endl;
	cout<<"degree distribution(node perspective):"<<endl;
	cout<<"Varible Node:"<<endl;
	for(i=0;i<vdd_cnt;i++) {
		cout<<vdd_n[i][0]<<":\t"<<vdd_n[i][1]<<endl;
	}
	cout<<"Check Node:"<<endl;
	for(i=0;i<cdd_cnt;i++) {
		cout<<cdd_n[i][0]<<":\t"<<cdd_n[i][1]<<endl;
	}
	cout<<endl;

	for(i=0;i<vdd_cnt;i++) {
		vdd_sum += vdd_n[i][0]*vdd_n[i][1];
	}
	cout<<"Average Column Weight:"<<vdd_sum<<endl;
	for(i=0;i<cdd_cnt;i++) {
		cdd_sum += cdd_n[i][0]*cdd_n[i][1];
	}
	cout<<"Average Row Weight:"<<cdd_sum<<endl<<endl;
}

void check_degree_dis() {
	int i;
	double vdd_n_sum=0;
	double cdd_n_sum=0;
	double vdd_e_sum=0;
	double cdd_e_sum=0;
	
	for (i=0;i<vdd_cnt;i++) {
		vdd_n_sum += vdd_n[i][1];
	}
	for (i=0;i<cdd_cnt;i++) {
		cdd_n_sum += cdd_n[i][1];
	}	
	for (i=0;i<vdd_cnt;i++) {
		vdd_e_sum += vdd_e[i][1];
	}
	for (i=0;i<cdd_cnt;i++) {
		cdd_e_sum += cdd_e[i][1];
	}
	if (vdd_n_sum<0.99 || vdd_n_sum>1.01 || cdd_n_sum<0.99 || cdd_n_sum>1.01
	  ||vdd_e_sum<0.99 || vdd_e_sum>1.01 || cdd_e_sum<0.99 || cdd_e_sum>1.01)
		cout<<"Warning: The sum of Degree Distribution not equal to 1!"<<endl<<endl;
}

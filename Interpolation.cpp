#include "Interpolation.h"
#include "eigen/Eigen/Dense"
#include "BivariatePolynomial.h"
#include <fstream>

using namespace std;
using namespace Eigen;


Interpolation::~Interpolation()
{
}
Interpolation::Interpolation(uint64_t b, MatrixXd B):k(B.rows()),deg(b),rank(0)
{
	d=((b+1)*(b+2)/2)-1;
	A=B;
	this->P = new double*[d+1];
	for (int i = 0; i < d+1; ++i) {
		this->P[i] = new double[d+1];
		for (int j = 0; j < d+1; ++j) {
        		this->P[i][j] = 0;
        	}
	}	
}


double Interpolation::powerOf(double Num, uint64_t power)
{	
	double p=1.0;
	if(power!=0){
		for (int i=0; i<power; ++i) {
		        p=p*Num;
		}
	}
	return p;	
}

double Interpolation::entry(uint64_t i, uint64_t a, uint64_t b)
{	double x,y,m,n;
	x=A(i,0);
	y=A(i,1);
	m=this->powerOf(x,a);
	n=this->powerOf(y,b);
	m=m*n;
	if(m==-0.0) m=0.0;
	return m;
}

MatrixXd & Interpolation::ComputeM()
{
	this->M = MatrixXd(k,k+1);
	int i,j,a,b,sum,y;
	i=0;
	while(i<k){
		y=0;
		sum=0;
		a=sum;
		while(y!=k+1){
			b=sum-a;
			M(i,y)=this->entry(i,a,b);
			y++;
			if(b==sum){
				sum++;
				a=sum;
			}
			else a--;
		}
		i++;
	}
	return this->M;
}


int Interpolation::Check_k()
{
	if(k==d){
		//cout<<"The system of equations is well-constrained."<<endl;
		return 1;
	}
	else if(k>d){
		cout<<"Overconstrained system, usually no solution."<<endl;
		return 0;
	}
	else{
		cout<<"Underdefined, usually infinite number of solutions."<<endl;
		return 0;
	}
}


BivariatePolynomial* Interpolation::find(ofstream & equationsTxt)
{
	if(this->Check_k()!=0){
		BivariatePolynomial* Bp;
		this->ComputeM();
		cout << "M is : "<<endl;
		cout << M << endl;
		FullPivLU<MatrixXd> lu_decomp(M);
		rank=lu_decomp.rank();
//		cout<<"The system of equations is well-constrained."<<endl;
		//cout << "The rank of M is " << rank << endl;
		if(rank==k){
			MatrixXd ker = M.fullPivLu().kernel();
			double t=ker(k);
			ker=ker/t;
			cout << "The kernel of M is "<< endl << ker << endl;
			int count=0;
			this->P[0][0] = ker(0);
			int a,b,sum=0;
			a=sum;
			cout<< "Polynomial : ";
			while(count!=k+1){
				b=sum-a;
				t=ker(count);
				this->P[b][a]=t;
				if(t>0.0 && count>0){
					equationsTxt<<"+";
					cout<<"+";
				}
				if(count>0) {
					equationsTxt<<"";
					cout<<"";
				}
				if(t!=0.0){
					equationsTxt<<t;
					cout<<t;
				if(a==1) {
					equationsTxt<<"*x";
					cout<<"*x";
				}
				if(b==1){
					equationsTxt<<"*y";
					cout<<"*y";
				}
				if(a>1) {
					equationsTxt<<"*x^"<<a;
					cout<<"*x^"<<a;
				}
				if(b>1){
					equationsTxt<<"*y^"<<b;
					cout<<"*y^"<<b;
				}
			}
			count++;
			if(b==sum){
				sum++;
				a=sum;
			}
			else a--;
		}
			equationsTxt<<endl;
			cout<<endl;
			Bp = new BivariatePolynomial(deg,P);
		}
		else {
			cout<<"The system is infeasible or the solution is numerically unstable."<<endl;
			return NULL;
		}
		return Bp;
	}
//	cout<<"K is not proper"<<endl;
	return NULL;

}





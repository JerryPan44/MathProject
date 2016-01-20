#include <cstdlib>
#include <iostream>
#include "BivariatePolynomial2.h"

using namespace std;

int main(){
	int i=2;
	char a[]="2.73333333+2.25*x+0.50*y-4.8333333*x^2+0.100*y^2";
	BivariatePolynomial* Bp1 = new BivariatePolynomial(a,i);
	Bp1->Print();
	return 1;
}

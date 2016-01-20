#ifndef BHREJ_H
#define BHREJ_H

#include <Rcpp.h>
#include <algorithm>
#include "priorFuns.h"

using namespace Rcpp;

std::vector<bool> MedWVY (NumericMatrix & , std::vector<int> & ) ;
std::vector<bool> MADWVY (NumericMatrix & , std::vector<int> & ) ; 
std::vector<double> MatRankit (NumericMatrix & ) ; 

class weave {
public:
	weave (NumericMatrix & ) ;
	std::vector<double> rankit () ;
	std::vector<bool> med_y (std::vector<int> & ) ;
	std::vector<bool> MAD_y (std::vector<int> & ) ;
private:
	NumericMatrix mat;
	int nrow = 0;
	int ncol = 0;
	int nweave = 0;
};
#endif

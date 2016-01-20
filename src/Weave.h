#ifndef BHREJ_H
#define BHREJ_H

#include <Rcpp.h>
#include <algorithm>
#include "priorFuns.h"

using namespace Rcpp;

std::vector<bool> MedWVY (const NumericMatrix & , const std::vector<int> & ) ;
std::vector<bool> MADWVY (const NumericMatrix & , const std::vector<int> & ) ; 
std::vector<double> MatRankit (const NumericMatrix & ) ; 

class weave {
public:
	weave (const NumericMatrix & ) ;
	std::vector<double> rankit () ;
	std::vector<bool> med_y (const std::vector<int> & ) ;
	std::vector<bool> MAD_y (const std::vector<int> & ) ;
private:
	NumericMatrix mat;
	int nrow;
	int ncol;
	int nweave;
};
#endif

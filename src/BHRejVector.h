#ifndef BHREJVECTOR_H
#define BHREJVECTOR_H

#include <Rcpp.h>
#include "priorFuns.h"
#include "BHRej.h"

using namespace Rcpp;

std::vector<double> Cpp_fvalue (const std::vector<std::vector<double> > & ) ;
std::vector<double> Cpp_fvalues (const std::vector<NumericMatrix> & );

class BHRejVector {
public:
	// Initialize the object with a list of matrix. 
	BHRejVector (std::vector<NumericMatrix> &);
	// Calculate the p-values based on imput matrix.
	std::vector<double> f_values () ;
	// Return the total columns and rows, and number of group;
	std::vector<int> dim () ;
private:
	std::vector<NumericMatrix> data;
	int n;
	int nrow;
	std::vector<int> ncol;
};

#endif

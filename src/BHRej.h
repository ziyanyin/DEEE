#ifndef BHREJ_H
#define BHREJ_H

#include <algorithm>
#include <Rcpp.h>
#include "priorFuns.h"

using namespace Rcpp;

std::vector<int> Cpp_FDR(std::vector<double> &, const float);

// Analyse a list of p-values.
class BHRej {
public:
	// Make a pair heap, consist of the p-values and their IDs.
	BHRej (std::vector<double> & ) ;
	// Pop out the minimum.
	std::pair<int, double> pop_heap() ;
	// Number of p-values
    int size() ;
	// Given FDR, return the rejecting list.
	std::vector<int> rej_list (float) ;
	// Original data
	std::vector<double> p_values() ;
	double median() ;
private:
    std::vector<std::pair<int, double> > Pdata;
	std::vector<double> pv;
};

// Provide the comparison function need in heap sort in BHRej.
struct comp {
    bool operator() (std::pair<int, double> & a, std::pair<int, double> & b) const {
        return a.second > b.second;
    }
};
#endif

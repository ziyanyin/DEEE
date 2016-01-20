#include "BHRej.h"
#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// Declare the functions exported.
double Cpp_abs (const double) ;
double Cpp_sum (const std::vector<double> & ) ;
double Cpp_var (const std::vector<double> & ) ;
bool comppair1(const std::pair<int, double> & , const std::pair<int, double> & ) ;
bool comppair2(const std::pair<int, double> & , const std::pair<int, double> & ) ;

// Exported for internal usage.
std::vector<std::vector<float> > Cpp_plotMat (const int , const int , const float , const float ) ;
double Cpp_median(const std::vector<double> ) ;
double Cpp_MAD(const std::vector<double> ) ;

// Exported for external usage.
std::vector<int> Cpp_order (const std::vector<double> & ) ;
std::vector<int> Cpp_rank (const std::vector<double> & ) ;

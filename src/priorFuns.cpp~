#include "priorFuns.h"
//// Function Cpp_sum, Cpp_var, Cpp_abs.
double Cpp_sum (const std::vector<double> & l) {
    double res = 0;
    for (auto i = l.begin(); i != l.end(); i++) {
        res += *i;
    }
    return res;
}
double Cpp_var (const std::vector<double> & l) {
    double res = 0, m = Cpp_sum(l) / l.size();
    for (auto i = l.begin(); i != l.end(); i++) {
        res += ((*i) - m) * ((*i) - m);
    }
    return res / (l.size() - 1);
}
double Cpp_abs (const double t) {
	return (t > 0) ? t : (0 - t);
}

//// Functions exported for internal use.
// [[Rcpp::export]]
std::vector<std::vector<float> > Cpp_plotMat (const int nfig, const int figcol, const float xshrink = 1.0, const float yshrink = 1.0) {
    int nfigrow = nfig / figcol + (nfig % figcol != 0);
    int nfigcol = nfig / nfigrow + (nfig % nfigrow != 0);
    std::vector<std::vector<float> > res;

    for (int i = 0; i != nfig; i++) {
        int colNum = i / nfigrow + 1;
        int rowNum = ((i + 1) % nfigrow == 0) ? nfigrow : (i + 1) % nfigrow;
        std::vector<float> tmploc {(colNum - 1) / (float) nfigcol * xshrink, colNum / (float) nfigcol * xshrink, (nfigrow - rowNum) / (float) nfigrow * yshrink, (nfigrow - rowNum + 1) / (float) nfigrow * yshrink};
        res.push_back(tmploc);
    }

    return res;
}

// [[Rcpp::export]]
double Cpp_median(std::vector<double> l) {
	int len = l.size();
	std::vector<double> resvec;
	std::make_heap(l.begin(), l.end());
	for (int i = 0; i != len / 2 + 1; i++) {
		std::pop_heap(l.begin(), l.end());
		resvec.push_back(l.back());
		l.pop_back();
	}
	return (resvec[len / 2] + resvec[(len - 1) / 2]) / 2;
}

// [[Rcpp::export]]
double Cpp_MAD(std::vector<double> l) {
	double m1 = Cpp_median(l);
	std::vector<double> ml;
	for (auto i = l.begin(); i != l.end(); i++) {
		ml.push_back(Cpp_abs(*i - m1));
	}
	return Cpp_median(ml) * 1.4826;
}

//// Functions exported for external use.
//' The order of a vector
//'
//' \code{Cpp_order} returns the order of a given vector. It has the same purpose as 
//' \code{\link[base]{order}}, but is implemented with Cpp and much faster.
//'
//' @param vec A vector of numerics.
//' @return A vector of integers containing the order the given vector is returned.
//' @export
// [[Rcpp::export]]
std::vector<int> Cpp_order (std::vector<double> & vec) {
	int len = vec.size();
	std::vector<int> res;
	std::vector<std::pair<int, double> > mv;
	for (int ir = 0; ir != len; ir++) {
		mv.push_back(std::pair<int, double>(ir, vec[ir]));
	}
	std::sort(mv.begin(), mv.end(), [] (const std::pair<int, double> & a, const std::pair<int, double> & b) {return a.second < b.second; });
	for (auto i = mv.begin(); i != mv.end(); i++) {
		res.push_back((*i).first);
	}
	std::sort(mv.begin(), mv.end(), [] (const std::pair<int, double> & a, const std::pair<int, double> & b) {return a.first < b.first; });
	return res;
}

//' The rank of a vector
//'
//' \code{Cpp_rank} returns the order of a given vector. It has the same purpose as
//' \code{\link[base]{rank}}, but is implemented with Cpp and much faster.
//'
//' @param vec A vector of numerics.
//' @return A vector of integers containing the rank the given vector is returned.
//' @export
// [[Rcpp::export]]
std::vector<int> Cpp_rank (std::vector<double> & vec) {
	int len = vec.size();
	std::vector<int> res(len);
	std::vector<int> to = Cpp_order(vec);
	for(int i = 0; i != len; i++) {
		res[to[i]] = i;
	}
	return res;
}

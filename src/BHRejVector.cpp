#include "BHRejVector.h"

BHRejVector::BHRejVector (std::vector<NumericMatrix> & lm) {
	data = lm;
	n = lm.size();
	nrow = (lm[0]).nrow();
	for (std::vector<NumericMatrix>::iterator i = data.begin(); i != data.end(); i++) {
		ncol.push_back((*i).ncol());
	}
}

std::vector<double> BHRejVector::f_values () {
	std::vector<double> fvalues;
	for (int ir = 0; ir != nrow; ir++) {
		std::vector<std::vector<double> > subdata;
		for (int ik = 0; ik != n; ik++) {
			NumericMatrix subm = data[ik];
			std::vector<double> sgcol;
				for(int j = 0; j != ncol[ik]; j++) {
					sgcol.push_back(subm(ir, j));
				}
			subdata.push_back(sgcol);
		}
		fvalues.push_back(Cpp_fvalue(subdata)[0]);
	}
	return fvalues;
}

std::vector<int> BHRejVector::dim () {
    std::vector<int> res;
    int totaln = 0;
    for (std::vector<int>::iterator i = ncol.begin(); i != ncol.end(); i++) {
        totaln += (*i);
    }
    res.push_back(ncol.size());
    res.push_back(totaln);
    res.push_back(nrow);
    return res;
}

// Produce the f value. Return p-value, samples total sizes and group number.
// [[Rcpp::export]]
std::vector<double> Cpp_fvalue (const std::vector<std::vector<double> > & myData) {
    // Total numbers of groups and the number of treatments.
    int k = myData.size(), n = 0;

    // Theoretical details could be found in Welleck's book.
    double totalsum = 0.0, denumvar = 0.0;
    for (std::vector<std::vector<double> >::const_iterator i = myData.begin(); i != myData.end(); i++) {
        n += (*i).size();
        totalsum += Cpp_sum(*i);
        denumvar += Cpp_var(*i) * ((*i).size() - 1);
    }
    double xbar = totalsum / n, denum = 1 / (double) (n - k) * denumvar, num = 0.0;
    for (std::vector<std::vector<double> >::const_iterator i = myData.begin(); i != myData.end(); i++) {
        num += (*i).size() / (double) n * k * (Cpp_sum(*i) / (*i).size() - xbar) * (Cpp_sum(*i) / (*i).size() - xbar);
    }
	std::vector<double> res;
	res.push_back((denum == 0) ? 0 : num / denum * n / k / (k - 1));
	res.push_back(n);
	res.push_back(k);
    return res;
}

// [[Rcpp::export]]
std::vector<double> Cpp_fvalues (std::vector<NumericMatrix> & lm) {
    BHRejVector tmp = BHRejVector(lm);
    std::vector<double> res;
	std::vector<int> tmp2 = tmp.dim();
    for (std::vector<int>::iterator i = tmp2.begin(); i != tmp2.end(); i++) {
        res.push_back(*i);
    }
    std::vector<double> fvalues = tmp.f_values();
    res.insert(res.end(), fvalues.begin(), fvalues.end());
    return res;
}

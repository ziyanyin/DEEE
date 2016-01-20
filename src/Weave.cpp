#include "Weave.h"

weave::weave (const NumericMatrix & m) {
	mat = m;
	nrow = mat.nrow();
	ncol = mat.ncol();
	nweave = ncol;
}

std::vector<double> weave::rankit () {
	std::vector<double> mv;
	for (int ir = 0; ir != nrow; ir++) {
		double rsum = 0;
		for (int ic = 0; ic != ncol; ic++) {
			rsum += mat(ir, ic);
		}
		mv.push_back(rsum);
	}
	std::vector<int> mr = Cpp_rank(mv);
	std::vector<double> res;
	for (std::vector<int>::iterator i = mr.begin(); i != mr.end(); i++) {
		res.push_back(((*i) + 1.0) / (nrow + 1));
	}
	return res;
}

std::vector<bool> weave::med_y (const std::vector<int> & refl) {
	std::vector<int> ref = refl;
	std::vector<bool> res;
	std::sort(ref.begin(), ref.end());
	std::vector<int> total;
	for (int i = 0; i != ncol; i++) {
		total.push_back(i);
	}
	for (int ir = 0; ir != nrow; ir++) {
		std::vector<double> ref_x;
		std::vector<double> ref_y;
		std::vector<int>::iterator pref = ref.begin();
		std::vector<int>::iterator ptotal = total.begin();
		while (pref != ref.end() && ptotal != total.end()) {
			if ((*ptotal) == (*pref)) {
				ref_x.push_back(mat(ir, *pref));
				pref++;
				ptotal++;
			} else {
				ref_y.push_back(mat(ir, *ptotal));
				ptotal++;
			}
		}
		while (ptotal != total.end()) {
			ref_y.push_back(mat(ir, *ptotal));
			ptotal++;
		}
	res.push_back(Cpp_median(ref_x) > Cpp_median(ref_y));
	}
	return res;
}

std::vector<bool> weave::MAD_y (const std::vector<int> & refl) {
	std::vector<int> ref = refl;
	std::vector<bool> res;
	std::sort(ref.begin(), ref.end());
	std::vector<int> total;
	for (int i = 0; i != ncol; i++) {
		total.push_back(i);
	}
	for (int ir = 0; ir != nrow; ir++) {
		std::vector<double> ref_x;
		std::vector<double> ref_y;
		std::vector<int>::iterator pref = ref.begin();
		std::vector<int>::iterator ptotal = total.begin();
		while (pref != ref.end() && ptotal != total.end()) {
			if ((*ptotal) == (*pref)) {
				ref_x.push_back(mat(ir, *pref));
				pref++;
				ptotal++;
			} else {
				ref_y.push_back(mat(ir, *ptotal));
				ptotal++;
			}
		}
		while (ptotal != total.end()) {
			ref_y.push_back(mat(ir, *ptotal));
			ptotal++;
		}
	res.push_back(Cpp_MAD(ref_x) > Cpp_MAD(ref_y));
	}
	return res;
}

// [[Rcpp::export]]
std::vector<bool> MedWVY (const NumericMatrix & mat, const std::vector<int> & ref) {
	weave wm = weave(mat);
	return wm.med_y(ref);
}

// [[Rcpp::export]]
std::vector<bool> MADWVY (const NumericMatrix & mat, const std::vector<int> & ref) {
	weave wm = weave(mat);
	return wm.MAD_y(ref);
}

// [[Rcpp::export]]
std::vector<double> MatRankit (const NumericMatrix & mat) {
	weave wm = weave(mat);
	return wm.rankit();
}









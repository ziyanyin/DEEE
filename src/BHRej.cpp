#include "BHRej.h"

BHRej::BHRej (std::vector<double> & d) {
    pv = d;
    for (int i = 0; i != d.size(); i++) {
        Pdata.push_back(std::make_pair(i, d[i]));
        std::push_heap(Pdata.begin(), Pdata.end(), comp());
    }
}

std::pair<int, double> BHRej::pop_heap() {
    std::pop_heap(Pdata.begin(), Pdata.end(), comp());
    std::pair<int, double> res = Pdata.back();
    Pdata.pop_back();
    return res;
}

int BHRej::size() {
    return Pdata.size();
}

std::vector<double> BHRej::p_values() {
    return (*this).pv;
}

std::vector<int> BHRej::rej_list (const float FDR) {
    std::vector<int> rl;
    int tn = (*this).size(), i = 0;
    while (i++ < tn) {
        std::pair<int, double> tmpmin = (*this).pop_heap();
        if (tmpmin.second > i * FDR / tn) {
            return rl;
        } else {
            rl.push_back(tmpmin.first);
        }
    }
    return rl;
}

// [[Rcpp::export]]
std::vector<int> Cpp_FDR(std::vector<double> & p, const float FDR = 0.1) {
    std::vector<int> res = (BHRej(p)).rej_list(FDR);
    for(auto i = res.begin(); i != res.end(); i++) {
        (*i)++;
    }
    return res;
}


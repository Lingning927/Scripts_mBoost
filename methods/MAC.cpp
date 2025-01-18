#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix compute_distances(const NumericMatrix& X, const NumericVector& id) {
    int n = X.nrow(), k = X.ncol(), m = id.size();
    NumericMatrix dists(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            double sum = 0;
            for (int l = 0; l < k; ++l) {
                sum += pow(X(i, l) - X(id[j] - 1, l), 2);
            }
            dists(i, j) = sqrt(sum);
        }
    }
    return dists;
}

double Sij(const NumericMatrix& DW, const NumericMatrix& DV, const NumericVector& list_i, const NumericVector& list_j, int n, int i, int j) {
    int k = 0;
    int m = DW.nrow() - n;
    bool PIW, PIV;
    double tmp = 0, R = 0;
    double P[4] = {0};
    double Q[4] = {0};
    for (k = 0; k < n; ++k) {
        PIW = (DW(k, i) <= DW(list_j[j]-1, i));
        PIV = (DV(k, j) <= DV(list_i[i]-1, j));

        if (PIW & PIV) P[0]++;
        if (PIV) P[1]++;
        if (PIW) P[2]++;
    }
    P[1] -= P[0];
    P[2] -= P[0];
    P[3] = n - P[0] - P[1] - P[2];
    for (k = n; k < (n + m); ++k) {
        PIW = (DW(k, i) <= DW(list_j[j]-1, i));
        PIV = (DV(k, j) <= DV(list_i[i]-1, j));

        if (PIW & PIV) Q[0]++;
        if (PIV) Q[1]++;
        if (PIW) Q[2]++;
    }
    Q[1] -= Q[0];
    Q[2] -= Q[0];
    Q[3] = m - Q[0] - Q[1] - Q[2];
    tmp = 0;
    for(k = 0; k < 4; ++k){
        R = 1.0*(P[k] + Q[k])/(n + m);
        if(R > 0.00000001) {
            tmp += (P[k] - n*R) * (P[k] - n*R) / (n*R);
            tmp += (Q[k] - m*R) * (Q[k] - m*R) / (m*R);
        }
    }
    return (tmp);
}

// [[Rcpp::export]]
double compute_MAC(const NumericMatrix& X_bind, const NumericMatrix& Y_bind, const NumericVector& list_i, const NumericVector& list_j, int n) {
    int k1 = list_i.size(), k2 = list_j.size();
    int m = X_bind.nrow() - n;
    double mac = 0, tmp = 0;
    NumericMatrix DW((n + m), k1);
    DW = compute_distances(X_bind, list_i);
    NumericMatrix DV((n + m), k2);
    DV = compute_distances(Y_bind, list_j);
    for(int i = 0; i < k1; ++i) {
        for(int j = 0; j < k2; ++j) {
            tmp = Sij(DW, DV, list_i, list_j, n, i, j);
            mac = (mac > tmp) ? mac : tmp;
        }
    }
    return(mac);
}
#include <Rcpp.h>  
using namespace Rcpp;


// [[Rcpp::export(rng = false)]]
double sscXYopt(NumericMatrix X, NumericMatrix Y, int n1, int n2, int py, int px){
    double summ = 0; 
    int i, ii, iii;
    
    for(i = 0; i < px; i++){
        double sumc = 0;
        for(ii = 0; ii < py; ii++){
            double sumv = 0;
            for(iii = 0; iii < n1; iii+=4){
                sumv += (Y(iii, ii) * X(iii, i)) + (Y(iii + 1, ii) * X(iii + 1, i)) + 
                    (Y(iii + 2, ii) * X(iii + 2, i)) + (Y(iii + 3, ii) * X(iii + 3, i));
            }
            if (n2 == 1) {
                sumv += Y(n1, ii) * X(n1, i);
            } else if (n2 == 2) {
                sumv += Y(n1, ii) * X(n1, i) + Y(n1+1, ii) * X(n1+1, i);
            } else if (n2 == 3) {
                sumv += Y(n1, ii) * X(n1, i) + Y(n1+1, ii) * X(n1+1, i) + Y(n1+2, ii) * X(n1+2, i);
            }
            sumc += sumv * sumv;
        }
        summ += sumc;
    }
    return summ; 
}

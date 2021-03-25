// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;



arma::mat shuffle(arma::mat Y, arma::uvec s){
  s = s - 1;
  return Y.rows(s);
}


Rcpp::List rrpp(Rcpp::List Fitted, Rcpp::List Residuals, arma::uvec s, int k){
  Rcpp::List out(k);
  int i;
  for(i = 0; i < k; i++) {
    arma::mat R = shuffle(Residuals[i], s);
    arma::mat F = Fitted[i];
    arma::mat Y = F + R;
    out[i] = Y;
  }
  return(out);
}


double vsum(arma::vec v){
  int n = v.n_elem;
  int n2 = n % 5;
  int n1 = n - n2;
  double sum = 0;
  int i;
  
  for(i=0; i < n1; i+=5){
    sum += (v[i] + v[i+1] + v[i+2] + v[i+3] + v[i+4]);
  }
  
  if(n2 == 1) {
    sum += v[n1];
  } else if(n2 == 2){
    sum += (v[n1] + v[n1+1]);
  } else if(n2 == 3){
    sum += (v[n1] + v[n1+1] + v[n1+2]);
  } else if(n2 == 4){
    sum += (v[n1] + v[n1+1] + v[n1+2] + v[n1+3]);
  } 
  
  return(sum);
}


// [[Rcpp::export]]
double sscpUY(arma::mat U, arma::mat Y){
  int i, ii;
  double sum = 0;
  int py = Y.n_cols;
  int pu = U.n_cols;
  int p2 = py % 4;
  int p1 = py - p2;
  
  for(i = 0; i < pu; i++) {
    arma::vec u = U.col(i);
    
    for(ii = 0; ii < p1; ii+=4){
      double v = vsum(u % Y.col(ii));
      double v1 = vsum(u % Y.col(ii+1));
      double v2 = vsum(u % Y.col(ii+2));
      double v3 = vsum(u % Y.col(ii+3));
      sum += (v * v + v1 * v1 + v2 * v2 + v3 *v3);
    }
      
    if(p2 == 1) {
      double v = vsum(u % Y.col(p1));
      sum += v * v;
    } else if(p2 == 2){
      double v = vsum(u % Y.col(p1));
      double v1 = vsum(u % Y.col(p1+1));
      sum += (v * v + v1 * v1);
    } 
    else if(p2 == 3){
      double v = vsum(u % Y.col(p1));
      double v1 = vsum(u % Y.col(p1+1));
      double v2 = vsum(u % Y.col(p1+2));
      sum += (v * v + v1 * v1 + v2 * v2); 
      }
    }
  return(sum);
}



Rcpp::List listSScpUY(Rcpp::List U, Rcpp::List Y, 
                      int k){
  int i;
  Rcpp::List out(k);
  for(i = 0; i < k; i++) {
    arma::mat u = U[i];
    arma::mat y = Y[i];
    double ss = sscpUY(u, y);
    out[i] = ss;
  }
  return(out);
}


double ssY(arma::mat Y){
  return(arma::accu(Y % Y));
}


Rcpp::List listRSS(arma::mat Uf, 
                   Rcpp::List Y,
                   int k){
  int i;
  Rcpp::List out(k);
  
  for(i = 0; i < k; i++) {
    arma::mat y = Y[i];
    out[i] = ssY(y) - sscpUY(Uf, y);
  }
  return(out);
}


Rcpp::List getStats(Rcpp::List Ur, Rcpp::List Uf, 
                    arma::mat Ufull, arma::mat Unull,
                    Rcpp::List Fitted, Rcpp::List Residuals,
                    arma::uvec s, arma::mat Y0,
                    int k){
  
  Rcpp::List out(k + 2);
  Rcpp::List Ylist = rrpp(Fitted, Residuals, s, k);
  Rcpp::List red = listSScpUY(Ur, Ylist, k);
  Rcpp::List full = listSScpUY(Uf, Ylist, k);
  Rcpp::List rss = listRSS(Ufull, Ylist, k);
  double yy = ssY(Y0);
  double tss = yy - sscpUY(Unull, Y0);
  double rssM = yy - sscpUY(Ufull, Y0);
  
  return Rcpp::List::create(Rcpp::Named("red") = red,
                            Rcpp::Named("full") = full,
                            Rcpp::Named("rss") = rss,
                            Rcpp::Named("tss") = tss,
                            Rcpp::Named("rssM") = rssM);
}


// [[Rcpp::export]]
Rcpp::List iterSS(Rcpp::List ind, Rcpp::List Ur, Rcpp::List Uf, 
                  arma::mat Ufull, arma::mat Unull,
                  Rcpp::List Fitted, Rcpp::List Residuals,
                  arma::mat Yh0, arma::mat R0, int k){
  
  int perms = ind.size();
  Rcpp::List out(perms);
  arma::mat Y0 = Yh0;
  int i;
  
  for(i = 0; i < perms; i++){
    arma::uvec s = ind[i];
    Y0 = Yh0 + shuffle(R0, s);
    out[i] = getStats(Ur, Uf, Ufull, Unull, 
                      Fitted, Residuals, s, Y0, k);
  }
  return(out); 
}

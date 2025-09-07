#include <Rcpp.h>
using namespace Rcpp;

// helper function for max of a set of ints
R_xlen_t max_len(IntegerVector vars){
  int n = vars.size();
  int top_size;

  if(n==0){
    top_size=0;
  }else{
    top_size = vars[0];
    for (int i=1; i<n; i++){
      if(vars[i] > top_size){
        top_size = vars[i];
      }
    }
  }
  return top_size;
}


///// ERLANG DIST /////

// PDF
double derlang_raw(double x, double k, double l, bool log=false) {
  double f;

  if(x < 0){
    f = log ? R_NegInf : 0.0;
  }else{
    double density = std::pow(l,k) * std::pow(x, k-1) * std::exp(-l * x) / std::tgamma(k);
    f = log ? std::log(density) : density;
  }
  return f;
}

// CDF
double perlang_raw(double x, int k, double l, bool lower_tail=true, bool log_p=false) {
  double F;

  if(x < 0){
    F = log_p ? R_NegInf : 0.0;
  }else{

    double lambda_x = l*x;
    double s = 0;
    for(int j=0; j<k; j++){
      double term = std::pow(l*x, j) / std::tgamma(j+1.0) ;
      s += term;
    }
    double p = 1.0 - std::exp(-lambda_x) * s;
    if(p<1e-15)p=1e-15;   //clamp??

    F = lower_tail ? p : 1-p;
    F = log_p ? std::log(F): F;

  }
  return F;
}

// [[Rcpp::export]]
NumericVector derlang_c(
  NumericVector x,
  NumericVector k,
  NumericVector l,
  bool log=false
){
  int xs = x.size();
  int ks = k.size();
  int ls = l.size();
  const R_xlen_t n = max_len({xs,ks,ls});

  NumericVector x2 = Rcpp::rep_len(x, n);
  NumericVector k2 = Rcpp::rep_len(k, n);
  NumericVector l2 = Rcpp::rep_len(l, n);
  NumericVector f(n);

  for(int i=0; i<n; i++){
    f[i] = derlang_raw(x2[i], k2[i], l2[i], log);
  }

  return f;
}

// [[Rcpp::export]]
NumericVector perlang_c(
    NumericVector x,
    NumericVector k,
    NumericVector l,
    bool lower_tail=true,
    bool log_p=false
){
  int xs = x.size();
  int ks = k.size();
  int ls = l.size();
  const R_xlen_t n = max_len({xs,ks,ls});

  NumericVector x2 = Rcpp::rep_len(x, n);
  NumericVector k2 = Rcpp::rep_len(k, n);
  NumericVector l2 = Rcpp::rep_len(l, n);
  NumericVector F(n);

  for(int i=0; i<n; i++){
    F[i] = perlang_raw(x2[i], k2[i], l2[i], lower_tail, log_p);
  }

  return F;
}

///// GAMMA GOMPERTZ //////

// PDF
double dgamgomp_raw(double x, double b, double s, double beta, bool log=false) {
  double f;
  if(x <0){
    f = log ? R_NegInf : 0.0;
  }else{

    double bx = b*x;
    double top = b * s * std::exp(bx) * pow(beta,s);
    double bot1 = beta - 1 + exp(bx);
    double bot = pow(bot1,s+1);
    double frac = top/bot;
    f = log ? std::log(frac) : frac;
  }
  return f;
}

// CDF
double pgamgomp_raw(double x, double b, double s, double beta, bool lower_tail=true, bool log_p=false) {
  double F;
  if(x <0){
    F = log_p ? R_NegInf : 0.0;

  }else{
    double bx = b*x;
    double bot = beta-1+std::exp(bx);
    double test = 1 - (pow(beta,s) / pow(bot,s));

    F = lower_tail ? test : 1-test;
    F = log_p ? std::log(F): F;
  }
  return F;
}

// [[Rcpp::export]]
NumericVector dgamgomp_c(
      NumericVector x,
      NumericVector b,
      NumericVector s,
      NumericVector beta,
      bool log=false
  ){
    int xs = x.size();
    int bs = b.size();
    int ss = s.size();
    int betas = beta.size();
    const R_xlen_t n = max_len({xs,bs,ss,betas});

    NumericVector x2 =    Rcpp::rep_len(x, n);
    NumericVector b2 =    Rcpp::rep_len(b, n);
    NumericVector s2 =    Rcpp::rep_len(s, n);
    NumericVector beta2 = Rcpp::rep_len(beta, n);

    NumericVector f(n);

    for(int i=0; i<n; i++){
      f[i] = dgamgomp_raw(x2[i], b2[i], s2[i], beta2[i], log);
    }

    return f;
  }

// [[Rcpp::export]]
NumericVector pgamgomp_c(
    NumericVector x,
    NumericVector b,
    NumericVector s,
    NumericVector beta,
    bool lower_tail=true,
    bool log_p=false
){
  int xs = x.size();
  int bs = b.size();
  int ss = s.size();
  int betas = beta.size();
  const R_xlen_t n = max_len({xs,bs,ss,betas});

  NumericVector x2 =    Rcpp::rep_len(x, n);
  NumericVector b2 =    Rcpp::rep_len(b, n);
  NumericVector s2 =    Rcpp::rep_len(s, n);
  NumericVector beta2 = Rcpp::rep_len(beta, n);

  NumericVector F(n);

  for(int i=0; i<n; i++){
    F[i] = pgamgomp_raw(x2[i], b2[i], s2[i], beta2[i], lower_tail, log_p);
  }

  return F;
}


///// LOG CAUCHY //////

// PDF
double dlogcauchy_raw(double x, double mu, double sigma, bool log=false) {
  double f;
  if(x <=0){
    f = log ? R_NegInf : 0.0;
  }else{
  const double pi = 3.14159265358979323846;
  double term1 = sigma / (x * pi);
  double term2 = std::log(x) - mu;
  double term3 = std::pow(term2,2) + std::pow(sigma,2);
  double term = term1 / term3;
    f = log ? std::log(term) : term;
  }
  return f;
}

// CDF
double plogcauchy_raw(double x, double mu, double sigma, bool lower_tail=true, bool log_p=false) {
  double F;
  if(x <=0){
    F = log_p ? R_NegInf : 0.0;
  }else{
    const double pi = 3.14159265358979323846;
    double norm = (std::log(x) - mu)/sigma;
    double term = .5 + (std::atan(norm) / pi);
    F = lower_tail ? term : 1-term;
    F = log_p ? std::log(F): F;
  }
  return F;
}

// [[Rcpp::export]]
NumericVector dlogcauchy_c(
    NumericVector x,
    NumericVector mu,
    NumericVector sigma,
    bool log=false
){
  int xs = x.size();
  int ms = mu.size();
  int ss = sigma.size();
  const R_xlen_t n = max_len({xs,ms,ss});

  NumericVector x2 =    Rcpp::rep_len(x, n);
  NumericVector m2 =    Rcpp::rep_len(mu, n);
  NumericVector s2 =    Rcpp::rep_len(sigma, n);

  NumericVector f(n);

  for(int i=0; i<n; i++){
    f[i] = dlogcauchy_raw(x2[i], m2[i], s2[i], log);
  }

  return f;
}

// [[Rcpp::export]]
NumericVector plogcauchy_c(
    NumericVector x,
    NumericVector mu,
    NumericVector sigma,
    bool lower_tail=true,
    bool log_p=false
){
  int xs = x.size();
  int ms = mu.size();
  int ss = sigma.size();
  const R_xlen_t n = max_len({xs,ms,ss});

  NumericVector x2 =    Rcpp::rep_len(x, n);
  NumericVector m2 =    Rcpp::rep_len(mu, n);
  NumericVector s2 =    Rcpp::rep_len(sigma, n);

  NumericVector F(n);

  for(int i=0; i<n; i++){
    F[i] = plogcauchy_raw(x2[i], m2[i], s2[i], lower_tail, log_p);
  }

  return F;
}


/////// Hypertabastic Dist

// helper functions
double sech(double x) {
  return 1.0 / std::cosh(x);
}

double coth(double x) {
  return std::cosh(x) / std::sinh(x);
}

double csch(double x) {
  return 1.0 / std::sinh(x);
}

// PDF
double dhypertab_raw(double x, double a, double b, bool log=false) {
  double f;
  if(x <=0){
    f = log ? R_NegInf : 0.0;
  }else{
    double wt = a*(1 - pow(x, b) * coth(pow(x, b)))/b;
    double chunk_a = a*pow(x,(2*b)-1)*pow(csch(pow(x, b)),2) ;
    double chunk_b = a*pow(x, b-1) * coth(pow(x,b));
    double chunk_c = sech(wt) * tanh(wt) * (chunk_a - chunk_b);
    f = log ? std::log(chunk_c) : chunk_c;
  }
  return f;
}

// CDF
double phypertab_raw(double x, double a, double b, bool lower_tail=true, bool log_p=false) {
  double F;
  if(x <=0){
    F = log_p ? R_NegInf : 0.0;
  }else{
    double xb = pow(x,b);
    double term1 = a * (1 - xb * coth(xb)) / b;
    double term2 = 1 - sech(term1);
    F = lower_tail ? term2 : 1-term2;
    F = log_p ? std::log(F): F;
  }
  return F;
}

// [[Rcpp::export]]
NumericVector dhypertab_c(
    NumericVector x,
    NumericVector a,
    NumericVector b,
    bool log=false
){
  int xs = x.size();
  int as = a.size();
  int bs = b.size();
  const R_xlen_t n = max_len({xs,as,bs});

  NumericVector x2 =    Rcpp::rep_len(x, n);
  NumericVector a2 =    Rcpp::rep_len(a, n);
  NumericVector b2 =    Rcpp::rep_len(b, n);

  NumericVector f(n);

  for(int i=0; i<n; i++){
    f[i] = dhypertab_raw(x2[i], a2[i], b2[i], log);
  }

  return f;
}

// [[Rcpp::export]]
NumericVector phypertab_c(
    NumericVector x,
    NumericVector a,
    NumericVector b,
    bool lower_tail=true,
    bool log_p=false
){
  int xs = x.size();
  int as = a.size();
  int bs = b.size();
  const R_xlen_t n = max_len({xs,as,bs});

  NumericVector x2 =    Rcpp::rep_len(x, n);
  NumericVector a2 =    Rcpp::rep_len(a, n);
  NumericVector b2 =    Rcpp::rep_len(b, n);

  NumericVector F(n);

  for(int i=0; i<n; i++){
    F[i] = phypertab_raw(x2[i], a2[i], b2[i], lower_tail, log_p);
  }

  return F;
}

///// INVERSE LINDLEY //////

// PDF
double dinvlind_raw(double x, double theta, bool log=false) {
  double f;
  if(x <=0){
    f = log ? R_NegInf : 0.0;
  }else{
    double term = pow(theta,2)/(1+theta) * ((1+x)/pow(x,3)) * std::exp(-theta/x);
    f = log ? std::log(term) : term;
  }
  return f;
}

// CDF
double pinvlind_raw(double x, double theta, bool lower_tail=true, bool log_p=false) {
  double F;
  if(x <=0){
    F = log_p ? R_NegInf : 0.0;
  }else{
    double term = (1 + ((theta/(1+theta))/x)) * std::exp(-theta/x);
    F = lower_tail ? term : 1-term;
    F = log_p ? std::log(F): F;
  }
  return F;
}

// [[Rcpp::export]]
NumericVector dinvlind_c(
    NumericVector x,
    NumericVector theta,
    bool log=false
){
  int xs = x.size();
  int ts = theta.size();
  const R_xlen_t n = max_len({xs,ts});

  NumericVector x2 =    Rcpp::rep_len(x, n);
  NumericVector t2 =    Rcpp::rep_len(theta, n);
  NumericVector f(n);

  for(int i=0; i<n; i++){
    f[i] = dinvlind_raw(x2[i], t2[i],  log);
  }

  return f;
}

// [[Rcpp::export]]
NumericVector pinvlind_c(
    NumericVector x,
    NumericVector theta,
    bool lower_tail=true,
    bool log_p=false
){
  int xs = x.size();
  int ts = theta.size();
  const R_xlen_t n = max_len({xs,ts});

  NumericVector x2 =    Rcpp::rep_len(x, n);
  NumericVector t2 =    Rcpp::rep_len(theta, n);
  NumericVector F(n);

  for(int i=0; i<n; i++){
    F[i] = pinvlind_raw(x2[i], t2[i], lower_tail, log_p);
  }

  return F;
}

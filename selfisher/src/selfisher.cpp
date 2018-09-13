#define TMB_LIB_INIT R_init_selfisher
#include <TMB.hpp>

namespace selfisher{

  template<class Type>
  bool isNA(Type x){
    return R_IsNA(asDouble(x));
  }

  extern "C" {
    /* See 'R-API: entry points to C-code' (Writing R-extensions) */
    double Rf_logspace_sub (double logx, double logy);
    void   Rf_pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
  }

  /* y(x) = logit_invcloglog(x) := log( exp(exp(x)) - 1 ) = logspace_sub( exp(x), 0 )

     y'(x) = exp(x) + exp(x-y) = exp( logspace_add(x, x-y) )

   */
  TMB_ATOMIC_VECTOR_FUNCTION(
                             // ATOMIC_NAME
                             logit_invcloglog
                             ,
                             // OUTPUT_DIM
                             1,
                             // ATOMIC_DOUBLE
                             ty[0] = Rf_logspace_sub(exp(tx[0]), 0.);
                             ,
                             // ATOMIC_REVERSE
                             px[0] = exp( logspace_add(tx[0], tx[0]-ty[0]) ) * py[0];
                             )
  template<class Type>
  Type logit_invcloglog(Type x) {
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return logit_invcloglog(tx)[0];
  }

  /* y(x) = logit_pnorm(x) := logit( pnorm(x) ) =
     pnorm(x, lower.tail=TRUE,  log.p=TRUE) -
     pnorm(x, lower.tail=FALSE, log.p=TRUE)

     y'(x) = dnorm(x) * ( (1+exp(y)) + (1+exp(-y)) )

  */
  double logit_pnorm(double x) {
    double log_p_lower, log_p_upper;
    Rf_pnorm_both(x, &log_p_lower, &log_p_upper, 2 /* both tails */, 1 /* log_p */);
    return log_p_lower - log_p_upper;
  }
  TMB_ATOMIC_VECTOR_FUNCTION(
                             // ATOMIC_NAME
                             logit_pnorm
                             ,
                             // OUTPUT_DIM
                             1,
                             // ATOMIC_DOUBLE
                             ty[0] = logit_pnorm(tx[0])
                             ,
                             // ATOMIC_REVERSE
                             Type zero = 0;
                             Type tmp1 = logspace_add(zero, ty[0]);
                             Type tmp2 = logspace_add(zero, -ty[0]);
                             Type tmp3 = logspace_add(tmp1, tmp2);
                             Type tmp4 = dnorm(tx[0], Type(0), Type(1), true) + tmp3;
                             px[0] = exp( tmp4 ) * py[0];
                             )
  template<class Type>
  Type logit_pnorm(Type x) {
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return logit_pnorm(tx)[0];
  }

}

enum valid_family {
  binomial_family = 100
};

enum valid_link {
  logit_link             = 0,
  probit_link            = 1,
  cloglog_link           = 2,
  loglog_link            = 3,
  richards_link          = 4
};

enum valid_covStruct {
  diag_covstruct = 0,
  us_covstruct   = 1,
  cs_covstruct   = 2,
  ar1_covstruct  = 3,
  ou_covstruct   = 4,
  exp_covstruct = 5,
  gau_covstruct = 6,
  mat_covstruct = 7,
  toep_covstruct = 8
};

enum valid_ppredictCode {
  response_ppredictcode = 0,
  selection_ppredictcode = 1,
  prob_ppredictcode = 2,
  ratio_ppredictcode = 3
};

template<class Type>
Type linkfun(Type p, Type etad, int link) {
  Type ans;
  switch (link) {
  case logit_link:
    ans = logit(p);
    break;
  case probit_link:
    ans = qnorm(p);
    break;
  case cloglog_link:
    ans = log(-log(Type(1)-p));
    break;
  case loglog_link:
    ans = -log(-log(p));
    break;
  case richards_link:
    ans = logit(pow(p, exp(etad)));
    break;
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

template<class Type>
Type inverse_linkfun(Type eta, Type etad, int link) {
  Type ans;
  switch (link) {
  case logit_link:
    ans = invlogit(eta);
    break;
  case probit_link:
    ans = pnorm(eta);
    break;
  case cloglog_link:
    ans = Type(1) - exp(-exp(eta));
    break;
  case loglog_link:
    ans = exp(-exp(-eta));
    break;
  case richards_link:
    ans = pow(exp(eta)/(Type(1)+exp(eta)), Type(1)/exp(etad));
    break;
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}

/* logit transformed inverse_linkfun without losing too much
 accuracy */
template<class Type>
Type logit_inverse_linkfun(Type eta, Type etad, int link) {
  Type ans;
  switch (link) {
  case logit_link:
    ans = eta;
    break;
  case probit_link:
    ans = selfisher::logit_pnorm(eta);
    break;
  case cloglog_link:
    ans = selfisher::logit_invcloglog(eta);
    break;
  case loglog_link:
    ans = -exp(-eta)-logspace_sub(Type(0), -exp(-eta));
    break;
  case richards_link:
    ans = (eta- logspace_add(Type(0), eta))/exp(etad) -
      logspace_sub(Type(0), (eta- logspace_add(Type(0), eta))/exp(etad));
    break;
  default:
    ans = logit( inverse_linkfun(eta, etad, link) );
  } // End switch
  return ans;
}
//template<class Type>
//Type logit_phifun(Type etar, Type etad, Type etap, int link, int cover) {
//  Type logit_phi;
//  if(!cover) { //trowser-trawl
//    Type p=invlogit(etap);
//    if(link==logit_link) {
//      logit_phi = log(p) + etar - log(1-p+exp(etar)-p*exp(etar));
//    } else {
//      Type r = inverse_linkfun(etar, etad, link);
//      logit_phi = log(p) + log(r) - log(Type(1.0)-p);
//    }
//  } else { //cover=1 covered codend
//    Type r = inverse_linkfun(etar, etad, link);
//    logit_phi = log(r) - log(Type(1.0)-r);
//  }
//	return logit_phi;
//}

//template<class Type>
//Type phifun(Type etar, Type etad, Type etap, int link, int cover) {
//  Type phi;
//  if(!cover) { //trowser-trawl
//    Type p=invlogit(etap);
//    if(link==logit_link) {
//      phi = p*exp(etar) /(1-p+exp(etar));
//    } else {
//      Type r = inverse_linkfun(etar, etad, link);
//      phi = invlogit(log(p) + log(r) - log(Type(1.0)-p));
//    }
//  } else { //cover=1 covered codend
//      phi = inverse_linkfun(etar, etad, link);
//  }
//	return phi;
//}

template <class Type>
struct per_term_info {
  // Input from R
  int blockCode;     // Code that defines structure
  int blockSize;     // Size of one block
  int blockReps;     // Repeat block number of times
  int blockNumTheta; // Parameter count per block
  matrix<Type> dist;
  vector<Type> times;// For ar1 case
  // Report output
  matrix<Type> corr;
  vector<Type> sd;
};

template <class Type>
struct terms_t : vector<per_term_info<Type> > {
  terms_t(SEXP x){
    (*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP y = VECTOR_ELT(x, i);    // y = x[[i]]
      int blockCode = (int) REAL(getListElement(y, "blockCode", &isNumericScalar))[0];
      int blockSize = (int) REAL(getListElement(y, "blockSize", &isNumericScalar))[0];
      int blockReps = (int) REAL(getListElement(y, "blockReps", &isNumericScalar))[0];
      int blockNumTheta = (int) REAL(getListElement(y, "blockNumTheta", &isNumericScalar))[0];
      (*this)(i).blockCode = blockCode;
      (*this)(i).blockSize = blockSize;
      (*this)(i).blockReps = blockReps;
      (*this)(i).blockNumTheta = blockNumTheta;
      // Optionally, pass time vector:
      SEXP t = getListElement(y, "times");
      if(!isNull(t)){
	RObjectTestExpectedType(t, &isNumeric, "times");
	(*this)(i).times = asVector<Type>(t);
      }
      // Optionally, pass distance matrix:
      SEXP d = getListElement(y, "dist");
      if(!isNull(d)){
	RObjectTestExpectedType(d, &isMatrix, "dist");
	(*this)(i).dist = asMatrix<Type>(d);
      }
    }
  }
};

template <class Type>
Type termwise_nll(array<Type> &U, vector<Type> theta, per_term_info<Type>& term, bool do_simulate = false) {
  Type ans = 0;
  if (term.blockCode == diag_covstruct){
    // case: diag_covstruct
    vector<Type> sd = exp(theta);
    for(int i = 0; i < term.blockReps; i++){
      ans -= dnorm(vector<Type>(U.col(i)), Type(0), sd, true).sum();
      if (do_simulate) {
        U.col(i) = rnorm(Type(0), sd);
      }
    }
    term.sd = sd; // For report
  }
  else if (term.blockCode == us_covstruct){
    // case: us_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> corr_transf = theta.tail(theta.size() - n);
    vector<Type> sd = exp(logsd);
    density::UNSTRUCTURED_CORR_t<Type> nldens(corr_transf);
    density::VECSCALE_t<density::UNSTRUCTURED_CORR_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == cs_covstruct){
    // case: cs_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    Type corr_transf = theta(n);
    vector<Type> sd = exp(logsd);
    Type a = Type(1) / (Type(n) - Type(1));
    Type rho = invlogit(corr_transf) * (Type(1) + a) - a;
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
	corr(i,j) = (i==j ? Type(1) : rho);
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == toep_covstruct){
    // case: toep_covstruct
    int n = term.blockSize;
    vector<Type> logsd = theta.head(n);
    vector<Type> sd = exp(logsd);
    vector<Type> parms = theta.tail(n-1);              // Corr parms
    parms = parms / sqrt(Type(1.0) + parms * parms );  // Now in (-1,1)
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
        corr(i,j) = (i==j ? Type(1) :
                     parms( (i > j ? i-j : j-i) - 1 ) );
    density::MVNORM_t<Type> nldens(corr);
    density::VECSCALE_t<density::MVNORM_t<Type> > scnldens = density::VECSCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = nldens.cov(); // For report
    term.sd = sd;             // For report
  }
  else if (term.blockCode == ar1_covstruct){
    // case: ar1_covstruct
    //  * NOTE: Valid parameter space is phi in [-1, 1]
    //  * NOTE: 'times' not used as we assume unit distance between consecutive time points.
    int n = term.blockSize;
    Type logsd = theta(0);
    Type corr_transf = theta(1);
    Type phi = corr_transf / sqrt(1.0 + pow(corr_transf, 2));
    Type sd = exp(logsd);
    for(int j = 0; j < term.blockReps; j++){
      ans -= dnorm(U(0, j), Type(0), sd, true);   // Initialize
      if (do_simulate) {
        U(0, j) = rnorm(Type(0), sd);
      }
      for(int i=1; i<n; i++){
        ans -= dnorm(U(i, j), phi * U(i-1, j), sd * sqrt(1 - phi*phi), true);
        if (do_simulate) {
          U(i, j) = rnorm( phi * U(i-1, j), sd * sqrt(1 - phi*phi) );
        }
      }
    }
    // For consistency with output for other structs we report entire
    // covariance matrix.
    if(isDouble<Type>::value) { // Disable AD for this part
      term.corr.resize(n,n);
      term.sd.resize(n);
      for(int i=0; i<n; i++){
	term.sd(i) = sd;
	for(int j=0; j<n; j++){
	  term.corr(i,j) = pow(phi, abs(i-j));
	}
      }
    }
  }
  else if (term.blockCode == ou_covstruct){
    // case: ou_covstruct
    //  * NOTE: this is the continuous time version of ar1.
    //          One-step correlation must be non-negative
    //  * NOTE: 'times' assumed sorted !
    int n = term.times.size();
    Type logsd = theta(0);
    Type corr_transf = theta(1);
    Type sd = exp(logsd);
    for(int j = 0; j < term.blockReps; j++){
      ans -= dnorm(U(0, j), Type(0), sd, true);   // Initialize
      if (do_simulate) {
        U(0, j) = rnorm(Type(0), sd);
      }
      for(int i=1; i<n; i++){
	Type rho = exp(-exp(corr_transf) * (term.times(i) - term.times(i-1)));
	ans -= dnorm(U(i, j), rho * U(i-1, j), sd * sqrt(1 - rho*rho), true);
        if (do_simulate) {
          U(i, j) = rnorm( rho * U(i-1, j), sd * sqrt(1 - rho*rho));
        }
      }
    }
    // For consistency with output for other structs we report entire
    // covariance matrix.
    if(isDouble<Type>::value) { // Disable AD for this part
      term.corr.resize(n,n);
      term.sd.resize(n);
      for(int i=0; i<n; i++){
	term.sd(i) = sd;
	for(int j=0; j<n; j++){
	  term.corr(i,j) =
	    exp(-exp(corr_transf) * CppAD::abs(term.times(i) - term.times(j)));
	}
      }
    }
  }
  // Spatial correlation structures
  else if (term.blockCode == exp_covstruct ||
           term.blockCode == gau_covstruct ||
           term.blockCode == mat_covstruct){
    int n = term.blockSize;
    matrix<Type> dist = term.dist;
    if(! ( dist.cols() == n && dist.rows() == n ) )
      error ("Dimension of distance matrix must equal blocksize.");
    // First parameter is sd
    Type sd = exp( theta(0) );
    // Setup correlation matrix
    matrix<Type> corr(n,n);
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
        switch (term.blockCode) {
        case exp_covstruct:
          corr(i,j) = (i==j ? Type(1) : exp( -dist(i,j) * exp(-theta(1)) ) );
          break;
        case gau_covstruct:
          corr(i,j) = (i==j ? Type(1) : exp( -pow(dist(i,j),2) * exp(-2. * theta(1)) ) );
          break;
        case mat_covstruct:
          corr(i,j) = (i==j ? Type(1) : matern( dist(i,j),
                                                exp(theta(1)) /* range */,
                                                exp(theta(2)) /* smoothness */) );
          break;
        default:
          error("Not implemented");
        }
      }
    }
    density::MVNORM_t<Type> nldens(corr);
    density::SCALE_t<density::MVNORM_t<Type> > scnldens = density::SCALE(nldens, sd);
    for(int i = 0; i < term.blockReps; i++){
      ans += scnldens(U.col(i));
      if (do_simulate) {
        U.col(i) = sd * nldens.simulate();
      }
    }
    term.corr = corr;   // For report
    term.sd.resize(n);  // For report
    term.sd.fill(sd);
  }
  else error("covStruct not implemented!");
  return ans;
}

template <class Type>
Type allterms_nll(vector<Type> &u, vector<Type> theta,
		  vector<per_term_info<Type> >& terms,
                  bool do_simulate = false) {
  Type ans = 0;
  int upointer = 0;
  int tpointer = 0;
  int nr, np = 0, offset;
  for(int i=0; i < terms.size(); i++){
    nr = terms(i).blockSize * terms(i).blockReps;
    // Note: 'blockNumTheta=0' ==> Same parameters as previous term.
    bool emptyTheta = ( terms(i).blockNumTheta == 0 );
    offset = ( emptyTheta ? -np : 0 );
    np     = ( emptyTheta ?  np : terms(i).blockNumTheta );
    vector<int> dim(2);
    dim << terms(i).blockSize, terms(i).blockReps;
    array<Type> useg( &u(upointer), dim);
    vector<Type> tseg = theta.segment(tpointer + offset, np);
    ans += termwise_nll(useg, tseg, terms(i), do_simulate);
    upointer += nr;
    tpointer += terms(i).blockNumTheta;
  }
  return ans;
}

template <class Type>
Type calcLprob(Type etar, Type L, Type b, Type etad, Type prob, int link) {
  Type Lprob = (linkfun(prob, etad, link) - etar + L*b)/b;
  return Lprob;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(Xr);
  DATA_SPARSE_MATRIX(Zr);
  DATA_MATRIX(Xp);
  DATA_SPARSE_MATRIX(Zp);
  DATA_MATRIX(Xd);
  DATA_VECTOR(yobs);
  DATA_VECTOR(total);
  DATA_VECTOR(offset);
//  DATA_VECTOR(retp); //retention probability to predict corresponding L

  // Define covariance structure for the selectivity model
  DATA_STRUCT(termsr, terms_t); //check this

  // Define covariance structure for the power model
  DATA_STRUCT(termsp, terms_t);

  // Parameters related to design matrices
  PARAMETER_VECTOR(betar);
  PARAMETER_VECTOR(betap);
  PARAMETER_VECTOR(betad);
  PARAMETER_VECTOR(br);
  PARAMETER_VECTOR(bp);

  // Joint vector of covariance parameters
  PARAMETER_VECTOR(thetar);
  PARAMETER_VECTOR(thetap);

  DATA_INTEGER(link);

  // Flags
  DATA_INTEGER(pPredictCode);
  DATA_INTEGER(doPredict);
//  DATA_INTEGER(Lindex);
  DATA_INTEGER(Lpflag);
  DATA_INTEGER(cover);
  DATA_IVECTOR(whichPredict);

  // Joint negative log-likelihood
  Type jnll = 0;

  // Random effects
  jnll += allterms_nll(br, thetar, termsr, this->do_simulate);
  jnll += allterms_nll(bp, thetap, termsp, this->do_simulate);

  // Linear predictor
  vector<Type> etar = Xr * betar + Zr * br + offset;
  vector<Type> etap = Xp * betap + Zp * bp;
  vector<Type> etad = Xd * betad;

  // Apply link
  vector<Type> r(etar.size());
  for (int i = 0; i < r.size(); i++)
    r(i) = inverse_linkfun(etar(i), etad(i), link);

  // Calculate binomial probability parameter (phi)
  //vector<Type> phi(r.size());
  vector<Type> logit_phi(r.size());
  vector<Type> p(etap.size());
  if(!cover) { //trowser-trawl, alternate haul, catch comparison
    p = invlogit(etap);
    //phi=p*r/(p*r+Type(1)-p);// as in eqn 3 of Wileman et al. 1996, not like in glmmTMB
    logit_phi = log(p) + log(r) - log(Type(1.0)-p);
  } else { //cover=1 covered codend
    logit_phi = logit_inverse_linkfun(etar,  etad, link); //logit(r)
    //phi=r;
  }

//  //Calculate binomial probability parameter (phi)
//  vector<Type> phi(yobs.size());
//  for (int i = 0; i < yobs.size(); i++)
//    phi(i) = phifun(etar(i), etad(i), etap(i), link, cover);

  // Observation likelihood
  for (int i=0; i < yobs.size(); i++){
    if ( !selfisher::isNA(yobs(i)) ) {
      jnll -= dbinom_robust(yobs(i) * total(i), total(i), logit_phi(i), true);
      SIMULATE{yobs(i) = rbinom(total(i), invlogit(logit_phi(i)));}
      //jnll -= dbinom(yobs(i) * total(i), total(i), phi(i), true);
      //SIMULATE{yobs(i) = rbinom(total(i), phi(i));}

    }
  }

  // Report / ADreport / Simulate Report
  vector<matrix<Type> > corrr(termsr.size());
  vector<vector<Type> > sdr(termsr.size());
  for(int i=0; i<termsr.size(); i++){
    // NOTE: Dummy terms reported as empty
    if(termsr(i).blockNumTheta > 0){
      corrr(i) = termsr(i).corr;
      sdr(i) = termsr(i).sd;
    }
  }
  vector<matrix<Type> > corrp(termsp.size());
  vector<vector<Type> > sdp(termsp.size());
  for(int i=0; i<termsp.size(); i++) {
    // NOTE: Dummy terms reported as empty
    if(termsp(i).blockNumTheta > 0){
      corrp(i) = termsp(i).corr;
      sdp(i) = termsp(i).sd;
    }
  }

  REPORT(corrr);
  REPORT(sdr);
  REPORT(corrp);
  REPORT(sdp);
  SIMULATE{ REPORT(yobs);}

  // For predict
  //vector<Type> r(etar.size());

  switch(pPredictCode) {
  case response_ppredictcode:
    r = invlogit(logit_phi); // Account for relative fishing power
    //r=phi;
    break;
  case selection_ppredictcode:
    for (int i = 0; i < r.size(); i++)
      r(i) = inverse_linkfun(etar(i), etad(i), link);
    break;
  case prob_ppredictcode:
    r = invlogit(etap); // Predict relative fishing power
    break;
  case ratio_ppredictcode:
    for (int i = 0; i < r.size(); i++)
      r(i) = inverse_linkfun(etar(i), etad(i), link)/(1-inverse_linkfun(etar(i), etad(i), link));
    break;
  default:
    error("Invalid 'PredictCode'");
  }

  whichPredict -= 1; // R-index -> C-index
  vector<Type> mu_predict = r(whichPredict);
  REPORT(mu_predict);
  // ADREPORT expensive for long vectors - only needed by predict()
  // method.
  if (doPredict) ADREPORT(mu_predict);

  if((Lpflag!=0) && (Xr.cols()==2) && (Zr.cols()==0) && (Xd.cols()<2)) { //length only model (for now)
    vector<Type> retp(3); //retention probability
    vector<int> SRcalcs(2); //store indecies of .25 and .75 in retp
    switch(Lpflag) {
    case 0: //none
      break;
    case 1: //basic
      retp(0)=Type(0.25);
      retp(1)=Type(0.5);
      retp(2)=Type(0.75);
      SRcalcs(0)=0;
      SRcalcs(1)=2;
      break;
    case 2: //full
      retp.resize(19);
      retp << Type(0.05),Type(0.10),Type(0.15),Type(0.20),Type(0.25),Type(0.30),Type(0.35),Type(0.40),Type(0.45),Type(0.50),Type(0.55),Type(0.60),Type(0.65),Type(0.70),Type(0.75),Type(0.80),Type(0.85),Type(0.90),Type(0.95);
      SRcalcs(0)=4;
      SRcalcs(1)=14;
      break;
    case 3: //100
      retp.resize(100);
      for(int i=0; i<100;i++) {
        retp(i) = Type(0.01+0.01*i);
      }
      SRcalcs(0)=24;//index of retp==0.25
      SRcalcs(1)=74;//index of retp==0.75
      break;
    default:
      error("Invalid 'Lpflag'");
    }
    vector<Type> Lp(retp.size());
    for(int i=0; i<retp.size(); i++) {
      Lp(i) = (linkfun(retp(i), etad(0), link) - betar(0))/betar(1);
    }
	  Type SR = Lp(SRcalcs(1))-Lp(SRcalcs(0));
    ADREPORT(Lp);
    REPORT(retp);
    ADREPORT(SR);
	}

/* TODO This section calculates an L and SR for each observation.
This will be useful for letting L50 and SR depend on covariates,
but it requires postprocessing...
			(1) match back to predictor columns of origianl data
			(2) summarize somehow accounting for sources of uncertainty

  if(Lpflag!=0) {
    matrix<Type> L(etar.size(), retp.size());
    for(int i=0; i<etar.size(); i++) {
      for(int j=0; j<retp.size(); j++) {
        L(i,j) = calcLprob(etar(i), Xr(i, Lindex), betar(Lindex), etad(i), retp(j), link);
      }
    }
      SR = L.col(SRcalcs(1))-L.col(SRcalcs(0));
      ADREPORT(L);
      ADREPORT(SR);
  }
*/

  if(Xp.cols()==1 && Zp.cols()==0) {
		Type p = invlogit(etap)(0);
    ADREPORT(p);
  }
  return jnll;
}


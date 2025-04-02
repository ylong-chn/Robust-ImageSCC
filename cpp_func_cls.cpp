//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 

using namespace arma;

//[[Rcpp::export]]
mat EachColSubtract(const mat& A, const vec& b) {
  mat OutputMat = A.each_col() - b;
  return OutputMat;
}

//[[Rcpp::export]]
mat EachColProd(const mat& A, const vec& b) {
  mat OutputMat = A.each_col() % b;
  return OutputMat;
}

//[[Rcpp::export]]
mat MatHuberDeriv(const double& HuberConst, const mat& ResidMat) {
  mat OutputMat = ResidMat;
  return OutputMat.clamp(-HuberConst, HuberConst);
}

// //[[Rcpp::export]]
// mat cholDecomp(mat X) {
//   return chol(X);
// }
// 
// //[[Rcpp::export]]
// SEXP cholDecompR(SEXP X) {
//   Rcpp::Function f("chol");
//   return f(X, Rcpp::Named("pivot")=FALSE);
// }

// The code for conjugate gradient algorithm below is adapted from the source code of R package `cPCG`
// the code is modified to take arbitrary initial solutions
// [[Rcpp::export]]
arma::vec cgsolve(arma::mat A, arma::vec b, arma::vec x, float tol = 1e-6, int maxIter = 1000) {
  // Function for solving linear equations Ax = b using conjugate gradient
  // get number of columns of A
  // int C = A.n_cols ;
  int R = A.n_rows ;
  // initiate solution x as zeros
  // arma::vec x(C) ;
  // x.zeros() ;
  
  arma::vec r = b - A * x;
  arma::vec p = r;
  double rs_old = sum( r % r );
  double rs_new=1;
  
  arma::vec Ap(R);
  double alpha, beta;
  // vector version of alpha
  // arma::vec alphaVec(1);
  
  //start of iteration 
  for(int iter = 0; (iter < maxIter) && (rs_new > tol); iter++){
    
    Ap = A * p;
    
    alpha = rs_old / sum(p % Ap);
    // alphaVec.fill(alpha); 
    x += alpha * p;
    r -= alpha * Ap;
    rs_new = sum(r % r);
    beta = rs_new / rs_old; 
    
    p = r + beta * p;
    rs_old = rs_new;
    if (iter >= maxIter){
      Rcpp::Rcout << "cg did not converge." << endl;
    }
  }
  
  return x;
  
} 

class BpstHuberReg {
private:
  mat WeightMat, ResidualMat, PsiValMat;
  vec ParUpdate, DiffUpdate;
  int N, n, IterCounter;
  
  // member function 
  // update weight
  void WeightUpdate(vec& ThetaUpdate, double& HuberConst) {
    ResidualMat = EachColSubtract(YMat, BTildeMat*ThetaUpdate); // calculate residuals
    PsiValMat = MatHuberDeriv(HuberConst, ResidualMat); // calculate psi func values
    WeightMat = PsiValMat/ResidualMat; // update weights
  }
  
  void WeightUpdateWeighted(vec& ThetaUpdate, double& HuberConst, vec& SampleWeight) {
    ResidualMat = EachColSubtract(YMat, BTildeMat*ThetaUpdate); // calculate residuals
    PsiValMat = MatHuberDeriv(HuberConst, ResidualMat); // calculate psi func values
    WeightMat = PsiValMat/ResidualMat; // update weights
    WeightMat = EachColProd(WeightMat.t(), SampleWeight); 
    WeightMat = WeightMat.t();
  }
  
  // compute M regression estimate for a single lambda 
  void MEstSingleLambdaInternal(vec& ThetaUpdate, double& LambVal, double& HuberConst, int& MaxIter, double& Tolerance) {
    
    IterCounter = 0;
    
    while(IterCounter < MaxIter) {
      IterCounter++;
      
      // iterations for M estimation
      WeightUpdate(ThetaUpdate, HuberConst); 
      ParUpdate = cgsolve(
        BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambVal * DMat,
        BTildeTransMat * sum(WeightMat % YMat, 1),
        ThetaUpdate,
        1e-10,
        1000); // update coefficients 
      
      // check stopping criterion
      DiffUpdate = ParUpdate - ThetaUpdate;
      if(sqrt(dot(DiffUpdate, DiffUpdate)/DiffUpdate.size()) > Tolerance) {
        ThetaUpdate = ParUpdate;
      } else{break;}
    }
  }
  
  // compute M regression estimate for a single lambda given sample weight vector
  void MEstSingleLambdaInternalWeighted(
      vec& ThetaUpdate, 
      vec& SampleWeight, 
      double& LambVal, 
      double& HuberConst, 
      int& MaxIter, 
      double& Tolerance) {
    
    IterCounter = 0;
    
    while(IterCounter < MaxIter) {
      IterCounter++;
      
      // iterations for M estimation
      WeightUpdateWeighted(ThetaUpdate, HuberConst, SampleWeight); 
      ParUpdate = cgsolve(
        BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambVal * DMat,
        BTildeTransMat * sum(WeightMat % YMat, 1),
        ThetaUpdate,
        1e-10,
        1000); // update coefficients 
      
      // check stopping criterion
      DiffUpdate = ParUpdate - ThetaUpdate;
      if(sqrt(dot(DiffUpdate, DiffUpdate)/DiffUpdate.size()) > Tolerance) {
        ThetaUpdate = ParUpdate;
      } else{break;}
    }
  }
  
public:
  mat BTildeMat, BTildeTransMat, DMat, YMat, HatMat;
  
  // constructor
  BpstHuberReg(mat& BTildeMat_, mat& DMat_, mat& YMat_) {
    BTildeMat = BTildeMat_;
    DMat = DMat_;
    YMat = YMat_; // N by n, which is the transpose of original Y matrix
    BTildeTransMat = BTildeMat.t();
    N = YMat.n_rows;
    n = YMat.n_cols;
  }
  
  // public functions 
  
  // external function to compute M regression estimate for a single lambda 
  Rcpp::List MEstSingleLambda(double& LambVal, double& HuberConst, int& MaxIter, double& Tolerance) {
    
    vec ThetaUpdate = solve(BTildeTransMat * BTildeMat + 2 * LambVal * DMat, BTildeTransMat * mean(YMat, 1));
    
    MEstSingleLambdaInternal(ThetaUpdate, LambVal, HuberConst, MaxIter, Tolerance);
    
    return Rcpp::List::create(
      Rcpp::Named("Iteration") = IterCounter,
      Rcpp::Named("Coeff") = ParUpdate
    );
  }
  
  // perform cv given a lambda grid, using weighted gcv 
  Rcpp::List MEstCrossValidation(vec& LambVec, double& HuberConst, int& MaxIter, double& Tolerance) {
    
    int LambCount = LambVec.n_elem;
    double GCVScoreRun, GCVScoreOut = -1.0, LambdaSelected = 0.0, LambCur, HatMatTrace;
    vec IterationCounts = zeros(size(LambVec)),
      GCVOut = zeros(size(LambVec)), 
      ThetaUpdate  = solve(BTildeTransMat * BTildeMat + 2 * LambVec(0) * DMat, BTildeTransMat * mean(YMat, 1)), 
      ThetaOut = ThetaUpdate;
    
    for(int i = 0; i < LambCount; ++i) {
      
      LambCur = LambVec(i);
      MEstSingleLambdaInternal(ThetaUpdate, LambCur, HuberConst, MaxIter, Tolerance);
      
      IterationCounts(i) = IterCounter;
      
      WeightUpdate(ThetaUpdate, HuberConst); 
      
      HatMatTrace = trace(
        solve(BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambCur * DMat, 
              BTildeTransMat * diagmat(sum(WeightMat, 1))) * BTildeMat)/n/N;
      
      GCVScoreRun = mean(sum(WeightMat % pow(ResidualMat, 2)/(1-HatMatTrace)/(1-HatMatTrace)));
      
      GCVOut(i) = GCVScoreRun;
      if(GCVScoreRun < GCVScoreOut || GCVScoreOut < 0) {
        GCVScoreOut = GCVScoreRun;
        LambdaSelected = LambCur;
        ThetaOut = ParUpdate;
      }
    }
    
    // to output weight matrix
    WeightUpdate(ThetaOut, HuberConst); 
    
    return Rcpp::List::create(
      Rcpp::Named("Iteration") = IterationCounts,
      Rcpp::Named("GCVall") = GCVOut,
      Rcpp::Named("Coeff") = ThetaOut,
      Rcpp::Named("GCV") = GCVScoreOut,
      Rcpp::Named("LambdaCV") = LambdaSelected,
      Rcpp::Named("WeightMatFinal") = WeightMat
    );
  }
  
  // perform cv given a lambda grid, using weighted gcv -- weighted version 
  Rcpp::List MEstCrossValidationWeighted(
      vec& SampleWeight,
      vec& LambVec, 
      double& HuberConst, 
      int& MaxIter, 
      double& Tolerance) {
    
    int LambCount = LambVec.n_elem;
    double GCVScoreRun, GCVScoreOut = -1.0, LambdaSelected = 0.0, LambCur, HatMatTrace;
    vec IterationCounts = zeros(size(LambVec)),
      GCVOut = zeros(size(LambVec)), 
      ThetaUpdate  = solve(BTildeTransMat * BTildeMat + 2 * LambVec(0) * DMat, BTildeTransMat * mean(YMat, 1)), 
      ThetaOut = ThetaUpdate;
    
    for(int i = 0; i < LambCount; ++i) {
      
      LambCur = LambVec(i);
      MEstSingleLambdaInternalWeighted(ThetaUpdate, SampleWeight, LambCur, HuberConst, MaxIter, Tolerance);
      
      IterationCounts(i) = IterCounter;
      
      WeightUpdateWeighted(ThetaUpdate, HuberConst, SampleWeight); 
      
      HatMatTrace = trace(
        solve(BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambCur * DMat, 
              BTildeTransMat * diagmat(sum(WeightMat, 1))) * BTildeMat)/n/N;
      
      GCVScoreRun = mean(sum(WeightMat % pow(ResidualMat, 2)/(1-HatMatTrace)/(1-HatMatTrace)));
      
      GCVOut(i) = GCVScoreRun;
      if(GCVScoreRun < GCVScoreOut || GCVScoreOut < 0) {
        GCVScoreOut = GCVScoreRun;
        LambdaSelected = LambCur;
        ThetaOut = ParUpdate;
      }
    }
    
    WeightUpdateWeighted(ThetaUpdate, HuberConst, SampleWeight); 
    
    return Rcpp::List::create(
      Rcpp::Named("Iteration") = IterationCounts,
      Rcpp::Named("GCVall") = GCVOut,
      Rcpp::Named("Coeff") = ThetaOut,
      Rcpp::Named("GCV") = GCVScoreOut,
      Rcpp::Named("LambdaCV") = LambdaSelected,
      Rcpp::Named("WeightMatFinal") = WeightMat
    );
  }
};

// weighted least squares 
class BpstWLSReg {
private:
  mat WeightMat, ResidualMat;
  vec ParUpdate;
  int N, n;
  
  // member function 
  
public:
  mat BTildeMat, BTildeTransMat, DMat, YMat, HatMat;
  
  // constructor
  BpstWLSReg(mat& BTildeMat_, mat& DMat_, mat& YMat_, mat& WeightMat_) {
    BTildeMat = BTildeMat_;
    DMat = DMat_;
    YMat = YMat_; // N by n, which is the transpose of original Y matrix
    BTildeTransMat = BTildeMat.t();
    WeightMat = WeightMat_;
    N = YMat.n_rows;
    n = YMat.n_cols;
  }
  
  // public functions 
  
  // external function to compute WLS regression estimate for a single lambda 
  Rcpp::List WLSEstSingleLambda(double& LambVal) {
    
    ParUpdate = solve(
      BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambVal * DMat,
      BTildeTransMat * sum(WeightMat % YMat, 1));
    
    return Rcpp::List::create(
      Rcpp::Named("Coeff") = ParUpdate
    );
  }
  
  // perform cv given a lambda grid, using weighted gcv 
  Rcpp::List WLSEstCrossValidation(vec& LambVec) {
    
    int LambCount = LambVec.n_elem;
    double GCVScoreRun, GCVScoreOut = -1.0, LambdaSelected = 0.0, LambCur, HatMatTrace;
    vec IterationCounts = zeros(size(LambVec)),
      GCVOut = zeros(size(LambVec)), 
      ParUpdate  = solve(
        BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambVec(0) * DMat,
        BTildeTransMat * sum(WeightMat % YMat, 1)
      ),
      ThetaOut = ParUpdate;
    
    for(int i = 0; i < LambCount; ++i) {
      
      LambCur = LambVec(i);
      
      ParUpdate = cgsolve(
        BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambCur * DMat,
        BTildeTransMat * sum(WeightMat % YMat, 1),
        ParUpdate,
        1e-10,
        1000);
      
      HatMatTrace = trace(
        solve(BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambCur * DMat, 
              BTildeTransMat * diagmat(sum(WeightMat, 1))) * BTildeMat)/n/N;
      
      ResidualMat = EachColSubtract(YMat, BTildeMat*ParUpdate); 
      GCVScoreRun = mean(sum(WeightMat % pow(ResidualMat, 2)/(1-HatMatTrace)/(1-HatMatTrace)));
      
      GCVOut(i) = GCVScoreRun;
      if(GCVScoreRun < GCVScoreOut || GCVScoreOut < 0) {
        GCVScoreOut = GCVScoreRun;
        LambdaSelected = LambCur;
        ThetaOut = ParUpdate;
      }
    }
    
    
    
    
    return Rcpp::List::create(
      Rcpp::Named("GCVall") = GCVOut,
      Rcpp::Named("Coeff") = ParUpdate,
      Rcpp::Named("GCV") = GCVScoreOut,
      Rcpp::Named("LambdaCV") = LambdaSelected
    );
  }
  
  Rcpp::List WLSEstCrossValidation2(vec& LambVec) {
    
    int LambCount = LambVec.n_elem;
    double GCVScoreRun, GCVScoreOut = -1.0, LambdaSelected = 0.0, LambCur, HatMatTrace;
    vec IterationCounts = zeros(size(LambVec)),
      GCVOut = zeros(size(LambVec)), 
      ParUpdate  = solve(
        BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambVec(0) * DMat,
        BTildeTransMat * sum(WeightMat % YMat, 1)
      ),
      ThetaOut = ParUpdate;
    
    for(int i = 0; i < LambCount; ++i) {
      
      LambCur = LambVec(i);
      
      ParUpdate = solve(
        BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambCur * DMat,
        BTildeTransMat * sum(WeightMat % YMat, 1));
      
      HatMatTrace = trace(
        solve(BTildeTransMat * diagmat(sum(WeightMat, 1)) * BTildeMat + 2 * LambCur * DMat, 
              BTildeTransMat * diagmat(sum(WeightMat, 1))) * BTildeMat)/n/N;
      
      ResidualMat = EachColSubtract(YMat, BTildeMat*ParUpdate); 
      GCVScoreRun = mean(sum(WeightMat % pow(ResidualMat, 2)/(1-HatMatTrace)/(1-HatMatTrace)));
      
      GCVOut(i) = GCVScoreRun;
      if(GCVScoreRun < GCVScoreOut || GCVScoreOut < 0) {
        GCVScoreOut = GCVScoreRun;
        LambdaSelected = LambCur;
        ThetaOut = ParUpdate;
      }
    }
    
    
    
    
    return Rcpp::List::create(
      Rcpp::Named("GCVall") = GCVOut,
      Rcpp::Named("Coeff") = ParUpdate,
      Rcpp::Named("GCV") = GCVScoreOut,
      Rcpp::Named("LambdaCV") = LambdaSelected
    );
  }
};

//[[Rcpp::export]]
Rcpp::List MEstLambda(
    mat& BTilde,
    mat& TransY,
    mat& PenaltyMat,
    double& LambVal,
    double& HuberConst,
    int& MaxIter,
    double& Tolerance
) {
  
  BpstHuberReg BpstHuberObj = BpstHuberReg(BTilde, PenaltyMat, TransY);
  return BpstHuberObj.MEstSingleLambda(LambVal, HuberConst, MaxIter, Tolerance);
  
}

//[[Rcpp::export]]
Rcpp::List MEstLambdaPath(
    mat& BTilde,
    mat& TransY,
    mat& PenaltyMat,
    vec& LambVec,
    double& HuberConst,
    int& MaxIter,
    double& Tolerance
) {
  
  BpstHuberReg BpstHuberObj = BpstHuberReg(BTilde, PenaltyMat, TransY);
  return BpstHuberObj.MEstCrossValidation(LambVec, HuberConst, MaxIter, Tolerance);
  
}

//[[Rcpp::export]]
Rcpp::List MEstLambdaPathWeighted(
    mat& BTilde,
    mat& TransY,
    mat& PenaltyMat,
    vec& SampleWeight,
    vec& LambVec,
    double& HuberConst,
    int& MaxIter,
    double& Tolerance
) {
  
  BpstHuberReg BpstHuberObj = BpstHuberReg(BTilde, PenaltyMat, TransY);
  return BpstHuberObj.MEstCrossValidationWeighted(SampleWeight, LambVec, HuberConst, MaxIter, Tolerance);
  
}

//[[Rcpp::export]]
Rcpp::List WLSEstLambdaPath(
    mat& BTilde,
    mat& TransY,
    mat& PenaltyMat,
    mat& WeightMat,
    vec& LambVec
) {
  
  BpstWLSReg BpstWLSObj = BpstWLSReg(BTilde, PenaltyMat, TransY, WeightMat);
  return BpstWLSObj.WLSEstCrossValidation(LambVec);
  
}

//[[Rcpp::export]]
Rcpp::List WLSEstLambdaPath2(
    mat& BTilde,
    mat& TransY,
    mat& PenaltyMat,
    mat& WeightMat,
    vec& LambVec
) {
  
  BpstWLSReg BpstWLSObj = BpstWLSReg(BTilde, PenaltyMat, TransY, WeightMat);
  return BpstWLSObj.WLSEstCrossValidation2(LambVec);
  
}
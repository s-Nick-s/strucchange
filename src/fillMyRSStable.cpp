#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fillMyRSStable(NumericMatrix myRSStable, List & RSStriang, IntegerVector & myIndex, int m, float h) {
  
  NumericMatrix myRSStableCPY = clone(myRSStable);
  CharacterVector indexI(1), indexj(1);
  int n = myIndex.size();
  CharacterVector nms = rownames(myRSStableCPY);
  
  for(int myIndexI = 0; myIndexI < n; ++myIndexI) {
    int i = myIndex[myIndexI];
    
    IntegerVector potIndex =  seq((m - 1)*h, i - h);
    //findRSSbreaksC
    int np = potIndex.size();
    NumericVector breakRSS(np);
    
    for(int pin = 0; pin < np; ++pin) {
      int j = potIndex[pin];
      NumericVector RSStriangI = as<NumericVector>(RSStriang[j]);
      //Rcout << "value : " << RSStriangI[i - j] << "\n";
      
      indexj[0] = j;
      IntegerVector matchedIdxJ = match(indexj, nms);
      
      breakRSS[pin] = myRSStableCPY(matchedIdxJ[0] - 1, 1) + RSStriangI[i - j - 1];
    }
    // end findRSSbreaksC
    int opt = which_min(breakRSS);
    
    indexI[0] = i;
    IntegerVector matchedIdxI = match(indexI, nms);

    myRSStableCPY(matchedIdxI[0] - 1, 2) = potIndex[opt];
    myRSStableCPY(matchedIdxI[0] - 1, 3) = breakRSS[opt];
  }
  return myRSStableCPY;
}
#include <Rcpp.h>
using namespace Rcpp;

// [["Rcpp::export"]]
NumericVector ct_likelihood_fast(const NumericVector &expected, const NumericVector &obs, const NumericVector &pars){
    int total_cts = expected.size();
    NumericVector ret(total_cts);
    const double sd = pars["error"];
    const double max_ct = 40.0;
    const double log_const = log(0.5);
    const double den = sd*M_SQRT2;
    
    for(int i = 0; i < total_cts; ++i){
        // Most titres are between 0 and max_titre, this is the difference in normal cdfs
        if(obs[i] < max_ct && obs[i] >= 0.0){
            ret[i] = log_const + log((erf((obs[i] + 1.0 - expected[i]) / den) -
                erf((obs[i] - expected[i]) / den)));    
            // For titres above the maximum, 
        } else if(obs[i] >= max_ct) {
            ret[i] = log_const + log(erfc((max_ct - expected[i])/den));
        } else {
            ret[i] = log_const + log(1.0 + erf((0.0 - expected[i])/den));
        }
    }
    return(ret);
} 
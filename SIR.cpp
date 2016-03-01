#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]

List tauleapCpp(List params) {
  
    // chained operations are tricky in cpp
    // pull out list w/in list into its own object
    List init = params["init"];
    
    // use Rcpp as() function to "cast" R vector to cpp scalar
    int nsteps = as<int>(params["nsteps"]);

    // initialize each state vector in its own vector
    // set all vals to initial vals
    //
    // I use doubles (NumericVector) rather than 
    // ints (IntegerVector), since rpois returns double,
    // and the domain of double is a superset of int
    NumericVector SS(nsteps, init["S"]);
    NumericVector II(nsteps, init["I"]);
    NumericVector RR(nsteps, init["R"]);
    NumericVector NN(nsteps, init["pop"]);
    
    // fill time w/zeros
    NumericVector time(nsteps);

    // pull out params for easy reading 
    double nu = params["nu"];
    double mu = params["mu"];
    double beta = params["beta"];
    double gamma = params["gamma"];
    double tau = params["tau"];

    // Calculate the number of events for each step, update state vectors
    for (int istep = 0; istep < (nsteps-1); istep++) {
      
        // pull out this step's scalars for easier reading
        // and to avoid compiler headaches
        double iS = SS[istep];
        double iI = II[istep];
        double iR = RR[istep];
        double iN = NN[istep];
        
        /////////////////////////
        // State Equations
        /////////////////////////
        
        // R::rpois always returns a single value
        // to return multiple (e.g. Integer/NumericVector, 
        // use Rcpp::rpois(int ndraw, param) and friends
        double births = R::rpois(nu*iN*tau);
        
        // Prevent negative states
        double Sdeaths = std::min(iS, R::rpois(mu*iS*tau));
        double maxtrans = R::rpois(beta*(iI/iN)*iS*tau);
        double transmission = std::min(iS-Sdeaths, maxtrans);
        double Ideaths = std::min(iI, R::rpois(mu*iI*tau));
        double recovery = std::min(iI-Ideaths, R::rpois(gamma*iI*tau));
        double Rdeaths = std::min(iR, R::rpois(mu*iR*tau));
        
        // Calculate the change in each state variable
        double dS = births-Sdeaths-transmission;
        double dI = transmission-Ideaths-recovery;
        double dR = recovery-Rdeaths;
        
        // Update next timestep
        SS[istep+1] = iS + dS;
        II[istep+1] = iI + dI;
        RR[istep+1] = iR + dR;
        
        // Sum population
        NN[istep+1] = iS + iI + iR + dS + dI + dR;
        
        // time in fractional years (ie units parameters are given in)
        time[istep+1] = (istep+1)*tau;
    }
    
    // Return results as data.frame
    DataFrame sim = DataFrame::create(
        Named("time") = time,
        Named("S") = SS,
        Named("I") = II,
        Named("R") = RR,
        Named("N") = NN
    );
    return sim;
};
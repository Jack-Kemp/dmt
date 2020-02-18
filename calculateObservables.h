#ifndef CALCULATE_OBSERVABLES_H
#define CALCULATE_OBSERVABLES_H

#include "itensor/all.h"
#include "DMT.h"


using namespace itensor;

//DMT Functions
Complex
calculateExpectation( const char * op_name, const int& site_i,
		      const DMT& dmt){
  auto op_i = dmt.siteOp(op_name, site_i);
  auto exp = dmt.traceLeftOf(site_i)*op_i*dmt.rho(site_i)*dmt.traceRightOf(site_i);
  return eltC(exp)/dmt.trace();
}

Complex
calculateTwoPoint(const char * op_name_i, const int& site_i,
		   const char * op_name_j, const int& site_j,
		      const DMT& dmt){
  if (site_i == site_j)
    Error("Two-point correlator must be on different sites.");
  auto op_i = dmt.siteOp(op_name_i, site_i);
  auto op_j = dmt.siteOp(op_name_j, site_j);  
  auto left = dmt.traceLeftOf(site_i)*op_i*dmt.rho(site_i);
  for (int i = site_i +1; i < site_j; i++){
    left *= dmt.traceOf(i);
  }
  auto right = op_j*dmt.rho(site_j)*dmt.traceRightOf(site_j);
  return eltC(left*right)/dmt.trace();
}

ITensor
reducedDensityMatrix(const DMT& dmt, int siteStart, int siteEnd){
  auto left = dmt.traceLeftOf(siteStart);
  auto right = dmt.traceRightOf(siteEnd);
  for (int i = siteStart; i <= siteEnd; i++){
    left *= dmt.rho(i);
  }
  return left*right; 
}

Real secondRenyiEntropyHalfSystem(const DMT& dmt){
  auto rdm = reducedDensityMatrix(dmt, floor(dmt.len()/2)+1, dmt.len())/dmt.trace();
  auto rdm2_trace = rdm*dag(rdm);
  //auto rdm2_trace = rdm2*delta(dag(rdm2.inds()));
  return abs(log(eltC(rdm2_trace)));
}


  
  
  


//Psi Functions


#endif

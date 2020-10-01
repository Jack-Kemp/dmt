#ifndef CALCULATE_OBSERVABLES_H
#define CALCULATE_OBSERVABLES_H

#include <itensor/all.h>
#include "DMT.h"


using namespace itensor;

//DMT Functions


//Calculate expectation of op_name on site i
Complex
calculateExpectation( const char * op_name, const int& site_i,
		      const DMT& dmt){
  auto op_i = dmt.siteOp(op_name, site_i);
  auto exp = dmt.traceLeftOf(site_i)*op_i*dmt.rho(site_i)*dmt.traceRightOf(site_i);
  return eltC(exp)/dmt.trace();
}

//Calculate expectation of MPO op.

//As part of the calculation this must be converted to DMT basis,
//(a no op if not vector basis). If you are going to call with same
//MPO many times, an optimisation is to preconvert and then set
//convertToDMTBasis = False.
Complex
calculateExpectation( const MPO& op,
		      const DMT& dmt,
		      bool convertToDMTBasis = true){
  return dmt.trace(op, convertToDMTBasis)/dmt.trace();
}

//Calculate expectation <(op_name on site i) * (op_name on site j)>
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

//Calculate the reduced matrix from siteStart to siteEnd inclusive,
//i.e. trace out all sites NOT between siteStart to siteEnd inclusive.
//This returns a TENSOR, not an MPO, so if the region
//not traced out is large this will be very LARGE.
ITensor
reducedDensityMatrix(const DMT& dmt, int siteStart, int siteEnd){
  auto left = dmt.traceLeftOf(siteStart);
  auto right = dmt.traceRightOf(siteEnd);
  for (int i = siteStart; i <= siteEnd; i++){
    left *= dmt.rho(i);
  }
  return left*right; 
}

//Because this takes the square of the density matrix, it is signifcantly
//more taxing to compute than the above quantities. Only measure if you need
//it.
Real secondRenyiEntropyHalfSystem(const DMT& dmt){
  auto centre = floor(dmt.len()/2)+1;
  auto left = dmt.traceLeftOf(centre);
  left *= prime(left);
  for (int i = centre; i <= dmt.len(); i++){
    left *= dmt.rho(i)*prime(dag(dmt.rho(i)), "Link");
  }
  return abs(log(eltC(left)))/(dmt.trace()*dmt.trace());
}


  
  
  


//Psi Functions


#endif

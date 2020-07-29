#ifndef CALCULATE_OP_OBSERVABLES_H
#define CALCULATE_OP_OBSERVABLES_H

#include "itensor/all.h"
#include "DMT.h"


using namespace itensor;

std::vector<Real> sizeDistribution(DMT& dmt){
  
  auto d = dim(dmt.sites()(1));
  auto L = dmt.len();
  
  std::vector<Complex> Fs;
  for (int k = 0; k < L+1; k++)
    {
      auto z = exp(Complex(0.0, 2 * M_PI * k / (L+1)));
      auto cumul = ITensor(1);
      for (int i =1; i <= L;  i++)
      {
	auto idPart = (1-z)/d*dmt.traceOf(i);

	idPart *= prime(idPart);
	auto Wpart = z*(dmt.rho(i)*prime(swapPrime(dmt.rho(i),0,1, "Site"), "Link"));
	if (i == 1 or i == 2){
	PrintData(idPart);
	PrintData(Wpart);
	}
	cumul *= idPart + Wpart;
      }
      PrintData(cumul/pow(d,L));
      Fs.push_back(eltC(cumul)/pow(d,L));
  }
  std::vector<Real> pls;
  for(int l = 0; l < L+1; l++)
    {
      auto cumul = Complex(0.0, 0.0);
      for (int k = 0; k < L+1; k++)
	cumul += Fs[k]*exp(Complex(0.0, -2 * M_PI * k* l / (L+1)));
      pls.push_back(1/(L+1) * cumul.real());
    }
  return pls;
}

//DMT Functions
Complex
calculateSingleOpTrace( const char * op_name, const int& site_i,
		      const DMT& dmt){
  auto op_i = dmt.siteOp(op_name, site_i);
  auto exp = dmt.traceLeftOf(site_i)*op_i*dmt.rho(site_i)*dmt.traceRightOf(site_i);
  return eltC(exp);///pow(d,L);
}

#endif

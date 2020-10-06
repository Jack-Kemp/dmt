#include <itensor/all.h>

using namespace itensor;

//Calculate expectation of op_name on site i
Complex
calculateExpectation( const char * op_name, const int& site_i,
		      MPS & psi, const SiteSet & sites){
  auto op_i = op(sites, op_name, site_i);
  psi.position(site_i);
  auto ket = psi(site_i);
  auto bra = dag(prime(ket,"Site"));
  return eltC(bra*op_i*ket);
}

//Calculate expectation of MPO op.

//As part of the calculation this must be converted to DMT basis,
//(a no op if not vector basis). If you are going to call with same
//MPO many times, an optimisation is to preconvert and then set
//convertToDMTBasis = False.
Complex
calculateExpectation( const MPO& op,
		      const MPS& psi){
  return innerC(psi, op, psi);
}


//Calculate expectation <(op_name on site i) * (op_name on site j)>
//If i > j, swaps the terms. If operators on spatially seperated sites
//do not commute, (i.e. fermions), this swap is important to understand.
Complex
calculateTwoPoint(const char * op_name_i, int site_i,
		  const char * op_name_j, int site_j,
		  MPS & psi, const SiteSet & sites){


  auto op_i = op(sites,op_name_i, site_i);
  auto op_j = op(sites,op_name_j, site_j);
  if (site_i > site_j){
    std::swap(op_i, op_j);
    std::swap(site_i, site_j);
  }

  //'gauge' the MPS to site i
  //any 'position' between i and j, inclusive, would work here
  psi.position(site_i); 

  //Create the bra/dual version of the MPS psi
  auto psidag = dag(psi);

  //Prime the link indices to make them distinct from
  //the original ket links
  psidag.prime("Link");

  //index linking i-1 to i:
  auto li_1 = leftLinkIndex(psi,site_i);

  auto C = prime(psi(site_i),li_1)*op_i;
  C *= prime(psidag(site_i),"Site");
  for(int k = site_i+1; k < site_j; ++k)
    {
      C *= psi(k);
      C *= psidag(k);
    }
  //index linking j to j+1:
  auto lj = rightLinkIndex(psi,site_j);

  C *= prime(psi(site_j),lj)*op_j;
  C *= prime(psidag(site_j),"Site");

  return eltC(C);
}

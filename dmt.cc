//#include "calculate_observables.h"
#if 1
#include "input_output_utilities.h"
#include "density_matrix.h"
#include "mpo_tebd.h"
#include "itensor/all.h"
#include <itensor/mps/bondgate.h>
#include<chrono>
#include<functional>



using namespace itensor;
using std::vector;

int main(int argc, char* argv[])
{
  struct Flags;
  int N = 10; //number of sites
  Real tstep = 0.01; //time step (smaller is generally more accurate)
  Real ttotal = 1.00; //total time to evolve
  Real cutoff = 1E-8; //truncation error cutoff when restoring MPS form
  int maxDim = 64;
  bool vectorized = true;
  //Define a site set object "sites" which lets us
  //easily obtain Site indices defining our Hilbert space
  //and S=1/2 single-site operators
  auto sites = SpinHalf(N, {"ConserveQNs=",false});

  //Make initial MPS psi to be in the Neel state
  auto state = InitState(sites);
  for(auto j : range1(N))
    {
      state.set(j,j%2==1?"Up":"Dn");
    }
  auto psi = MPS(state);
  DMTDensityMatrix dmt;
  MPO & rho = dmt.rho();
  rho = projector(psi);
  dmt.presRange(1);
  dmt.sites(sites);
  if (vectorized)
    dmt.vec();
  PrintData(siteInds<MPO>(rho,1));

  //Create a std::vector (dynamically sizeable array)
  //to hold the Trotter gates
  auto mpsgates = vector<BondGate>();
  auto dmtgates = vector<BondGate>();

  auto ampo = AutoMPO(sites);
 
  //Create the gates exp(-i*tstep/2*hterm)
  //and add them to gates
  for(int b = 1; b < N; ++b)
    {
      auto hterm = op(sites,"Sz",b)*op(sites,"Sz",b+1);
      ampo += "Sz", b, "Sz", b+1;
      hterm += 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
      ampo += 0.5, "S+", b, "S-", b+1;
      hterm += 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);
      ampo += 0.5, "S-", b, "S+", b+1;
      
      mpsgates.push_back(BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm));
      dmtgates.push_back(dmt.calcGate(hterm, tstep, b));
      //dmtgates.emplace_back(sites, b, b+1, dmtgate);
    }
  auto H = toMPO(ampo);
  //Create the gates exp(-i*tstep/2*hterm) in reverse
  //order (to get a second order Trotter breakup which
  //does a time step of "tstep") and add them to gates
  for(int b = N-1; b >= 1; --b)
    {
      auto hterm = op(sites,"Sz",b)*op(sites,"Sz",b+1);
      hterm += 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
      hterm += 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);
      mpsgates.push_back(BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm));
      dmtgates.push_back(dmt.calcGate(hterm, tstep, b));
      //dmtgates.emplace_back(sites, b, b+1, dmtgate);
    }

  //Save initial state;
  auto psi0 = psi;
  //PrintData(psi);

  //Time evolve, overwriting psi when done
  gateTEvol(dmtgates,ttotal,tstep,dmt,{"Cutoff",cutoff,
				       "Verbose",true,
				       "MaxDim", maxDim,
				       "UseSVD", true,
				       "DoNormalize", false});
  gateTEvol(mpsgates,ttotal,tstep,psi,{"Cutoff=",cutoff,"Verbose=",true});
  

  printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));
  printfln("Maximum MPO bond dimension after time evolution is %d",maxLinkDim(dmt.rho()));

  //Print overlap of final state with initial state
  //(Will be complex so using innerC which can return complex);
  if (vectorized)
    dmt.unvec();
  auto overlap = innerC(psi,psi0);
  Print(std::norm(overlap));
  Print(traceC(rho,swapPrime(projector(psi0),0,1)));
  auto initialE = innerC(psi0, H, psi0);
  Print(initialE);
  auto E = innerC(psi, H, psi);
  Print(E);
  auto Erho = traceC(rho, H);
  Print(Erho);

  return 0;
}
#endif

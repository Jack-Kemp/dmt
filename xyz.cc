#if 1
#include "calculateObservables.h"
#include "inputOutputUtilities.h"
#include "DMT.h"
#include "DMTObserver.h"
#include "gateTevol.h"
#include "itensor/all.h"
#include <itensor/mps/bondgate.h>
#include<chrono>
#include<functional>



using namespace itensor;
using std::vector;

int main(int argc, char* argv[])
{

  //Get input parameters and options --------------------------------
  auto input = InputGroup(argv[1],"input");

  std::map<std::string, Real> params_real;
  std::map<std::string, int> params_int;
  std::map<std::string, bool> params_bool;

  std::vector<std::string> param_names_real = {"tstep",
					       "ttotal",
					       "cutoff",
					       "presCutoff"};
  
  std::vector<std::string> param_names_int = {"N", 
					      "maxDim",
					      "presRange"};

  std::vector<std::string> param_names_bool = {"vectorized",
					       "conserveQNs",
					       "normalize",
					       "useSVD",
					       "writeMagz",
					       "writeMagx",
					       "writeCorrzz",
					       "writeCorrxx",
					       "writeS2"
					       "writeMaxDim"};
  
  for(auto const& x : param_names_real)
      params_real[x] = input.getReal(x);
  for(auto const& x : param_names_int)
    params_int[x] = input.getInt(x);
  for(auto const& x : param_names_bool)
    params_bool[x] = input.getYesNo(x);

  const int N = params_int["N"];
  const Real tstep = params_real["tstep"];
  const bool vectorized = params_bool["vectorized"];

  std::string outputDir = input.getString("outputDir"); 
  std::string outputName = input.getString("outputName");

  
  //End input-output options -----------------------------------------

  

  //Set up initial state----------------------------------------------

  auto sites = SpinHalf(N, {"ConserveQNs=", params_bool["conserveQNs"]});

  DMT dmt;
  MPO & rho = dmt.rhoRef();
  dmt.presRange(params_int["presRange"]);
  dmt.sites(sites);
  
  
  /* From projector of wavefunction

  auto state = InitState(sites);
  for(auto j : range1(N))
    {
      state.set(j,j%2==1?"Up":"Dn");
    }
  auto psi = MPS(state);
  //Save initial state;
  auto psi0 = psi;
  rho = projector(psi);

  */

  //From product MPO of operators
  for(int j = 1; j <= N; ++j)
    {
      rho.ref(j) = sites.op("Id",j)+ 2*sites.op("Sz",j);
    }
  putMPOLinks(rho);
  
   
  if (vectorized)
    dmt.vec();


  //Set up Hamilitonian----------------------------------------------
  auto mpsgates = vector<BondGate>();
  auto dmtgates = vector<BondGate>();
  auto ampo = AutoMPO(sites);
 
  for(int b = 1; b < N; ++b)
    {
      auto hterm = op(sites,"Sz",b)*op(sites,"Sz",b+1);
      hterm += 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
      hterm += 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);
      
      //mpsgates.push_back(BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm));
      dmtgates.push_back(dmt.calcGate(hterm, tstep, b));
      
      //Construct the full Hamiltonian to calculate energy
      ampo += "Sz", b, "Sz", b+1;
      ampo += 0.5, "S+", b, "S-", b+1;
      ampo += 0.5, "S-", b, "S+", b+1;
    }
  auto H = toMPO(ampo);
  
  for(int b = N-1; b >= 1; --b)
    {
      auto hterm = op(sites,"Sz",b)*op(sites,"Sz",b+1);
      hterm += 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
      hterm += 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);
      mpsgates.push_back(BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm));
      dmtgates.push_back(dmt.calcGate(hterm, tstep, b));
    }

  //Set up measurements-----------------------------------------------
  std::vector<std::vector<Real>> magz, magx, corrzz, corrxx;
  std::vector<Real> S2, chi;

  auto measure = [&](DMT& dmt, Args const & args){
		     magz.push_back(std::vector<Real>());
		     magx.push_back(std::vector<Real>());
		     corrzz.push_back(std::vector<Real>());
		     corrxx.push_back(std::vector<Real>());	     
		     for(int i = 1; i <= N; i++)
		       {
		       if(params_bool["writeMagz"])
			 magz.back().push_back(abs(calculateExpectation("Sz", i, dmt)));
		       if(params_bool["writeMagx"])
			 magx.back().push_back(abs(calculateExpectation("Sx", i, dmt)));
		       if(params_bool["writeCorrzz"])
			 corrzz.back().push_back(abs(calculateTwoPoint("Sz", i, "Sz", (i%N)+1, dmt)));
		       if(params_bool["writeCorrxx"])
			 corrxx.back().push_back(abs(calculateTwoPoint("Sx", i, "Sx", (i%N)+1, dmt)));
		       }
		     if(params_bool["writeS2"])
		       S2.push_back(secondRenyiEntropyHalfSystem(dmt));
		     if(params_bool["writeMaxDim"])
		       chi.push_back(maxLinkDim(dmt.rho())); 							   
		 };

  DMTObserver obs(measure);
  
  //Time evolve-------------------------------------------------------
  gateTEvol(dmtgates,params_real["ttotal"],tstep,dmt, obs,
	    {"Cutoff",params_real["cutoff"],
	     "Verbose",true,
	     "MaxDim", params_int["maxDim"],
	     "UseSVD", params_bool["useSVD"],
	     "DoNormalize", params_bool["normalize"]});

  //gateTEvol(mpsgates,params_real["ttotal"],tstep,psi,{"Cutoff=",params_real["cutoff"],"Verbose=",true});
  

  //printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));
  printfln("Maximum MPO bond dimension after time evolution is %d",maxLinkDim(dmt.rho()));

  std::string label = outputDir + "/" + construct_label(outputName, "",
							params_real,
							params_bool,
							params_int); 

  //Write results-----------------------------------------------------
   if(params_bool["writeMagz"])
     write2DToFile (label + "_magz", magz);
   if(params_bool["writeMagx"])
     write2DToFile (label + "_magx", magx);
   if(params_bool["writeCorrzz"])
     write2DToFile (label + "_corrzz", corrzz);
   if(params_bool["writeCorrxx"])
     write2DToFile (label + "_corrxx", corrxx);
   if(params_bool["writeS2"])
     writeToFile (label + "_S2", S2, true);
   if(params_bool["writeMaxDim"])
     writeToFile (label + "_maxDim", chi, true);


   
   //Some tests
  if (vectorized)
    dmt.unvec();
  //auto overlap = innerC(psi,psi0);
  //Print(std::norm(overlap));
  //Print(traceC(rho,swapPrime(projector(psi0),0,1)));
  //auto initialE = innerC(psi0, H, psi0);
  //Print(initialE);
  //auto E = innerC(psi, H, psi);
  //Print(E);
  auto Erho = traceC(rho, H);
  Print(Erho);

  return 0;
}
#endif

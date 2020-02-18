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
#include<filesystem>
#include"HermitianSpinHalf.h"



using namespace itensor;
using std::vector;
using namespace std::chrono;

int main(int argc, char* argv[])
{

  //Get input parameters and options --------------------------------
  const char* inputName;
  if (argc > 1) {
    inputName = argv[1];
  }
  else{ 
    printfln("Using default params.");
    inputName = "input_params.txt";
  }
  auto input = InputGroup(inputName,"input");
  

  std::map<std::string, Real> paramsReal;
  std::map<std::string, int> paramsInt;
  std::map<std::string, bool> paramsBool;

  std::vector<std::string> paramNamesReal = {"tStep",
					     "tTotal",
					     "cutoff",
					     "presCutoff",
					     "Jz",
					     "Jx",
					     "Jy",
					     "hz",
					     "hx",
					     "hy"};
  
  std::vector<std::string> paramNamesInt = {"N", 
					    "maxDim",
					    "presRange",
					    "nSweeps"};

  std::vector<std::string> paramNamesBool = {"vectorized",
					     "hermitianBasis",
					     "fromPureState",
					     "conserveQNs",
					     "normalize",
					     "useSVD",
					     "writeMagz",
					     "writeMagx",
					     "writeCorrzz",
					     "writeCorrxx",
					     "writeS2",
					     "writeMaxDim",
					     "writeTrace"};
  
  for(auto const& x : paramNamesReal)
    paramsReal[x] = input.getReal(x);
  for(auto const& x : paramNamesInt)
    paramsInt[x] = input.getInt(x);
  for(auto const& x : paramNamesBool)
    paramsBool[x] = input.getYesNo(x);

  const int N = paramsInt["N"];
  const Real tStep = paramsReal["tStep"];
  const Real tSweep = tStep/paramsInt["nSweeps"];
  const bool vectorized = paramsBool["vectorized"];
  const bool hermitianBasis = paramsBool["hermitianBasis"];
  const bool conserveQNs = paramsBool["conserveQNs"];

  Real hx = paramsReal["hx"];
  Real hy = paramsReal["hy"];
  Real hz = paramsReal["hz"];
  Real Jx = paramsReal["Jx"];
  Real Jy = paramsReal["Jy"];
  Real Jz = paramsReal["Jz"];
  Real delta = 0.5*Jx;
  

  std::string outputDir = input.getString("outputDir"); 
  std::string outputName = input.getString("outputName");

  //Sanity check conservation
  if (conserveQNs)
    if( hx != 0 or hy != 0 or (Jx != Jy))
      Error("Hamiltionian does not conserve Z. Check hx,hy = 0 and Jx = Jy");

  //Sanity check bases
  if (vectorized and hermitianBasis)
    Error("Cannot vectorize a vector basis!");
  if (hermitianBasis and paramsBool["fromPureState"])
    Error("Hermitian Basis and from pure state not implemented");

  
  //End input-output options -----------------------------------------

 
  //Set up initial state----------------------------------------------

  SiteSet sites;
  if(paramsBool["hermitianBasis"])
    sites =  HermitianSpinHalf(N);
  else
    sites = SpinHalf(N, {"ConserveQNs=", conserveQNs});

  DMT dmt(sites, {"presRange", paramsInt["presRange"],
		  "vectorBasis", paramsBool["hermitianBasis"],
		  "vectorized", paramsBool["vectorized"]});
  
  
  //From projector of wavefunction
  if(paramsBool["fromPureState"])
    {
      auto state = InitState(sites);
      for(auto j : range1(N))
	{
	  state.set(j,j%2==1?"Up":"Dn");
	}
      auto psi = MPS(state);
      dmt.fromPureState(psi);
    }
  //From product MPO of operators
  else
    {
      for(int j = 1; j <= N; ++j)
	{
	  dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*(2*(j %2)-1)*dmt.stateOp("Sz",j);
	}
    }

  dmt.finishConstruction({"normalize", paramsBool["normalize"]});
			   
  //Set up Hamilitonian----------------------------------------------
  auto mpsgates = vector<BondGate>();
  auto dmtgates = vector<BondGate>();
  auto ampo = AutoMPO(sites);
 
  for(int b = 1; b < N; ++b)
    {
      auto hterm = Jz*dmt.twoSiteOpH("Sz",b,"Sz",b+1);   
      hterm += 0.5*hz*dmt.twoSiteOpH("Sz",b,"Id",b+1);
      hterm += 0.5*hz*dmt.twoSiteOpH("Id",b,"Sz",b+1);

      if(conserveQNs){
	hterm += delta*dmt.twoSiteOpH("S+",b,"S-",b+1);
	hterm += delta*dmt.twoSiteOpH("S-",b,"S+",b+1);
      }
      else{
	hterm += 0.5*hx*dmt.twoSiteOpH("Sx",b,"Id",b+1);
	hterm += 0.5*hy*dmt.twoSiteOpH("Sy",b,"Id",b+1);
	hterm += 0.5*hx*dmt.twoSiteOpH("Id",b,"Sx",b+1);
	hterm += 0.5*hy*dmt.twoSiteOpH("Id",b,"Sy",b+1);
	hterm += Jx*dmt.twoSiteOpH("Sx",b,"Sx",b+1);
	hterm += Jy*dmt.twoSiteOpH("Sy",b,"Sy",b+1);	
      }

      //Single site terms at the boundary
      if (b==N-1)
	{
	  hterm += 0.5*hz*dmt.twoSiteOpH("Id",b,"Sz",b+1);
	  if (not conserveQNs)
	    {
	      hterm += 0.5*hy*dmt.twoSiteOpH("Id",b,"Sy",b+1);
	      hterm += 0.5*hx*dmt.twoSiteOpH("Id",b,"Sx",b+1);
	    }
	}

      if (b==1)
	{
	  hterm += 0.5*hz*dmt.twoSiteOpH("Sz",b,"Id",b+1);
	  if (not conserveQNs)
	    {
	      hterm += 0.5*hy*dmt.twoSiteOpH("Sx",b,"Id",b+1);
	      hterm += 0.5*hx*dmt.twoSiteOpH("Sy",b,"Id",b+1);
	    }
	}

      dmtgates.push_back(dmt.calcGate(hterm, tSweep, b));
      
      //mpsgates.push_back(BondGate(sites,b,b+1,BondGate::tReal,tstep/2.,hterm));
      
      
      // //Construct the full Hamiltonian to calculate energy
      // ampo += "Sz", b, "Sz", b+1;
      // ampo += 0.5, "S+", b, "S-", b+1;
      // ampo += 0.5, "S-", b, "S+", b+1;
    }
  
  for(int b = N-1; b >= 1; --b)
    {
      dmtgates.push_back(dmtgates[b-1]);
      //mpsgates.push_back(mpsgates[b-1]);
    }

  //auto H = toMPO(ampo);

  //Set up measurements-----------------------------------------------
  std::vector<std::vector<Real>> magz, magx, corrzz, corrxx;
  std::vector<Real> S2, chi, traces;

  auto measure = [&](DMT& dmt, Args const & args){
		   magz.push_back(std::vector<Real>());
		   magx.push_back(std::vector<Real>());
		   corrzz.push_back(std::vector<Real>());
		   corrxx.push_back(std::vector<Real>());	     
		   for(int i = 1; i <= N; i++)
		     {
		       if(paramsBool["writeMagz"])
			 magz.back().push_back(calculateExpectation("Sz", i, dmt).real());
		       if(paramsBool["writeMagx"])
			 magx.back().push_back(calculateExpectation("Sx", i, dmt).real());
		       if(paramsBool["writeCorrzz"])
			 corrzz.back().push_back(calculateTwoPoint("Sz", i, "Sz", (i%N)+1, dmt).real());
		       if(paramsBool["writeCorrxx"])
			 corrxx.back().push_back(calculateTwoPoint("Sx", i, "Sx", (i%N)+1, dmt).real());
		     }
		   if(paramsBool["writeS2"])
		     S2.push_back(secondRenyiEntropyHalfSystem(dmt));
		   if(paramsBool["writeMaxDim"])
		     chi.push_back(maxLinkDim(dmt.rho()));
		   if(paramsBool["writeTrace"])
		     traces.push_back(dmt.trace()); 
		 };

  DMTObserver obs(measure);

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  //Time evolve-------------------------------------------------------
  gateTEvol(dmtgates,paramsReal["tTotal"],tStep,dmt, obs,
	    {"Cutoff",paramsReal["cutoff"],
	     "PresCutoff",paramsReal["presCutoff"],
	     "Verbose",true,
	     "MaxDim", paramsInt["maxDim"],
	     "UseSVD", paramsBool["useSVD"],
	     "DoNormalize", paramsBool["normalize"],
	     "nSweeps", paramsInt["nSweeps"] });

  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

  printfln("DMT evolution took %f seconds.", time_span.count());

  //gateTEvol(mpsgates,paramsReal["ttotal"],tstep,psi,{"Cutoff=",paramsReal["cutoff"],"Verbose=",true});
  

  //printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));
  printfln("Maximum MPO bond dimension after time evolution is %d",maxLinkDim(dmt.rho()));

  std::string label = outputDir + "/" + construct_label(outputName, "",
							paramsReal,
							paramsBool,
							paramsInt); 

  //Write results-----------------------------------------------------
  std::filesystem::create_directory(outputDir);
  if(paramsBool["writeMagz"])
    write2DToFile (label + "_magz", magz);
  if(paramsBool["writeMagx"])
    write2DToFile (label + "_magx", magx);
  if(paramsBool["writeCorrzz"])
    write2DToFile (label + "_corrzz", corrzz);
  if(paramsBool["writeCorrxx"])
    write2DToFile (label + "_corrxx", corrxx);
  if(paramsBool["writeS2"])
    writeToFile (label + "_S2", S2, true);
  if(paramsBool["writeMaxDim"])
    writeToFile (label + "_maxDim", chi, true);
  if(paramsBool["writeTrace"])
    writeToFile (label + "_trace", traces, true);


   
  // //Some tests
  // if (vectorized)
  //   dmt.unvec();
  // //auto overlap = innerC(psi,psi0);
  // //Print(std::norm(overlap));
  // //Print(traceC(rho,swapPrime(projector(psi0),0,1)));
  // //auto initialE = innerC(psi0, H, psi0);
  // //Print(initialE);
  // //auto E = innerC(psi, H, psi);
  // //Print(E);
  // auto Erho = traceC(rho, H);
  // Print(Erho);

  return 0;
}
#endif

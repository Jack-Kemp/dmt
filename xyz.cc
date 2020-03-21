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
#include"TrotterConstructor.h"



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
					     "Cutoff",
					     "FirstSVDCutoff",
					     "ThirdSVDCutoff",
					     "Jz",
					     "Jx",
					     "Jy",
					     "J2z",
					     "J2y",
					     "J2x",
					     "hz",
					     "hx",
					     "hy"};
  
  std::vector<std::string> paramNamesInt = {"N", 
					    "MaxDim",
					    "PresRadius",
					    "nSweeps"};

  std::vector<std::string> paramNamesBool = {"Vectorize",
					     "HermitianBasis",
					     "FromPureState",
					     "NextNearest",
					     "ConserveQNs",
					     "Normalize",
					     "OnlyPreserveEnergyDensity",
					     "UseSVD",
					     "UseSVDThird",
					     "AbsoluteCutoff",
					     "AbsolutePresCutoff",
					     "Verbose",
					     "WriteMagz",
					     "WriteMagx",
					     "WriteMagy",
					     "WriteCorrzz",
					     "WriteCorrxx",
					     "WriteCorryy",
					     "WriteS2",
					     "WriteMaxDim",
					     "WriteTrace"
					    };
  
  for(auto const& x : paramNamesReal)
    paramsReal[x] = input.getReal(x);
  for(auto const& x : paramNamesInt)
    paramsInt[x] = input.getInt(x);
  for(auto const& x : paramNamesBool)
    paramsBool[x] = input.getYesNo(x);

  const int N = paramsInt.at("N");
  const Real tStep = paramsReal.at("tStep");
  const Real tSweep = tStep/paramsInt.at("nSweeps");
  const bool vectorize = paramsBool.at("Vectorize");
  const bool hermitianBasis = paramsBool.at("HermitianBasis");
  const bool conserveQNs = paramsBool.at("ConserveQNs");

  Real hx = paramsReal.at("hx");
  Real hy = paramsReal.at("hy");
  Real hz = paramsReal.at("hz");
  Real Jx = paramsReal.at("Jx");
  Real Jy = paramsReal.at("Jy");
  Real Jz = paramsReal.at("Jz");
  Real J2x = paramsReal.at("J2x");
  Real J2y = paramsReal.at("J2y");
  Real J2z = paramsReal.at("J2z");
  Real delta = 0.5*Jx;

  auto svdMethod = input.getString("SVDMethod");
  
  auto outputDir = input.getString("OutputDir"); 
  auto outputName = input.getString("OutputName");

  //Sanity check conservation
  if (conserveQNs)
    if( hx != 0 or hy != 0 or (Jx != Jy))
      Error("Hamiltionian does not conserve Z. Check hx,hy = 0 and Jx = Jy");

  //End input-output options -----------------------------------------

 
  //Set up physical basis and initial state---------------------------

  auto sites = SpinHalf(N, {"ConserveQNs=", conserveQNs});
  auto vectorBasis = {"Id", "Sx", "Sy", "Sz"};
  
  auto dmt = DMT(sites, vectorBasis, {"PresRadius", paramsInt.at("PresRadius"),
			       "HermitianBasis", hermitianBasis,
			       "Vectorize", vectorize});
  
  //From projector of wavefunction
  if(paramsBool.at("FromPureState"))
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
			   
  //Set up Hamilitonian----------------------------------------------


  TrotterConstructor trott;
  
  trott.addNearestNeighbour("Sz", "Sz", Jz);
  trott.addSingleSite("Sz", hz);
  
  if(conserveQNs)
    {
      trott.addNearestNeighbour("S+","S-", delta);
      trott.addNearestNeighbour("S-","S+", delta);
    }
  else
    {
      trott.addSingleSite("Sx", hx);
      trott.addSingleSite("Sy", hy);
      trott.addNearestNeighbour("Sx","Sx", Jx);
      trott.addNearestNeighbour("Sy","Sy", Jy);
    }

  if(paramsBool.at("NextNearest"))
    {
      trott.addLongRange("Sz", "Sz", 2, J2z);
      trott.addLongRange("Sy", "Sy", 2, J2y);
      trott.addLongRange("Sx", "Sx", 2, J2x);
      
    }
    
  
  auto dmtgates = trott.twoSiteGates2ndOrderSweep(dmt, sites, tSweep, {"Verbose", true});

  //Set preserved operators and finish DMT set-up.

  if(paramsBool.at("OnlyPreserveEnergyDensity"))
    dmt.addPresOperator(trott.localEnergyDensity(sites), trott.maxRange()+1, true);

  dmt.finishConstruction();
  

  //Set up measurements-----------------------------------------------
  std::vector<std::vector<Real>> magz, magx,magy, corrzz, corrxx, corryy;
  std::vector<Real> S2, chi, traces;

  auto measure = [&](DMT& dmt, Args const & args){
		   magz.push_back(std::vector<Real>());
		   magx.push_back(std::vector<Real>());
		   magy.push_back(std::vector<Real>());
		   corrzz.push_back(std::vector<Real>());
		   corrxx.push_back(std::vector<Real>());
		   corryy.push_back(std::vector<Real>());
		   for(int i = 1; i <= N; i++)
		     {
		       if(paramsBool.at("WriteMagz"))
			 magz.back().push_back(calculateExpectation("Sz", i, dmt).real());
		       if(paramsBool.at("WriteMagx"))
			 magx.back().push_back(calculateExpectation("Sx", i, dmt).real());
		       if(paramsBool.at("WriteMagy"))
			 magy.back().push_back(calculateExpectation("Sy", i, dmt).real());
		       if(paramsBool.at("WriteCorrzz"))
			 corrzz.back().push_back(calculateTwoPoint("Sz", i, "Sz", (i%N)+1, dmt).real());
		       if(paramsBool.at("WriteCorrxx"))
			 corrxx.back().push_back(calculateTwoPoint("Sx", i, "Sx", (i%N)+1, dmt).real());
		       if(paramsBool.at("WriteCorryy"))
			 corryy.back().push_back(calculateTwoPoint("Sy", i, "Sy", (i%N)+1, dmt).real());
		     }
		   if(paramsBool.at("WriteS2"))
		     S2.push_back(secondRenyiEntropyHalfSystem(dmt));
		   if(paramsBool.at("WriteMaxDim"))
		     chi.push_back(maxLinkDim(dmt.rho()));
		   if(paramsBool.at("WriteTrace"))
		     traces.push_back(dmt.trace()); 
		 };

  DMTObserver obs(measure);

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  //Time evolve-------------------------------------------------------
  gateTEvol(dmtgates,paramsReal.at("tTotal"),tStep,dmt, obs,
	    {"Cutoff",paramsReal.at("Cutoff"),
	     "FirstSVDCutoff", paramsReal.at("FirstSVDCutoff"),
	     "ThirdSVDCutoff", paramsReal.at("ThirdSVDCutoff"),
	     "AbsoluteCutoff", paramsBool.at("AbsoluteCutoff"),
	     "AbsolutePresCutoff", paramsBool.at("AbsolutePresCutoff"),
	     "Verbose", paramsBool.at("Verbose"),
	     "MaxDim", paramsInt.at("MaxDim"),
	     "UseSVD", paramsBool.at("UseSVD"),
	     "UseSVDThird", paramsBool.at("UseSVDThird"),
	     "DoNormalize", paramsBool.at("Normalize"),
	     "nSweeps", paramsInt.at("nSweeps"),
	     "SVDMethod", svdMethod});

  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

  printfln("DMT evolution took %f seconds.", time_span.count());

  //gateTEvol(mpsgates,paramsReal.at("ttotal"),tstep,psi,{"Cutoff=",paramsReal.at("cutoff"),"Verbose=",true});
  

  //printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));
  printfln("Maximum MPO bond dimension after time evolution is %d",maxLinkDim(dmt.rho()));

  std::string label = outputDir + "/" + construct_label(outputName, "",
							paramsReal,
							paramsBool,
							paramsInt); 

  //Write results-----------------------------------------------------
  std::filesystem::create_directory(outputDir);
  if(paramsBool.at("WriteMagz"))
    write2DToFile (label + "_magz", magz);
  if(paramsBool.at("WriteMagx"))
    write2DToFile (label + "_magx", magx);
  if(paramsBool.at("WriteMagy"))
    write2DToFile (label + "_magy", magy);
  if(paramsBool.at("WriteCorrzz"))
    write2DToFile (label + "_corrzz", corrzz);
  if(paramsBool.at("WriteCorrxx"))
    write2DToFile (label + "_corrxx", corrxx);
  if(paramsBool.at("WriteCorryy"))
    write2DToFile (label + "_corryy", corryy);
  if(paramsBool.at("WriteS2"))
    writeToFile (label + "_S2", S2, true);
  if(paramsBool.at("WriteMaxDim"))
    writeToFile (label + "_maxDim", chi, true);
  if(paramsBool.at("WriteTrace"))
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

#if 1
#include "itensor/all.h"
#include <itensor/mps/bondgate.h>
#include<chrono>
#include<limits>
#include<functional>
#include<filesystem>
#include "allDMT.h"

using namespace itensor;
using std::vector;
using namespace std::chrono;

typedef std::vector<Real> VecReal;
typedef std::vector<VecReal> MatrixReal;
typedef std::vector<std::string> VecStr;


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
  auto args = readInMandatoryArgs(input);

  const int N = input.getInt("N");
  const int centre = N %2 == 0 ?  N/2 : (N+1)/2;
  
  const Real tStep = input.getReal("tStep");
  const Real tSweep = tStep/input.getInt("nSweeps");
  const Real tTotal = input.getReal("tTotal");

  Real hx = input.getReal("hx");
  Real hy = input.getReal("hy");
  Real hz = input.getReal("hz");
  Real Jx = input.getReal("Jx");
  Real Jy = input.getReal("Jy");
  Real Jz = input.getReal("Jz");
  Real J2x = input.getReal("J2x");
  Real J2y = input.getReal("J2y");
  Real J2z = input.getReal("J2z");
  Real delta = 0.5*Jx;

  const bool conserveQNs = input.getYesNo("ConserveQNs");
  
  auto outputDir = input.getString("OutputDir"); 
  auto outputName = input.getString("OutputName");

  //Sanity check conservation
  if (conserveQNs)
    if( hx != 0 or hy != 0 or (Jx != Jy))
      Error("Hamiltionian does not conserve Z. Check hx,hy = 0 and Jx = Jy");

  //End input-output options -----------------------------------------

 
  //Set up physical basis and initial state---------------------------

  auto sites = SpinHalf(N, {"ConserveQNs=", conserveQNs});
  std::vector<std::string> vectorBasis = {"Id", "Sx", "Sy", "Sz"};
  
  auto dmt = DMT(sites, vectorBasis, args);
  
  //From projector of wavefunction
  if(input.getYesNo("FromPureState"))
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
	  //Alternating
	  dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*(2*(j %2)-1)*dmt.stateOp("Sz",j);

	  //Unpolarized except for centre spin
	  //dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*(j == centre)*dmt.stateOp("Sz",j);

	  //Single Domain Wall
	  //dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*( 2*(j > centre) -1 )*dmt.stateOp("Sz",j);
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

  if(input.getYesNo("NextNearest"))
    {
      trott.addLongRange("Sz", "Sz", 2, J2z);
      trott.addLongRange("Sy", "Sy", 2, J2y);
      trott.addLongRange("Sx", "Sx", 2, J2x);      
    }
    
  
  auto dmtgates = trott.twoSiteGates2ndOrderSweep(dmt, sites, tSweep, {"Verbose", true});
  auto hamiltonian = trott.hamiltonian(sites);
  auto maxRange = trott.maxRange();

  std::vector<ITensor> localEnergyDensity, HSzComm;
  if (input.getYesNo("WriteEnergyDensity", false))
    for(int i=1; i<=N; ++i)
      localEnergyDensity.push_back(trott.localEnergyDensity(i, dmt.dmtSites()));

   if (input.getYesNo("WriteHSzComm", false))
    for(int i=1; i<=N; ++i)
      HSzComm.push_back(dmt.convertToSiteOp(trott.commuteWithSingleSite(i, "Sz", sites), i- maxRange, i+maxRange)); 

  //Set preserved operators and finish DMT set-up.

  if(input.getYesNo("OnlyPreserveEnergyDensity"))
    dmt.addPresOperator(trott.localEnergyDensity(sites), trott.maxRange()+1, true);

  dmt.finishConstruction();
  

  //Set up measurements-----------------------------------------------

  std::map<std::string, MatrixReal> data2D;
  std::map<std::string, VecReal> data;


  VecStr dataNames2D = {
		  "Sz",
		  "Sx",
		  "Sy",
		  "SzSzNN",
		  "SxSxNN",
		  "SySyNN",
		  "EnergyDensity",
		  "HSzComm",
		  "SpinCurrent"
  };

  VecStr dataNames = {"S2",
		  "Trace",
		  "Energy",
		  "MaxDim",
		  "TruncErr"
  };

  for (const auto & name : dataNames2D)
    if(input.getYesNo("Write" + name, false))
      data2D.emplace(name, MatrixReal());
  for (const auto & name : dataNames)
    if(input.getYesNo("Write" + name, false))
      data.emplace(name, VecReal());
  data.emplace("t", VecReal());

  auto measure = [&](DMT& dmt, Args const & args){
		   for (auto & [key, value] : data2D)
		     value.push_back(VecReal());
		   for(int i = 1; i <= N; i++)
		     {
		       for (auto & [key, value] : data2D)
			 {
			   Real ret = std::numeric_limits<Real>::quiet_NaN();
			   switch( hash(key.c_str()) ){
			   case "Sz"_: ret  = calculateExpectation("Sz", i, dmt).real(); break;
			   case "Sx"_: ret  = calculateExpectation("Sx", i, dmt).real(); break;
			   case "Sy"_: ret  = calculateExpectation("Sy", i, dmt).real(); break;
			   case "SzSzNN"_: ret  = calculateTwoPoint("Sz", i, "Sz", (i%N)+1, dmt).real(); break;
			   case "SxSxNN"_: ret  = calculateTwoPoint("Sx", i, "Sx", (i%N)+1, dmt).real(); break;
			   case "SySyNN"_: ret  = calculateTwoPoint("Sy", i, "Sy", (i%N)+1, dmt).real(); break;
			   case "EnergyDensity"_: ret  = calculateExpectation(localEnergyDensity[i-1], i-maxRange, i+maxRange, dmt).real(); break;
			   case "HSzComm"_: ret  = calculateExpectation(HSzComm[i-1], i-maxRange, i+maxRange, dmt).real(); break;
			   case "SpinCurrent"_: ret =  calculateTwoPoint("Sx", i, "Sy", (i%N)+1, dmt).real() - calculateTwoPoint("Sy", i, "Sx", (i%N)+1, dmt).real(); break;
			   }
			   value.back().push_back(ret);
			 }
		     }
		   for (auto & [key, value] : data)
		     {
		       Real ret = std::numeric_limits<Real>::quiet_NaN();
		       switch( hash(key.c_str()) ){
		       case "S2"_: ret  = secondRenyiEntropyHalfSystem(dmt); break;
		       case "MaxDim"_: ret = maxLinkDim(dmt.rho()); break;
		       case "Trace"_: ret = dmt.trace(); break;
		       case "Energy"_: ret = calculateExpectation(hamiltonian, dmt).real(); break;
		       case "t"_: ret = args.getReal("Time"); break;
		       case "TruncErr"_: ret = args.getReal("TruncError"); break;
		       }
		       value.push_back(ret);
		     }
		 };

  
  DMTObserver obs(measure);

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  //Time evolve-------------------------------------------------------
  gateTEvol(dmtgates, tTotal, tStep, dmt, obs, args);

  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);

  printfln("DMT evolution took %f seconds.", timeSpan.count());

  //gateTEvol(mpsgates,paramsReal.at("ttotal"),tstep,psi,{"Cutoff=",paramsReal.at("cutoff"),"Verbose=",true});
  

  //printfln("Maximum MPS bond dimension after time evolution is %d",maxLinkDim(psi));
  printfln("Maximum MPO bond dimension after time evolution is %d",maxLinkDim(dmt.rho()));

  std::map<std::string, double> runInfo = {{"TimeTaken", timeSpan.count()},
					   {"MaxBondDim", maxLinkDim(dmt.rho())}};

  std::string label = outputDir + "/" +outputName; 

  //Write results-----------------------------------------------------
  std::filesystem::create_directory(outputDir);

  writeDataToFile(label, data, data2D, runInfo, inputName);

  return 0;
}
#endif

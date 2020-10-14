#if 1
#include "allDMT.h"
#include "itensor/all.h"
#include <itensor/mps/bondgate.h>
#include<chrono>
#include<limits>
#include<functional>
#include<filesystem>




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
    inputName = "xyzLadder_example_input_params.txt";
  }
  
  auto input = InputGroup(inputName,"input");
  auto args = readInMandatoryArgs(input);

  //***Code for ladder assumes L is a multiple of 2***
  const int L = input.getInt("L");
  const int N = 2*L;
  const int centre = L;
  
  const Real tStep = input.getReal("tStep");
  const Real tSweep = tStep/input.getInt("nSweeps");
  const Real tTotal = input.getReal("tTotal");


  //Lattice and couplings:
  //      Jl
  // --1------3-----5-----
  //   |      |     |
  // Jr|      |     |
  //   |      |     |
  // --2------4-----6-----
	 
    
  Real Jlx = input.getReal("Jlx");
  Real Jly = input.getReal("Jly");
  Real Jlz = input.getReal("Jlz");
  Real Jrx = input.getReal("Jrx");
  Real Jry = input.getReal("Jry");
  Real Jrz = input.getReal("Jrz");
  Real eta = input.getReal("eta");
  Real rlRatio = input.getReal("rlRatio");
  
  auto outputDir = input.getString("OutputDir"); 
  auto outputName = input.getString("OutputName");

  //End input-output options -----------------------------------------

 
  //Set up physical basis and initial state---------------------------

  auto sites = SpinHalf(N, {"ConserveQNs=", false});
  VecStr vectorBasis = {"Id", "Sx", "Sy", "Sz"};
  
  auto dmt = DMT(sites, vectorBasis, args);
  
  for(int j = 1; j <= N; ++j)
	{
	  //Single Domain Wall
	  dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*eta*( 2*(j <= centre) -1 )*dmt.stateOp("Sz",j);
	}

  dmt.finishConstruction();
  
  //Set up Hamilitonian----------------------------------------------


  TrotterConstructor trott;
  
  trott.addNearestNeighbour("Sx","Sx", VecReal{Jrx*rlRatio, 0.}, {"Repeat", true});
  trott.addNearestNeighbour("Sy","Sy", VecReal{Jry*rlRatio, 0.}, {"Repeat", true});
  trott.addNearestNeighbour("Sz", "Sz", VecReal{Jrz*rlRatio, 0.}, {"Repeat", true});
  trott.addLongRange("Sx", "Sx", 2, Jlx);
  trott.addLongRange("Sy", "Sy", 2, Jly);
  trott.addLongRange("Sz", "Sz", 2, Jlz);
    
  
  auto dmtgates = trott.twoSiteGates2ndOrderSweep(dmt, sites, tSweep, {"Verbose", true});
  auto hamiltonian = trott.hamiltonian(sites);
  auto maxRange = trott.maxRange();

  
  

  //Set up measurements-----------------------------------------------

  std::map<std::string, MatrixReal> data2D;
  std::map<std::string, VecReal> data;

   VecStr dataNames2D = {
		  "Sz",
		  "Sx",
		  "Sy",
		  "SxMidSx",
		  "SyMidSy",
		  "SzMidSz",
		  "EnergyDensity"
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

  std::vector<ITensor> localEnergyDensity;
  if (input.getYesNo("WriteEnergyDensity", false))
    for(int i=1; i<=N; ++i)
      localEnergyDensity.push_back(trott.localEnergyDensity(i, dmt.dmtSites()));



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
			   case "SzMidSz"_: ret  = calculateTwoPoint("Sz", centre, "Sz", i, dmt).real(); break;
			   case "SxMidSx"_: ret  = calculateTwoPoint("Sx", centre, "Sx", i, dmt).real(); break;
			   case "SyMidSy"_: ret  = calculateTwoPoint("Sy", centre, "Sy", i, dmt).real(); break;
			   case "EnergyDensity"_: ret  = calculateExpectation(localEnergyDensity[i-1], i-maxRange, i+maxRange, dmt).real(); break;  
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
  printfln("Maximum MPO bond dimension after time evolution is %d",maxLinkDim(dmt.rho()));



  //Write results-----------------------------------------------------
  std::map<std::string, double> runInfo = {{"TimeTaken", timeSpan.count()},
					   {"MaxBondDim", maxLinkDim(dmt.rho())}};

  std::string label = outputDir + "/" +outputName; 

  
  std::filesystem::create_directory(outputDir);

  writeDataToFile(label, data, data2D, runInfo, inputName);

  return 0;
}
#endif

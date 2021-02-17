#if 1
#include "allDMT.h"
#include "itensor/all.h"

#include<chrono>
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
    inputName = "xyzAdvanced_input_params.txt";
  }
  
  auto input = InputGroup(inputName,"input");
  auto args = readInMandatoryArgs(input);

  const int N = input.getInt("N");
  const int centre = N %2 == 0 ?  N/2 : (N+1)/2;
  
  const Real tStep = args.getReal("tStep");
  const Real tSweep = tStep/args.getInt("nSweeps");
  const Real tTotal = args.getReal("tTotal");
  Real tStart = args.getReal("tStart");

  Real hx = input.getReal("hx");
  Real hy = input.getReal("hy");
  Real hz = input.getReal("hz");
  Real Jx = input.getReal("Jx");
  Real Jy = input.getReal("Jy");
  Real Jz = input.getReal("Jz");
  Real J2x = input.getReal("J2x");
  Real J2y = input.getReal("J2y");
  Real J2z = input.getReal("J2z");
  
  auto outputDir = args.getString("OutputDir");
  std::filesystem::create_directory(outputDir);
  auto outputName = args.getString("OutputName");
  auto checkpointName = args.getString("CheckpointName");

  //End input-output options -----------------------------------------

 
  //Set up physical basis and initial state---------------------------
  
  auto sites = SpinHalf(N, args);
  std::vector<std::string> vectorBasis = {"Id", "Sx", "Sy", "Sz"};

  std::string checkDMTName;
  bool loadFromCheckpoint = false;
  if(args.getBool("Checkpoint") and
     findCheckpointFile(checkpointName, outputDir, ".dmt", tStart, checkDMTName)
     )
    {
      loadFromCheckpoint = true;
      args.add("tStart", tStart);
      printfln("Found checkpoint file %s.", checkDMTName);
      readFromFile(checkDMTName + ".sites", sites);
    }

  auto dmt = DMT(sites, vectorBasis, args);

  if(loadFromCheckpoint)
    {
      dmt.readFromFile(checkDMTName);
    }
  else
    {
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
	      //dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*(2*(j %2)-1)*dmt.stateOp("Sz",j);
	      
	      //Unpolarized except for centre spin
	      //dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*(j == centre)*dmt.stateOp("Sz",j);

	      //Single Domain Wall
	      dmt.rhoRef(j) = dmt.stateOp("Id",j)+4*0.1*( 2*(j <= centre) -1 )*dmt.stateOp("Sz",j);
	    }
	}
    }
			   
  //Set up Hamilitonian----------------------------------------------


  TrotterConstructor trott;
  
  trott.addNearestNeighbour("Sz", "Sz", Jz);
  trott.addSingleSite("Sz", hz);
  
  trott.addSingleSite("Sx", hx);
  trott.addSingleSite("Sy", hy);
  trott.addNearestNeighbour("Sx","Sx", Jx);
  trott.addNearestNeighbour("Sy","Sy", Jy);

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


  

  //Set up checkpointing-----------------------------------------------------

  auto write_checkpoint = [&](DMT& dmt, Args const & args){
			    std::map<std::string, double> runInfo = {{"TimeTaken", args.getReal("WallTime")},
								     {"MaxBondDim", maxLinkDim(dmt.rho())}};
			    dmt.writeToFile(args.getString("WriteFilename") + ".dmt");
			    writeDataToFile(args.getString("WriteFilename") + ".dat", data, data2D, runInfo, inputName);
			  };

  
  DMTObserver obs(measure, write_checkpoint, args);

  
  //Time evolve-------------------------------------------------------
  gateTEvol(dmtgates, tTotal-tStart, tStep, dmt, obs, args);

  //Write Final Results
  std::string label = outputDir + "/" + outputName;

  write_checkpoint(dmt, {"WriteFilename=", label, "WallTime=", obs.getRunTime()});

  return 0;
}
#endif

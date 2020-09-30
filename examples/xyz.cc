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

int main(int argc, char* argv[])
{

  //Get input parameters and options --------------------------------
  const char* inputName;
  if (argc > 1) {
    inputName = argv[1];
  }
  else{ 
    printfln("Using default parameter file.");
    inputName = "xyz_example_input_params.txt";
  }
  
  auto input = InputGroup(inputName,"input");
  auto args = readInMandatoryArgs(input);

  const int N = input.getInt("N");
  const Real tStep = input.getReal("tStep");
  const Real tSweep = tStep/input.getInt("nSweeps");
  const Real tTotal = input.getReal("tTotal");

  Real hx = input.getReal("hx");
  Real hy = input.getReal("hy");
  Real hz = input.getReal("hz");
  Real Jx = input.getReal("Jx");
  Real Jy = input.getReal("Jy");
  Real Jz = input.getReal("Jz");

  auto outputDir = input.getString("OutputDir"); 
  auto outputName = input.getString("OutputName");

  //End input-output options -----------------------------------------

 
  //Set up physical basis and initial state---------------------------

  auto sites = SpinHalf(N, {"ConserveQNs=", false});
  auto vectorBasis = {"Id", "Sx", "Sy", "Sz"};
  
  auto dmt = DMT(sites, vectorBasis, args);
  
  //Initialise from product of MPO operators
  for(int j = 1; j <= N; ++j)
	{
	  //Alternating Spins
	  dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*(2*(j %2)-1)*dmt.stateOp("Sz",j);

	  //Unpolarized except for centre spin
	  //dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*(j == centre)*dmt.stateOp("Sz",j);

	  //Single Domain Wall
	  //dmt.rhoRef(j) = dmt.stateOp("Id",j)+2*( 2*(j > centre) -1 )*dmt.stateOp("Sz",j);
	}
			   
  //Set up Hamilitonian----------------------------------------------


  TrotterConstructor trott;

  trott.addSingleSite("Sx", hx);
  trott.addSingleSite("Sy", hy);
  trott.addSingleSite("Sz", hz);
 
  trott.addNearestNeighbour("Sx","Sx", Jx);
  trott.addNearestNeighbour("Sy","Sy", Jy);
  trott.addNearestNeighbour("Sz", "Sz", Jz);
        
  auto dmtgates = trott.twoSiteGates2ndOrderSweep(dmt, sites, tSweep, args);

  dmt.finishConstruction();
  

  //Set up measurements-----------------------------------------------

  std::map<std::string, MatrixReal> data2D;
  std::map<std::string, VecReal> data;
  
  data2D.emplace("Sz", MatrixReal());
  data2D.emplace("SzSzNN", MatrixReal());

  data.emplace("S2", VecReal());
  data.emplace("MaxDim", VecReal());

  auto measure = [&](DMT& dmt, Args const & args){
		   for (auto & [key, value] : data2D)
		     value.push_back(VecReal());
		   for(int i = 1; i <= N; i++)
		     {
			 data2D["Sz"].back().push_back(calculateExpectation("Sz", i, dmt).real());
			 data2D["SzSzNN"].back().push_back(calculateTwoPoint("Sz", i, "Sz", (i%N)+1, dmt).real());
		     }
		   data["S2"].push_back(secondRenyiEntropyHalfSystem(dmt));
		   data["MaxDim"].push_back(maxLinkDim(dmt.rho()));
		 };

  
  DMTObserver obs(measure);

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  //Time evolve-------------------------------------------------------
  gateTEvol(dmtgates, tTotal, tStep, dmt, obs, args);

  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  duration<double> timeSpan = duration_cast<duration<double>>(t2 - t1);

  printfln("DMT evolution took %f seconds.", timeSpan.count());
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

#ifndef DEFAULTARGS_H
#define DEFAULTARGS_H

#include<itensor/all.h>

#include<map>
using namespace itensor;


//Read the mandatory arguments from input, using defaults
//if not set for arguments with sensible defaults.
Args readInMandatoryArgs(InputGroup & input){

  Args args;


  std::vector<std::string> mandParamsReal = {"tStep", "tTotal"};
  std::vector<std::string> mandParamsInt = {"nSweeps", "MaxDim"};
  
  
  
  
  std::map<std::string, Real> defParamsReal = { {"Cutoff", 1e-16},
  					     {"FirstSVDCutoff", 1e-16},
  					     {"ThirdSVDCutoff", 1e-10}
  					     };
  
  std::map<std::string, int> defParamsInt = { {"PresRadius", 1} };

  std::map<std::string, bool> defParamsBool = {
  					    {"Vectorize", true},
  					    { "HermitianBasis", true},
  					    { "ConserveQNs", false},
  					    { "DoNormalize", true},
					    { "CacheTrace", true},
  					    { "OnlyPreserveEnergyDensity", false},
  					    { "UseSVD", true},
  					    { "UseSVDThird", true},
  					    { "AbsoluteCutoff", true},
  					    { "AbsolutePresCutoff", true},
  					    { "Verbose", true},
					    { "FromPureState", false},
					    { "ConserveQNs", false} 
  };

  std::map<std::string, std::string> defParamsString = {{"SVDMethod", "gesdd"}};


  for(auto const& k : mandParamsReal)
    args.add(k, input.getReal(k));
  for(auto const& k : mandParamsInt)
    args.add(k, input.getInt(k));
  
  for(auto const& [k,v] : defParamsReal)
    args.add(k, input.getReal(k,v));
  for(auto const& [k,v] : defParamsInt)
    args.add(k, input.getInt(k,v));
  for(auto const& [k,v] : defParamsBool)
    args.add(k, input.getYesNo(k,v));
  for(auto const& [k,v] : defParamsString)
    args.add(k, input.getString(k,v));

   return args;

 }
#endif

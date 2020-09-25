#ifndef DEFAULTARGS_H
#define DEFAULTARGS_H

#include<itensor/all.h>

#include<map>
using namespace itensor;


//Read the mandatory arguments from input, using defaults from
//DefaultArgs.h if not set.
Args readInMandatoryArgs(InputGroup & input){

  Args args;
  
  std::map<std::string, Real> paramsReal = { {"Cutoff", 1e-16},
  					     {"FirstSVDCutoff", 1e-16},
  					     {"ThirdSVDCutoff", 1e-10}
  					     };
  
  std::map<std::string, int> paramsInt = { {"MaxDim", 64},
  					   {"PresRadius", 1} };

  std::map<std::string, bool> paramsBool = {
  					    {"Vectorize", true},
  					    { "HermitianBasis", true},
  					    { "ConserveQNs", false},
  					    { "DoNormalize", true},
  					    { "OnlyPreserveEnergyDensity", false},
  					    { "UseSVD", true},
  					    { "UseSVDThird", true},
  					    { "AbsoluteCutoff", true},
  					    { "AbsolutePresCutoff", true},
  					    { "Verbose", true},
					    { "FromPureState", false},
					    { "ConserveQNs", false} 
  };

  std::map<std::string, std::string> paramsString = {{"SVDMethod", "gesdd"}};
  
  for(auto const& [k,v] : paramsReal)
    args.add(k, input.getReal(k,v));
  for(auto const& [k,v] : paramsInt)
    args.add(k, input.getInt(k,v));
  for(auto const& [k,v] : paramsBool)
    args.add(k, input.getYesNo(k,v));
  for(auto const& [k,v] : paramsString)
    args.add(k, input.getString(k,v));

   return args;

 }
#endif

#ifndef __ITENSOR_TROTTERCONSTRUCTOR_H
#define __ITENSOR_TROTTERCONSTRUCTOR_H

#include<itensor/all.h>
#include<limits>
#include<map>
#include<string>

using Real = double;

//Optional args: bool repeat, int start, int end.
//repeat=true or single value: couplings repeat between start and end inclusive.
//Otherwise provide a vector of length either of the system size or (end-start+1),
//for the value on each site between start and end, or the entire system.
class CouplingValues {
  using VecReal= std::vector<Real>;
  VecReal vals_;
  bool repeat_;
  int len_;
  int start_, end_;

public:
  
  CouplingValues(Real val, Args args = Args::global()):
    vals_({val}),
    repeat_(true),
    len_(1),
    start_(args.getInt("Start", 1)),
    end_(args.getInt("End",std::numeric_limits<int>::max()))													
  {   
  }

  CouplingValues(VecReal vals, Args args = Args::global()):
    vals_(vals),
    repeat_(args.getBool("Repeat", false)),
    len_(vals.size()),
    start_(args.getInt("Start", 1)),
    end_(args.getInt("End",std::numeric_limits<int>::max()))
  {
  }

  CouplingValues(){}

  Real operator () (int i) const{
    if(i >= start_ and i <= end_)
      return repeat_ ? vals_[(i - start_) % len_] : vals_[i - start_];
    return 0.0;
  }

  bool isUniform() const {
    return repeat_ and len_==1 and start_ == 1
      and end_ == std::numeric_limits<int>::max();
  }
			    
   
};

 

class TrotterConstructor {
  using MapCouplings = std::map<std::vector<std::string>, CouplingValues>;
  MapCouplings nearestNeighbour_;
  MapCouplings singleSite_;
  
public:


  template<typename ValueType>
  void addNearestNeighbour(const std::string & opnameL,
			   const std::string & opnameR,
			   ValueType values,
			   const Args & args = Args::global())
  {
    nearestNeighbour_[{opnameL, opnameR}] = CouplingValues(values, args);
  }

  template<typename ValueType>
  void addSingleSite(const std::string & opnameL,
			   ValueType values,
			   const Args & args = Args::global())
  {
    singleSite_[{opnameL}] = CouplingValues (values, args);
  }

  template<typename CalcGateClass>
  std::vector<BondGate> twoSiteGates2ndOrderSweep(const CalcGateClass & calc,
						  SiteSet sites,
						  Real tSweep, const Args& args = Args::global()) {
    bool verbose = args.getBool("Verbose", false);
    const int N  = length(sites);
    std::vector<BondGate> gates;
    for(int b = 1; b < N; ++b)
    {
      auto hterm = ITensor(op(sites,"Id",1).inds());
      for (const auto& [opnames, cvals] : nearestNeighbour_){
	hterm += cvals(b)*op(sites, opnames[0], b)*op(sites, opnames[1], b+1);
	if (verbose) printfln((opnames[0] + opnames[1] + " at %d,%d: %f").c_str(), b, b+1,  cvals(b));
      }

      for (const auto& [opnames, cvals] : singleSite_)
	{
	hterm += 0.5*cvals(b)*op(sites, opnames[0], b)*op(sites, "Id", b+1);
	hterm += 0.5*cvals(b+1)*op(sites, "Id", b)*op(sites, opnames[0], b+1);
	if (verbose) printfln((opnames[0] + " at %d: %f").c_str(), b, cvals(b));
	}
      
      if(b == 1){
	for (const auto& [opnames, cvals] : singleSite_)
	  hterm += 0.5*cvals(b)*op(sites, opnames[0], b)*op(sites, "Id", b+1);
      }
      if(b == N-1){
	for (const auto& [opnames, cvals] : singleSite_){
	  hterm += 0.5*cvals(b+1)*op(sites, "Id", b)*op(sites, opnames[0], b+1);
	  if (verbose) printfln((opnames[0] + " at %d: %f").c_str(), b+1,  cvals(b+1));
	}
      }
      
      gates.push_back(calc.calcGate(hterm, tSweep/2, b, args)); //Notice over 2!!!!
    }
    
  for(int b = N-1; b >= 1; --b)
    {
      gates.push_back(gates[b-1]);
    }

  return gates;

  }






};

#endif

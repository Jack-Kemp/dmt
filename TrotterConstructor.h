#ifndef __ITENSOR_TROTTERCONSTRUCTOR_H
#define __ITENSOR_TROTTERCONSTRUCTOR_H

#include<itensor/all.h>
#include<limits>
#include<map>
#include<string>

#include"ITensorUtilFunctions.h"


using Real = double;
namespace itensor{



  //Class CouplingValues to store information above the couplings.
  //Can either by single value or full list or list to be periodically
  //repeated. Helper for TrotterConstructor below to store and
  //construct Hamiltonians.

  //bool Repeat, int Start, int End.

  //If Repeat=true or single value is given: couplings repeat between site numbered by start
  //and site numbered by end inclusive.
  
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

  CouplingValues(const VecReal & vals, Args args = Args::global()):
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


BondGate
swapInputOutput(const BondGate & gate, const SiteSet & sites)
{
  if(order(findInds(gate.gate(), "2")) > 0)
    return BondGate(sites, gate.i1(), gate.i2(), swapPrime(swapPrime(gate.gate(),0,2),1,3));
  else
    return BondGate(sites, gate.i1(), gate.i2(), swapPrime(gate.gate(),0,1));
}



  //Class to store Coupling Values with lattice bonds, and construct
  //appropriate two site gates.

  //Use addSingleSite, addNearestNeighbour, addLongRange to add couplings
  //the Hamiltonian.

  //Then use construct2ndOrderSweep to get a list of gates which when applied
  //sweep once back and forth across the chain for a 2nd order Trotter decomposition
  // of the Hamiltonian.
  
class TrotterConstructor {
  using MapCouplings = std::map<std::vector<std::string>, CouplingValues>;
  MapCouplings nearestNeighbour_;
  std::map<int, MapCouplings> longRange_;
  MapCouplings singleSite_;
  int maxRange_ = 0;




  //Construct second order DMRG-like sweep of two-site gates
  //back and forth once. Results in 2nd order trotter decomposition.
  //Assumes open system and deals with the edge appropriately. 
  template<typename CalcGateClass>
  void 
  construct2ndOrderSweepNearestNeighbour_(std::vector<BondGate>& gates,
					  const CalcGateClass & calc,
					  SiteSet& sites,
					  Real tSweep,
					  const Args & args){
  bool verbose = args.getBool("Verbose", false);
  const int N  = length(sites);

  for(int b = 1; b < N; ++b)
    {
      auto hterm = ITensor(getId(sites,b,b+2).inds()).fill(0.0);
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
}


  //Constructs 2nd order sweep of two sites gates for arbitrary range interactions using swap gates.
  
  template<typename CalcGateClass>
  void 
  construct2ndOrderSweepLongRange_(std::vector<BondGate>& gates,
			      const CalcGateClass & calc,
			      SiteSet& sites,
			      Real tSweep,
			      const Args & argsIn) {

    Args args = argsIn;
    auto criter =  longRange_.rbegin();
    int range = criter->first;
    bool verbose = args.getBool("Verbose", false);
    const int N  = length(sites);
     

    for(int b = 1; b < N; ++b)
      {
	while (b + range > N) {
	  ++criter;
	  if (criter != longRange_.rend())
	    range = criter->first;
	  else
	    range = 1;
	}

	//First deal with nearest neighbour and single-site terms.
	auto hterm = ITensor(getId(sites,b,b+2).inds()).fill(0.0);
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
	
	//If we are next to the right edge, we are done, otherwise, long-range terms.
	if (range == 1)
	  {
	    args.add("SwapOutputs", false);
	    gates.push_back(calc.calcGate(hterm, tSweep/2, b, args)); //Notice over 2!!!!
	  }

	//For each seperation, construct the two-site gate, then swap again after
	//so we are ready to deal with the next seperation. Then apply the gates in
	//reverse order, swapping back. As each gate is applied twice, deltat =  tsweep/4.
	//Notice for there is no swapping after the longest-range gate which is only applied once.
	else
	  {
	    //Construct a swap gate after nearest neighbour terms.
	    int newgates = 0;
	    args.add("SwapOutputs", true);
	    gates.push_back(calc.calcGate(hterm, tSweep/4, b, args)); //Notice over 2 twice!!!!
	    newgates++;   
	    for (int sep = 2; sep < range; sep++){
	      auto couplings = longRange_.find(sep);
	      auto hterm =  ITensor(getId(sites,b+sep-1,b+sep+1).inds()).fill(0.0);
	      if (couplings != longRange_.end())
		for (const auto& [opnames, cvals] : couplings->second)
		  {
		    hterm += cvals(b)*op(sites, opnames[0], b+sep-1)*op(sites, opnames[1], b+sep);
		    if (verbose) printfln((opnames[0] + opnames[1] + " at %d,%d: %f").c_str(), b, b+sep,  cvals(b));
		  }
	      gates.push_back(calc.calcGate(hterm, tSweep/4, b+sep-1, args));
	      newgates++;
	    }
	    auto couplings = longRange_.find(range);
	    auto hterm =  ITensor(getId(sites,b+range-1,b+range+1).inds()).fill(0.0);
	    for (const auto& [opnames, cvals] : couplings->second)
	      {
		hterm += cvals(b)*op(sites, opnames[0], b+range-1)*op(sites, opnames[1], b+range);
		if (verbose) printfln((opnames[0] + opnames[1] + " at %d,%d: %f").c_str(), b, b+range,  cvals(b));
	      }
	    args.add("SwapOutputs", false);
	    gates.push_back(calc.calcGate(hterm, tSweep/2, b+range-1, args));
	    //Notice last gate added is not reversed and added back, but the rest are.
	    for(int i = 0; i < newgates; ++i)
	      {
		  gates.push_back(gates[gates.size()-2-i]);
	      }
	  }
      }
    for(int b = gates.size(); b > 0; --b)
      {
	gates.push_back(gates[b-1]);
      }
  }
  
public:


  int maxRange() const {return maxRange_;}
template<typename ValueType>
void
addLongRange(const std::string & opnameL,
		  const std::string & opnameR,
		  int latticeSeparation,
		  const ValueType & values,
		  const Args & args = Args::global())
{
  if (latticeSeparation == 1)
    Error("Please use addNearestNeighbour for separation = 1.");
  else if(latticeSeparation < 1)
    Error("Cannot have separation less than 1.");
  else
    {
    longRange_[latticeSeparation][{opnameL, opnameR}] = CouplingValues(values, args);
    maxRange_ = maxRange_ < (latticeSeparation) ? (latticeSeparation) : maxRange_;
    }
}


template<typename ValueType>
void
addNearestNeighbour(const std::string & opnameL,
			 const std::string & opnameR,
			 const ValueType & values,
			 const Args & args = Args::global())
{
  nearestNeighbour_[{opnameL, opnameR}] = CouplingValues(values, args);
  maxRange_ = maxRange_ < 1 ? 1 : maxRange_;
}

template<typename ValueType>
void
addSingleSite(const std::string & opnameL,
		   const ValueType & values,
		   const Args & args = Args::global())
{
  singleSite_[{opnameL}] = CouplingValues (values, args);
}

template<typename CalcGateClass>
std::vector<BondGate>
twoSiteGates2ndOrderSweep(const CalcGateClass & calc,
						SiteSet sites,
						Real tSweep, const Args& args = Args::global())
{
  std::vector<BondGate> gates;
  if(maxRange_ < 2)
    construct2ndOrderSweepNearestNeighbour_(gates,calc, sites, tSweep, args);
  else
    construct2ndOrderSweepLongRange_(gates, calc, sites, tSweep, args);
  // for (auto & gate : gates){
  //   PrintData(gate.i1());
  //   PrintData(gate.i2());
  // }
  
  return gates;


}

//Return the Hamiltonian as an MPO using ITensor's AutoMPO feature
//args controlling the precision of the MPO can be "Exact",
//and "CutOff" and "MaxDim".
  MPO
  hamiltonian(SiteSet sites, const Args & args = Args::global()){
    auto H = AutoMPO(sites);
    const int N  = length(sites);
    for (int b = 1; b <= N; ++b)
      {
	for (const auto& [opnames, cvals] : singleSite_)
	  {
	    H += cvals(b), opnames[0], b;
	  }
	if(b < N)
	  for (const auto& [opnames, cvals] : nearestNeighbour_)
	    {
	      H += cvals(b), opnames[0], b, opnames[1], b+1;
	    }
	for (int sep = 2; (sep <= maxRange_) and ((b+sep) <= N); sep++){
	  auto couplings = longRange_.find(sep);
	  if (couplings != longRange_.end())
	    for (const auto& [opnames, cvals] : couplings->second)
	      {
		 H += cvals(b), opnames[0], b, opnames[1], b+sep;
	      }
	}
      }
    return toMPO(H);
  }


  //Find the sum of the terms of the Hamiltonian soley
  //supported around site site_i, i.e site_i - maxRange to site_i + maxRange.
template<class Sites>
ITensor localEnergyDensity(int site_i, const Sites & sites){
  int start = std::max(site_i-maxRange_,1);
  int end = std::min(site_i+maxRange_, length(sites));

  auto hterm = ITensor(getId(sites, start,end+1).inds());
  for (int b = start; b <= end; ++b)
    {
      for (const auto& [opnames, cvals] : singleSite_)
	{
	  hterm += cvals(b)*getId(sites, start, b)*op(sites, opnames[0], b)*getId(sites, b+1, end+1);
	}
      if(b <= end - 1)
	{
	for (const auto& [opnames, cvals] : nearestNeighbour_)
	  {
	    hterm += cvals(b)*getId(sites, start, b)*op(sites, opnames[0], b)
	      *op(sites, opnames[1], b+1)*getId(sites, b+2, end+1);
	  }
	for (int sep = 2; sep <= std::min(end-b, maxRange_); sep++)
	  {
	  auto couplings = longRange_.find(sep);
	  if (couplings != longRange_.end())
	    for (const auto& [opnames, cvals] : couplings->second)
	      {
		hterm += cvals(b)*getId(sites, start, b)*op(sites, opnames[0], b)
		  *getId(sites, b+1, b+sep)*op(sites, opnames[1], b+sep)
		  *getId(sites, b+sep+1, end+1);
	      }
	  }
	}
    }
    return hterm;
  }



    //Find commutator [H, A_i] for operator A_i solely on site i. Cannot use DMTSiteSet
   //as requires operator products not traces
template<class Sites>
ITensor commuteWithSingleSite(int site_i, std::string commName, const Sites & sites){
  int start = std::max(site_i-maxRange_,1);
  int end = std::min(site_i+maxRange_, length(sites));
  auto commOp = op(sites, commName, site_i);
  auto hterm = ITensor(getId(sites, start,end+1).inds());
  for (int b = start; b <= end; ++b)
    {
      for (const auto& [opnames, cvals] : singleSite_)
	{
	  hterm += cvals(b)*getId(sites, start, b)*op(sites, opnames[0], b)*getId(sites, b+1, end+1);
	}
      if(b <= end - 1)
	{
	for (const auto& [opnames, cvals] : nearestNeighbour_)
	  {
	    hterm += cvals(b)*getId(sites, start, b)*op(sites, opnames[0], b)
	      *op(sites, opnames[1], b+1)*getId(sites, b+2, end+1);
	  }
	for (int sep = 2; sep <= std::min(end-b, maxRange_); sep++)
	  {
	  auto couplings = longRange_.find(sep);
	  if (couplings != longRange_.end())
	    for (const auto& [opnames, cvals] : couplings->second)
	      {
		hterm += cvals(b)*getId(sites, start, b)*op(sites, opnames[0], b)
		  *getId(sites, b+1, b+sep)*op(sites, opnames[1], b+sep)
		  *getId(sites, b+sep+1, end+1);
	      }
	  }
	}
    }
  auto acomm = multSiteOps(commOp, hterm);
  return 1_i*(dag(swapPrime(acomm,0,1)) - acomm);
  }
    

  //EXPERIMENTAL Find the sum of the terms of the Hamiltonian soley
  //supported within site 1 to site 1 + maxRange.
  ITensor localEnergyDensity(SiteSet sites){
    int maxSupp = maxRange_ + 1;
    auto hterm = ITensor(op(sites,"Id",1).inds());
    for (int b = 1; b <= maxSupp; ++b)
      {
	for (const auto& [opnames, cvals] : singleSite_)
	  {
	    hterm += cvals(b)*getId(sites, 1, b)*op(sites, opnames[0], b)*getId(sites, b+1, maxSupp+1);
	  }
	if(b <= maxSupp - 1)
	  for (const auto& [opnames, cvals] : nearestNeighbour_)
	    {
	      hterm += cvals(b)*getId(sites, 1, b)*op(sites, opnames[0], b)
		*op(sites, opnames[1], b+1)*getId(sites, b+2, maxSupp+1);
	    }
	for (int sep = 2; sep <= std::min(maxSupp-sep, maxRange_); sep++){
	  auto couplings = longRange_.find(sep);
	  if (couplings != longRange_.end())
	    for (const auto& [opnames, cvals] : couplings->second)
	      {
		hterm += cvals(b)*getId(sites, 1, b)*op(sites, opnames[0], b)
		  *getId(sites, b+1, b+sep)*op(sites, opnames[1], b+sep)
		  *getId(sites, b+sep+1, maxSupp+1);
	      }
	}
      }
    return hterm;
  }

};


}
#endif

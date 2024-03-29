#ifndef __ITENSOR_DMTSITESET_H
#define __ITENSOR_DMTSITESET_H
#include <itensor/mps/siteset.h>
#include <map>

namespace itensor
{


  //A wrapper for SiteSet so we can build vectorized SiteSets
  //on top of predefined ones. The base class is not vectorized.
  
  //Provides functions op, stateOp and bareOp. bareOp = op of
  //underlying SiteSet in base and derived classes.
  
  //op = bareOp for the base class, and should represent on observable or
  //transformation of a state. Has indices <Out>' <In>
  
  //stateOp represents the operator as a density matrix, with indices
  //<Out> <In>' defined such that op*stateOp = Tr(op stateOp).
  
  class DMTSiteSet
  {
    
  protected:
    SiteSet sites_;
    using StrVec = std::vector<std::string>;
    using StrITensorMap = std::map<std::string, ITensor>;
    

  public:

    DMTSiteSet(SiteSet sites): sites_(sites){}
    DMTSiteSet(){}

    const SiteSet & sites() const {return sites_;}
    const int length() const{return sites_.length();}

    virtual ITensor
    op(std::string const& opname, int i, Args const& args = Args::global()) const
    {
      return sites_.op(opname,i,args);
    }

    virtual ITensor
    stateOp(std::string const& opname, int i, Args const& args = Args::global()) const
    {
      return swapPrime(sites_.op(opname,i,args),0,1);
    }

    ITensor
    bareOp(std::string const& opname, int i, Args const& args = Args::global()) const
    {
      return sites_.op(opname,i,args);
    }


    virtual ITensor
    convertToSiteOp(const ITensor & ten, int start, int end)
    {
      return ten;
    }

    virtual ITensor
    convertToStateOp(const ITensor & ten, int start, int end)
    {
      return swapPrime(ten,0,1,"Site");
    }

    Index
    siteInd(int i)
    {
      return sites_(i);
    }

    virtual bool isVectorBasis() const { return false;}
    virtual bool isHermitianBasis() const { return false;}
    virtual StrVec vecOpNames() const {throw std::domain_error("Not vector basis");}
    virtual IndexSet vecInds() const {throw std::domain_error("Not vector basis");}
    virtual Index vecInd(const int & site_i) const {throw std::domain_error("Not vector basis");}
    virtual const std::vector<ITensor> & vecCombs () const {throw std::domain_error("Not vector basis");}
    virtual ITensor vecComb(const int & site_i) const {throw std::domain_error("Not vector basis");}
    virtual ITensor basisChange(const int & site_i) const {throw std::domain_error("Not vector basis");}

  };


  
 
  //Vectorized SiteSet. Supply vecOpNames, a list of names of the bareOps
  //of SiteSet that form a basis for the operator space. These operators
  //are cached as vectors; other operators are vectorized on the fly.
  
  //Vectorizes the index of SiteSet with {vecComb, vecInd} = combiner(ind, prime(ind)).
  //This means that unprimed indices run faster than primed.
  //By convention, *rows* in SiteSet are unprimed, (but obnoxiously
  //ordered first, counter to normal matrix convention), so
  //this means bareOps are *row* vectorized, while stateOps are
  //*column* vectorized. Column vectoization is the usual mathematical
  //definition of vec() : 

  //Definition of basisChange:
  //stateOp = swapPrime(bareOp,0,1)*vecComb*basisChange
  //op = bareOp*vecComb*conj(basisChange)

  //These definitions keep op*stateOp = tr(op stateOp) as required.
  //Can be passed in as Matrix (EXPERIMENTAL).

  //HermitianBasis: basisChange is chosen such that the stateOps are
  //of the form norm_1*(1,0,0..), norm_2*(0,1,0...), norm_3*(0,0,1...)
  //where norm_i is the square of the sum of the elements of the bareOp.

  //Only supports uniform SiteSets at the moment!
  
  class VecSiteSet : public DMTSiteSet
  {
    StrVec  vecOpNames_{};
    StrITensorMap vecOps_{};
    StrITensorMap stateOps_{};
    IndexSet vecInds_{};
    bool hermitianBasis_ = false;
    std::vector<ITensor> vecCombs_{};
    ITensor basisChange_ = ITensor(1);
    bool changeBasis_ = false;
    Index t_,dagt_,pt_;

    void constructVecInds()
    {
      const int N = length();
      vecInds_.resize(N);
      vecCombs_.resize(N);
      for(auto i : range1(N))
	std::tie(vecCombs_[i-1],vecInds_[i-1]) = combiner({dag(sites_(i)),prime(sites_(i))},
							{"Tags", "Site, n="+ std::to_string(i)});
      t_ = vecInds_[0];
      dagt_ = dag(t_);
      pt_ = prime(t_);
    }

     void constructVecOps(bool changeBasis)
    {
      if((int) vecOpNames_.size() != dim(t_))
	Error("Incorrect number  operators supplied to form complete, independent vector basis");
      for (auto & name: vecOpNames_){
	vecOps_[name] = sites_.op(name, 1)*vecCombs_[0];
	stateOps_[name] = swapPrime(sites_.op(name, 1),0,1)*vecCombs_[0];
	if (changeBasis)
	  {
	  stateOps_[name] = noPrime(stateOps_[name]*basisChange_);
	  vecOps_[name] = noPrime(vecOps_[name]*conj(basisChange_));
	  }
      }
    }

    void constructHermitianVecOps()
    {
      hermitianBasis_ = true;
      changeBasis_ = true;
      if((int) vecOpNames_.size() != dim(t_))
	Error("Incorrect number  operators supplied to form complete, independent vector basis");
      basisChange_ = ITensor({dagt_,pt_});
      for (int i = 1; i <= dim(t_); i++){
	auto name = vecOpNames_[i-1];
	auto vop = sites_.op(name, 1)*vecCombs_[0];
	stateOps_[name] = ITensor(t_);
	Real nrm = norm(vop);
	stateOps_[name].set(i,  nrm);

	//Notice here that op and stateOp are the same because
	//the implicit change in basis being transposed compensates
	//for the change in sign.
	vecOps_[name] = stateOps_[name];
	
	for (int j = 1; j <= dim(t_); j++)
	  basisChange_.set(j, i, eltC(vop, j)/nrm);
      }
    }
       
  public:

    StrVec vecOpNames() const override {return vecOpNames_;}
    IndexSet vecInds() const override {return vecInds_;}
    Index vecInd(const int & site_i) const override {return vecInds_[site_i-1];}
    const std::vector<ITensor> & vecCombs () const override {return vecCombs_;}
    ITensor vecComb(const int & site_i) const override {return vecCombs_[site_i-1];}
    bool isVectorBasis() const override {return true;}
     bool isHermitianBasis() const override { return hermitianBasis_;}
    
    
    ITensor basisChange(const int & site_i) const override {
      return site_i == 1 or not changeBasis_ ?
	basisChange_  :
	basisChange_*delta(dagt_,dag(vecInds_[site_i-1]))
	                 *delta(pt_, prime(vecInds_[site_i-1]));
    }

    VecSiteSet(const SiteSet & sites,
	       StrVec vecOpNames,
	       CMatrix matBasisChange, const Args& args): DMTSiteSet(sites),
					vecOpNames_(vecOpNames)			
    {
      constructVecInds();
      basisChange_ = matrixITensor(matBasisChange, {dagt_,pt_});
      changeBasis_ = true;
      constructVecOps(true);
      
    }

    VecSiteSet(const SiteSet & sites,
	       StrVec vecOpNames,
	       const Args& args) : DMTSiteSet(sites),
					     vecOpNames_(vecOpNames)
    {
      constructVecInds();
      if(args.getBool("HermitianBasis"))
	constructHermitianVecOps();
      else
	constructVecOps(false);
    }


    ITensor
    convertToSiteOp(const ITensor & ten, int start, int end) override
    {
      auto ret = ten;
      for (int i = std::max(start,1); i <= std::min(end, this->length()); ++i)
	{
	  ret *= vecComb(i);
	  ret *= conj(basisChange(i));
	  }
      if(hermitianBasis_){
	Real imagNorm = norm(imagPart(ret));
	if (imagNorm > 1e-12)
	    Error("Hermitian Basis but operator with AH part!");
	ret = realPart(ret);
      }
      return noPrime(ret, "Site");
    }

    ITensor
    convertToStateOp(const ITensor & ten, int start, int end) override
    {
      auto ret = swapPrime(ten,0,1,"Site");
      for (int i = start; i <= end; ++i)
	{
	  ret *= vecComb(i);
	  ret *= basisChange(i);
	  }
      if(hermitianBasis_)
	ret = realPart(ret);
      return noPrime(ret, "Site");
    }

    
    ITensor
    op(std::string const& opname, int i, Args const& args = Args::global()) const override
    {
      auto it = vecOps_.find(opname);
      return it != vecOps_.end() ?
	(i == 1 ? it->second : it->second*delta(t_,vecInds_[i-1]))
	: noPrime((sites_.op(opname,i,args)*vecComb(i))*conj(basisChange(i)));
    }


    ITensor
    stateOp(std::string const& opname, int i, Args const& args = Args::global()) const override
    {
      //Either find the Op in list of basis names or construct manually from bare Op.
      auto it = stateOps_.find(opname);
      return it != stateOps_.end() ?
	(i == 1 ? it->second : it->second*delta(t_,vecInds_[i-1]))
	: noPrime((swapPrime(sites_.op(opname,i,args),0,1)*vecComb(i))*basisChange(i));
    }

  };

  //Global namespace function overloads for DMTSiteSet
    ITensor inline
op(DMTSiteSet const& sites,
   std::string const& opname,
   int i,
   Args const& args = Args::global())
    {
    return sites.op(opname,i,args);
    }


      ITensor inline
bareOp(DMTSiteSet const& sites,
   std::string const& opname,
   int i,
   Args const& args = Args::global())
    {
    return sites.bareOp(opname,i,args);
    }

  int inline length(DMTSiteSet const & sites) {return sites.length();}

  
}

#endif

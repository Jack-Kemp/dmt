#ifndef DMT_H
#define DMT_H

#include<itensor/all.h>
#include <algorithm>
#include <cmath>
#include "DMTSiteSet.h"
#include "ITensorUtilFunctions.h"



//Class for performing DMT calculations.

//Calculations can either be done on MPO, or 'vectorized' into an
//MPS: they are equivalent. The Hermitian basis, which is real,
//can only be used vectorized.

//The main substance of the class is the custom truncation routine
//DMT::svdBond, which truncates a bond while preserving given operators.

//Conventions: pres prefix refers to preserved operators,
//those that are preserved by DMT truncation.
//'c' prefix refers to cached (precalculated) quantities.
// _ after a variable is private.

//StateOps have indices <Out> <In>' as MPO and <Out> as MPS.
//SiteOps have indices <Out>' <In>  and <Out> as MPS.
//This means siteOp*stateOp == tr(siteOp stateOp) in either language.
//(See DMTSiteSet.h for more information).

namespace itensor{

  class DMT
  {

    //The density matrix as a stateOp.
    MPO rho_;

    //The physical basis. Wraps a base itensor SiteSet, handling
    //vectorization if necessary.
    std::unique_ptr<DMTSiteSet> dmtSites_;
    bool vectorized_ = false;
    bool vectorBasis_;
    bool hermitianBasis_;

    
    //The number of sites preserved from a bond in a single direction.
    int presRadius_;

    //If presAll = False, preserve only operators in presOps_ [EXPERIMENTAL]
    //[Only works if unit cell is one site (open b.c.s are fine [CHECK])]
    bool presAll_ = true;
    std::map<int, std::vector<ITensor> > presOps_{};

    //Cache variables for trace and preservation rotation matrices Q.
    bool cacheTrace_ = true;
    Real ctrace_;   
    
    std::vector<ITensor> ctraceLeftOf_{};
    std::vector<ITensor> ctraceRightOf_{};

    std::vector<ITensor> cQpres{};
    std::vector<Index> cQpresInds{};
    std::vector<int> cnPresOps{};
    

    //Internal use for initialisation; have we added link indices yet?
    bool hasLinks_ = false;
    

    //Trace out the sites left of (not including) site presL
    ITensor
    traceLeftOf_(int presL) const
    {
      if (presL > 1) {
	if (not vectorized_)
	  {
	    return traceSubsection(rho_, 1, presL);
	  }
	else
	  {
	    auto ret = siteOp("Id",1)*rho_(1); 
	    for(int i=2; i<presL; i++)
	      {
		ret *= siteOp("Id",i)*rho_(i); 
	      }
	    return ret;
	  }
      }
      return ITensor(1);
    }

    //Trace out the sites left of (not including) site presR
    ITensor
    traceRightOf_(int presR) const
    {
      if (presR < length(rho_)){
	int lr = length(rho_);
	if (not vectorized_)
	  {
	    return traceSubsection(rho_, presR+1, lr+1);
	  }
	else
	  {
	    auto ret = siteOp("Id",lr)*rho_(lr); 
	    for(int i= lr-1; i>presR; i--)
	      {
		ret *= siteOp("Id",i)*rho_(i); 
	      }
	    return ret;
	  }
      }
      return ITensor(1);
    } 

    Real
    trace_() const{
      if (not vectorized_)
	  {
	    return traceC(rho_).real();
	  }
      else
	{
	  auto ret = siteOp("Id", 1)*rho_(1); 
	  for(int i=2; i<=length(rho_); i++)
	    {
	      ret *= siteOp("Id",i)*rho_(i); 
	    }
	  return eltC(ret).real();
	}
    }

    //Calculate and cache the preservation rotation operators Q.
    //Called on initialisation.
    void
    updateQpresCache()
    {
      if (presAll_)
	{
	  for (int radius=1; radius <= presRadius_; ++radius)
	    {
	      auto id = getId(1, radius+1);
	      auto sinds = inds(id);
	      auto [C,csiteInds] = combiner(dag(sinds));
	      //Set up a dummy index so matrix QR can be used on vector
	      Index dummyInd = hasQNs(csiteInds) ? dag(Index(qn(csiteInds,1),1)) : Index(1);
	      ITensor dummyT = ITensor(dummyInd);
	      dummyT.set(1,1.0);
	      id *= dummyT;
	      
	      Args qrArgs = Args{"Complete", true};

	      //QR the identity on pres sites to get basis with identity first el.
	      auto [Qid, Rid] = qr(C*id, csiteInds, qrArgs);
	      cQpres.push_back(Qid);
	      cQpresInds.push_back(csiteInds);
	    }
	}
      else
	{
	  for (int radius=1; radius <= presRadius_; ++radius)
	    {
	      std::vector<ITensor> radiusOps;
	      auto id = getId(1, radius+1);
	      auto sinds = inds(id);
	      auto [C,csiteInds] = combiner(dag(sinds));
	      radiusOps.push_back(C*id);
	      for (int supp = 1; supp <= radius; supp++){
		auto ops = presOps_.find(supp);
		if (ops != presOps_.end())
		  {
		    for (int start = 1; start <= radius - supp + 1; ++start)
		      {
			auto suppInds = inds(getId(start, start + supp));
			for (auto op : ops->second)
			  {
			    auto opInds = inds(op);
			    for(int i = 0; i < supp; i++)
			      if (opInds[i] != suppInds[i])
				op *= delta(opInds[i], suppInds[i]);
			    PrintData(op);
			    PrintData((getId(1, start)*op*getId(start+supp, radius+1)));
			    radiusOps.push_back(C*(getId(1, start)*op*getId(start+supp, radius+1)));
			  }
		      }
		  }
	      }
	      
	      int nPresOps = radiusOps.size();
	      //Set up a column index
	      Index colInd = hasQNs(csiteInds) ? dag(Index(qn(csiteInds,1),nPresOps)) : Index(nPresOps);
	      ITensor presTotal(csiteInds, colInd);
	      for (unsigned int i = 0; i < radiusOps.size(); i++)
		for (const auto & indvals : iterInds(radiusOps[i]))
		  presTotal.set(indvals[0].val, i+1, eltC(radiusOps[i],indvals));
			      
	      Args qrArgs = Args{"Complete", true};

	      //QR pres ops on pres sites to get basis with identity first el.
	      // and pres ops on following els.
	      auto [Qp, Rp] = qr(presTotal, csiteInds, qrArgs);
	      cQpres.push_back(Qp);
	      cQpresInds.push_back(csiteInds);
	      cnPresOps.push_back(nPresOps);
	    }
	}
    }
    
  public:


    //EXPERIMENTAL -- set preserved operators rather than calculating all.

    //Add state operator presOp with support on lattice support to preserved operators.
    void
    addPresOperator(ITensor presOp, int support, bool convertToStateOp = false)
    {
      presAll_ = false;
      if(support > presRadius_)
	{
	  Error("Support of preserved operator greater than preserved radius!");
	}
      if (convertToStateOp)
	presOps_[support].push_back(dmtSites_->convertToStateOp(presOp, 1, support+1));
      else
	presOps_[support].push_back(presOp);
    }
    
    //Number of preserved operators on radius number of sites.
    int
    presDim(const int & radius) const
    {
      return cnPresOps[radius-1];
    }

    //BASIS FUNCTIONS -- for conversion into DMT basis---------------------------------

    //Get the identity matrix as siteOp from [siteStart, siteEnd).
    ITensor
    getId(int siteStart, int siteEnd)
    {
      ITensor id = ITensor(1);
      for (int i = siteStart; i<siteEnd; i++)
	id *= siteOp("Id", i); 
      return id;
    }

    //Get the operator op_name on site_i as a siteOp (<Out>'<In>)
    ITensor
    siteOp(const char* op_name, int site_i) const
    {
      return op(*dmtSites_, op_name, site_i);
    }

    //Get the operator op_name on site_i as a stateOp (<Out><In>')
    ITensor
    stateOp(const char* op_name, int site_i) const
    {
      return dmtSites_->stateOp(op_name, site_i);
    }

    //Get the operator op supported on [site_start, site_end) as a siteOp
    ITensor
    convertToSiteOp(ITensor op, int site_start, int site_end) const
    {
      return dmtSites_->convertToSiteOp(op, site_start, site_end);
    }

    //Get the operator op supported on [site_start, site_end) as a stateOp
    ITensor
    convertToStateOp(ITensor op, int site_start, int site_end) const
    {
      return dmtSites_->convertToStateOp(op, site_start, site_end);
    }

    //Convert MPO mpo to DMT basis. An optimisation if you need to
    //do many calculations with mpo and rho.
    MPO convertToDMTBasis(const MPO & mpo) const
    {
      if(not vectorBasis_)
	return mpo;
      MPO ret(mpo);
      for(auto i : range1(length(ret)))
	{
	  ret.ref(i) = dmtSites_->convertToSiteOp(ret(i), i, i+1);
	}
      return ret;
    }


    
    //TRACE FUNCTIONS: For updating the trace cache and returning the trace---------------


    //Fully update the cache for tracing out from the left.
    void
    updateTraceCacheLeft()
    {
      ctraceLeftOf_[0] = ITensor(1);
      for (int i = 1; i < len(); i++)
	{
	  ctraceLeftOf_[i] = traceOf(i)*ctraceLeftOf_[i-1];
	}
      ctrace_ = eltC(ctraceLeftOf_[len()-1]*traceOf(len())).real(); 
    }

    //Fully update the cache for tracing out from the right.
    void
    updateTraceCacheRight()
    {
      ctraceRightOf_[len()-1] = ITensor(1);
      for (int i = len()-2; i >= 0; i--)
	{
	  ctraceRightOf_[i] = traceOf(i+2)*ctraceRightOf_[i+1];  
	}
    }

    void
    updateTraceCache()
    {
      updateTraceCacheLeft();
      updateTraceCacheRight();
    }
    

    //Update the appropriate trace cache for a one-site update at site
    //given a direction.
    void
    updateTraceCacheOneSite(int site, Direction dir)
    {
      if(dir == Fromleft)
	ctraceLeftOf_[site] = traceOf(site)*ctraceLeftOf_[site-1];
      else
	ctraceRightOf_[site-1] = traceOf(site+1)*ctraceRightOf_[site];
    }

    ITensor
    traceLeftOf(int presL) const
    {
      return cacheTrace_ ? ctraceLeftOf_[presL-1] : traceLeftOf_(presL); 
    }

    ITensor
    traceRightOf(int presR) const
    {
      return cacheTrace_ ? ctraceRightOf_[presR-1] : traceRightOf_(presR); 
    } 

    Real
    trace() const {
      return cacheTrace_ ? ctrace_ : trace_(); 
    }

    ITensor
    traceOf(int site_i) const{
      if (vectorized_)
	return op(*dmtSites_,"Id", site_i)*rho_(site_i);
      else 
	return rho_(site_i) * delta(dag(siteInds(rho_, site_i)));
    }

    //Tr(rho A). If A is already in DMT basis, convert = false.
    Complex
    trace(const MPO & A, bool convert = true) const{
      if (convert)
	return traceC(rho_, convertToDMTBasis(A));
      else
	return traceC(rho_, A);
    }

    //GETTER-SETTER METHODS -- notice "Ref" methods return modifiable references---------

     //Return the rotation Q needed to put rho in truncatable form
    ITensor
    Qpres(const int & radius, const Index & csiteInds) const
    {
      return cQpres[radius-1]*delta(cQpresInds[radius-1], csiteInds);
    }

    const DMTSiteSet &
    dmtSites () const { return *dmtSites_;}

    const SiteSet &
    sites () const { return dmtSites_->sites();}

    bool
    vectorized() const {return vectorized_;}

    bool
    vectorBasis() const {return vectorBasis_;}

    int
    len() const { return length(rho_);}

    MPO &
    rhoRef(){ return rho_; }

    const MPO &
    rho() const{ return rho_; }

    const ITensor &
    rho(const int & i) const{ return rho_(i);}

    ITensor & rhoRef(const int & i){ return rho_.ref(i); }

    int
    presRadius() const { return presRadius_; }

    //Given a two-site siteOp hterm starting at site b, calculate the corresponding gate
    //for DMT update time tDelta given DMT basis.

    //If SwapOutputs = true, constructs a swap gate after hterm gate.
    BondGate
    calcGate (ITensor hterm, double tDelta, int leftSite, const Args& args = Args::global()) const
    {
      int b = leftSite;
      auto & s = *dmtSites_;
      bool verbose = args.getBool("Verbose", false);
      bool swapOutputs = args.getBool("SwapOutputs", false);
      
      if (vectorized_)
	{
	  //If it is vectorized, we must kronecker the two Hamiltonian indices
	  //with the identity to construct the gate to perform time evolution.

	  //If (x) is kronecker product, we use vec(ABC) = C^T (x) A vec(B) and
	  //dp/dt = - i [H,p] = -i(HpI + IpH)

	  if (verbose) printfln("Vectorizing gate.");
	  auto idterm  = s.bareOp("Id",b)*s.bareOp("Id",b+1);
	  auto siteInds = IndexSet(s.siteInd(b), s.siteInd(b+1));
	  auto vecInds = IndexSet(s.vecInd(b), s.vecInd(b+1));
	  auto hsupterm = kron(idterm, hterm, siteInds, vecInds)
	    - kron(swapPrime(hterm,0,1), idterm, siteInds, vecInds);

	  //Change basis via similarity transform.
	  auto Udag = s.basisChange(b)*s.basisChange(b+1);
	  hsupterm = mapPrime(swapPrime(Udag,0,1),0,3) *(hsupterm* mapPrime(conj(Udag),1,2));
	  hsupterm.replaceTags("2","0").replaceTags("3","1");
	  
	  auto gate = vecBondGate(sites(),
				   kron(idterm, idterm, siteInds, vecInds),
				   b,b+1,BondGate::tReal,tDelta,hsupterm);

	  //Test code -- does kronned term = commutator of non-vectorised?
	  // auto testop = s.stateOp("Id",b)*s.op("Sy",b+1);
	  // auto btestop = s.bareOp("Id",b)*s.bareOp("Sy",b+1);
	  // PrintData(hsupterm*testop +
	  // 	    (mapPrime(mapPrime(hterm,1,2)*swapPrime(btestop,0,1),2,0)
	  // 	     - mapPrime(mapPrime(hterm,0,3)*swapPrime(btestop,0,1),3,1))
	  // 	    *s.vecComb(b)*s.vecComb(b+1)*Udag);
	  
	  if (hermitianBasis_){
	    Real imagNorm = norm(imagPart(gate.gate()));
	    if (verbose) printfln("Hermitian basis, truncating imaginary part norm: %f", imagNorm);
	    if (imagNorm > 1e-12)
	      Error("Hermitian Basis but imaginary gate.");
	    gate = BondGate(sites(),b,b+1, realPart(gate.gate()));
	  }
	  if(swapOutputs){
	    gate = BondGate(sites(),b,b+1,
			    mapPrime(gate.gate()
				     *delta(prime(vecInds[0]), prime(vecInds[1],2))
				     *delta(prime(vecInds[1]), prime(vecInds[0],2)),
				     2, 1)
			    );
	  }
	  return gate;
	}
      else
	{
	  //If not vectorised, we evolve via p(t) = e^(-iHt) p e^(iHt)
	  auto gp = BondGate(sites(),b,b+1,BondGate::tReal,tDelta,hterm);
	  auto gm = BondGate(sites(),b,b+1,BondGate::tReal,-tDelta,hterm);

	  if(swapOutputs)
	    {
	    gp =  BondGate(sites(),b,b+1,
			   mapPrime(gp.gate()
			   *delta(prime(dag(s.siteInd(b))), prime(s.siteInd(b+1),2))
				    *delta(prime(dag(s.siteInd(b+1))), prime(s.siteInd(b),2)),
				    2,1)
			   );
	    gm =   BondGate(sites(),b,b+1,
			   mapPrime(gm.gate()
			   *delta(prime(dag(s.siteInd(b))), prime(s.siteInd(b+1),2))
				    *delta(prime(dag(s.siteInd(b+1))), prime(s.siteInd(b),2)),
				    2,1)
			   );
	    }
	  
	  return BondGate(sites(),b,b+1,
			     mapPrime(gp.gate(),1,2) * mapPrime(gm.gate(),0,3));
	}
    }


    //The DMT truncation algorithm. Truncates a bond from b in direction dir,
    //preserving the appropriate operators.

    //Contains three SVDs with their own possible options:

    //---A truncationless "first" SVD to get singular values.
    //
    //---The main, truncating SVD in a basis defined by the preservation
    //---rotation operators Qpres.
    //
    //---A "third" SVD to bring MPO to canonical form after basis change back
    //---to physical basis, with minimal, but non-zero, truncation.

    //Ideas for good names for these are appreciated!


    //Options:
    
    //Cutoff refers to truncation cutoff (total of truncated s.v.s).

    //Absolute cutoff refers to discarding any s.v. below that value.

    //MaxDim is max bond dimension.
  
    template <typename BigMatrixT>
    Spectrum
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir,
	    BigMatrixT const & PH,
            Args args = Args::global())
    {

      //Get config params
      //auto noise = args.getReal("Noise",0.);
      //auto cutoff = args.getReal("Cutoff",MIN_CUT);


      //Use SVD or lower accuracy denmatdecomp method.
      auto useSVD = args.getBool("UseSVD",false);
      auto useSVDThird =  args.getBool("UseSVDThird",false);

      //Max Bond dimension
      int maxDim =  args.getInt("MaxDim", MAX_DIM);
      
      
      //Repack the SVD values for "non-truncating" first and third SVDs
      //The SVD value for the truncating SVD are passed through args.
      auto absolutePresCutoff = args.getBool("AbsolutePresCutoff", true);  
      auto firstCutoff = args.getReal("FirstSVDCutoff",1e-16);
      auto thirdCutoff = args.getReal("ThirdSVDCutoff",1e-16);
      auto svdMethod = args.getString("SVDMethod", "gesdd");
      args.add("SVDMethod", svdMethod);
      auto firstSVDArgs = Args{"Truncate", true, "Cutoff=",firstCutoff,
			       "AbsoluteCutoff", absolutePresCutoff, "SVDMethod", svdMethod};
      auto thirdSVDArgs = Args{"Truncate", true,"Cutoff=",thirdCutoff, "MaxDim=", maxDim,
			       "AbsoluteCutoff", absolutePresCutoff, "RespectDegenerate", true,
			       "SVDMethod", svdMethod};

      auto original_link_tags = tags(linkIndex(rho_, b));
	  
      //Check ortho centre on bond
      if(dir == Fromleft && b-1 > rho_.leftLim())
	{
	  printfln("b=%d, l_orth_lim_=%d",b,rho_.leftLim());
	  Error("b-1 > l_orth_lim_");
	}
      if(dir == Fromright && b+2 < rho_.rightLim())
	{
	  printfln("b=%d, r_orth_lim_=%d",b,rho_.rightLim());
	  Error("b+2 < r_orth_lim_");
	}

      if(presRadius_ < 1)
	Error("presRadius must be > 0 for DMT.");


      //First SVD to get singular values -- no trunc. above machine prec.
      Spectrum res;
      ITensor D;
      res = svd(AA,rho_.ref(b),D,rho_.ref(b+1), firstSVDArgs);
      auto indDL = commonIndex(rho_(b),D);
      auto indDR = commonIndex(rho_(b+1),D);

      //Find the furthest preserved site (presRadius away from bond
      //unless at edge).
      int leftmostPres = std::max(1,   b - presRadius_ + 1);
      int rightmostPres = std::min(length(rho_), b + presRadius_);
      int presRadiusL = b - leftmostPres + 1;
      int presRadiusR = rightmostPres - b;
      
      //Multiply all preserved tensors into the left and right bases.
      auto basisL = rho_(leftmostPres);
      for (int i = 1; i < presRadiusL; i ++)
	basisL *= rho_(leftmostPres + i);

      auto basisR = rho_(rightmostPres);
      for (int i = 1; i < presRadiusR; i ++)
	basisR *= rho_(rightmostPres - i);

      //Trace out the rest of rho and complete the bases
      basisL *= this->traceLeftOf(leftmostPres);
      basisR *= this->traceRightOf(rightmostPres);

      //Get physical indices (i.e. not bond)
      auto siteIndsL = uniqueInds(basisL, {D});
      auto siteIndsR =  uniqueInds(basisR, {D});

      //Vectorise physical indices
      auto [CL,csiteIndsL] = combiner(dag(siteIndsL));
      auto [CR,csiteIndsR] = combiner(dag(siteIndsR));
      int sdimL, sdimR;

      if (presAll_)
	{
	  sdimL = dim(csiteIndsL);
	  sdimR = dim(csiteIndsR);
	}
      else
	{
	  sdimL = presDim(presRadiusL);
	  sdimR = presDim(presRadiusR);
	}
      
      //If the pres. dimension is greater than current, no trunc. needed.
      if(sdimL < dim(indDL) and sdimR < dim(indDR))
	{

	  Args qrArgs = Args{"Complete", true};
	  
	  //Get basis tranform to put preserved quantities in top left.
	  auto QidL = Qpres(presRadiusL, csiteIndsL);
	  auto QidR = Qpres(presRadiusR, csiteIndsR);

	  //QR the rho left and right bases to find DMT trunc. basis.
	  auto [QbasisL, RbasisL] = qr(QidL*(dag(CL)*basisL), indDL, qrArgs);
	  auto [QbasisR, RbasisR] = qr(QidR*(dag(CR)*basisR), indDR, qrArgs);

	  auto qrLinkL = commonIndex(QbasisL, RbasisL);
	  auto qrLinkR =  commonIndex(QbasisR, RbasisR);


	  D = QbasisL * D * QbasisR;
	  
	  auto connectedComp = (D*setElt(qrLinkL=1).dag())*(D*setElt(qrLinkR=1).dag())/eltC(D,1,1);
	  D -= connectedComp;
	  

	  auto subindL = reduceDimTop(qrLinkL, sdimL);
	  auto subindR = reduceDimTop(qrLinkR, sdimR); 
	 
	  ITensor subD{subindL, subindR};
	   
	  for (int i = sdimL+1; i <= dim(qrLinkL); i++)
	    for (int j = sdimR+1; j <= dim(qrLinkR); j++)
	      {
		auto el = eltC(D,i,j);
		if (abs(el) > 0)
		  subD.set(i-sdimL,j-sdimR, el);
	      }

	  int subMaxDim = maxDim;
	  subMaxDim -= sdimL + sdimR-1;
	  if (subMaxDim <= 0)
	    {
	    printfln("Warning: MaxDim <= preservation range in DMT.");
	    subD.fill(0);
	    }
	  else
	    {
	      args.add("MaxDim", subMaxDim);
	      args.add("RespectDegenerate",args.getBool("RespectDegenerate",true));
	      if(useSVD)
		{
		  //Need high accuracy, use svd which calls the
		  //accurate SVD method in the MatrixRef library
		  ITensor W(subindL),S,V;
		  res = svd(subD,W,S,V,args);
		  subD = W*S*V;
		}
	      else
		{
		  //If we don't need extreme accuracy
		  //or need to use noise term
		  //use density matrix approach
		  ITensor W, V;
		  res = denmatDecomp(subD,W,V,dir,PH,args);
		  subD = W*V;
		}
	  }
	  
	  for (int i = sdimL+1; i <= dim(qrLinkL); i++)
	    for (int j = sdimR+1; j <= dim(qrLinkR); j++)
	      {
	      auto el = eltC(D,i,j);
	      if (abs(el) > 0)
		D.set(i,j, eltC(subD,i-sdimL,j-sdimR));
	      }


	  D += connectedComp;

	  //Reverse the basis transformation
	  auto newAA = dag(QbasisL)*D*dag(QbasisR);

	  //Third SVD. Should only remove sing. values already removed
	  // in second SVD.
	  if(useSVDThird)
	    {
	      //Need high accuracy, use svd which calls the
	      //accurate SVD method in the MatrixRef library
	      ITensor Dv, Av(indDL), Bv(indDR);
	      res = svd(newAA,Av,Dv,Bv,thirdSVDArgs);
	      rho_.ref(b) *= Av;
	      rho_.ref(b+1) *= Bv;
	      //Normalize the ortho center if requested
	      if(args.getBool("DoNormalize",false))
		{
		  Dv *= 1./itensor::norm(Dv);
		}

	      //Push the singular values into the appropriate site tensor
	      if(dir == Fromleft) rho_.ref(b+1) *= Dv;
	      else                rho_.ref(b)   *= Dv;
	    }
	  else
	    {
	      //If we don't need extreme accuracy
	      //or need to use noise term
	      //use density matrix approach
	      ITensor Av(indDL), Bv(indDR);
	      res = denmatDecomp(newAA,Av,Bv,dir,PH,thirdSVDArgs);
	      rho_.ref(b) *= Av;
	      rho_.ref(b+1) *= Bv;
	      //Normalize the ortho center if requested
	      if(args.getBool("DoNormalize",false))
		{
		  ITensor& oc = (dir == Fromleft ? rho_.ref(b+1) : rho_.ref(b));
		  auto nrm = itensor::norm(oc);
		  if(nrm > 1E-16) oc *= 1./nrm;
		}
	    }
	}
      else
	{
	  //Push the singular values into the appropriate site tensor
	  if(dir == Fromleft) rho_.ref(b+1) *= D;
	  else                rho_.ref(b)   *= D;
	}

      // Put the old tags back onto the new index
      auto lb = commonIndex(rho_(b),rho_(b+1));
      rho_.ref(b).setTags(original_link_tags,lb);
      rho_.ref(b+1).setTags(original_link_tags,lb);


      if(dir == Fromleft)
	{
	  rho_.leftLim(b);
	  if(rho_.rightLim() < b+2) rho_.rightLim(b+2);
	}
      else //dir == Fromright
	{
	  if(rho_.leftLim() > b-1) rho_.leftLim(b-1);
	  rho_.rightLim(b+1);
	}
      updateTraceCacheOneSite(b, dir);
      return res;
    }

    //Overload for svdBond without noise matrix.
     Spectrum
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir,
            Args args = Args::global())
    {
      return this->svdBond(b,AA,dir, LocalOp(), args);
    }


    //Truncate bonds so that all bonds around i are in truncated form.
      MPO& 
    position(int i, Args args = Args::global()){
      while(rho_.leftLim() < i-1)
	{
	  if(rho_.leftLim() < 0) rho_.leftLim(0);
	  //setBond(rho_.leftLim()+1); Again, only fails for write to disk, I think
	  auto WF = rho_(rho_.leftLim()+1) * rho_(rho_.leftLim()+2);
	  auto original_link_tags = tags(linkIndex(rho_,rho_.leftLim()+1));
	  svdBond(rho_.leftLim()+1,WF,Fromleft,{args,"LeftTags=",original_link_tags});
	}
      while(rho_.rightLim() > i+1)
      {
	if(rho_.rightLim() > len()+1) rho_.rightLim(len()+1);
	//setBond(rho_.rightLim()-2);
	auto WF = rho_(rho_.rightLim()-2) * rho_(rho_.rightLim()-1);
	auto original_link_tags = tags(linkIndex(rho_,rho_.rightLim()-2));
	svdBond(rho_.rightLim()-2,WF,Fromright,{args,"RightTags=",original_link_tags});
      }
      return rho_;
    }

    //Vectorize the MPO rho via the vectorization given in dmtSites_.
    //Assumes the tensors in rho are stateOps.
    void
    vec()
    {
      if(not vectorBasis_)
	Error("Sites is not vector basis, cannot vectorize.");
      if(vectorized_)
	Error("Cannot vectorize already vectorized MPO!");
      vectorized_ = true;
      for(auto i : range1(length(rho_)))
	{
	  auto & A = rho_.ref(i);
	  A *= dmtSites_->vecComb(i);
	  A *= dmtSites_->basisChange(i);
	  A.noPrime();
	if(dmtSites_->isHermitianBasis())
	  {
	  Real imagNorm = norm(imagPart(A));
	  A = realPart(A);
	  if (imagNorm > 1e-12)
	    Error("Hermitian Basis but imaginary initial state!");
	  }
	}
    }

    //Unvectorize, inverse of DMT::vec()
    void
    unvec()
    {
      if(not vectorized_)
	Error("Cannot unvectorize a not vectorized MPO!");
      vectorized_ = false;
      for(auto i : range1(length(rho_)))
	rho_.ref(i) *= dmtSites_->vecComb(i);
    }


    //Initialise density matrix as projector onto pure state |psi>
     void
    fromPureState(const MPS& psi){
       vectorized_ = false;
      rho_ = projector(psi);
      hasLinks_ = true;
      if(vectorBasis_)
	this->vec();
     }
    
    //Called after initialiastion to set caches.
        void
    finishConstruction(){
      if(not hasLinks_)
	putMPOLinks(rho_);
      // if (args.getBool("normalize"))
      // 	{
      //   rho_ *= 1/trace_();
      // 	}
      updateTraceCache();
      updateQpresCache();
    }  


    //The constructor requires underlying physical basis sites and takes
    //options on whether to vectorize. vecBasis is a list of operator names
    //needed to form a basis (see DMTSiteSet), for not vectorized call overload
    //below.
    
    DMT(SiteSet sites, std::vector<std::string> vecBasis, Args const & args = Args::global()) {
      hermitianBasis_ = args.getBool("HermitianBasis", false);
      vectorBasis_ = hermitianBasis_ or args.getBool("Vectorize", false);
      presRadius_ = args.getInt("PresRadius", 1);
      cacheTrace_ = args.getBool("CacheTrace", true);
      
      if(vectorBasis_)
	dmtSites_ = std::make_unique<VecSiteSet>(sites, vecBasis, args);
      else
	dmtSites_ = std::make_unique<DMTSiteSet>(sites);
            
      rho_ = MPO(sites);
      if(vectorBasis_){
	this->vec();
      }
      if(cacheTrace_){
	ctraceLeftOf_.resize(len());
	ctraceRightOf_.resize(len());
      }
      
    }

    DMT(SiteSet sites, Args const& args = Args::global()) : DMT(sites, {" "}, args) {}

    
  };







}

#endif

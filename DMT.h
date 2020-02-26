#ifndef DMT_H
#define DMT_H

#include<itensor/all.h>
#include <algorithm>
#include <cmath>
#include "DMTSiteSet.h"
#include "ITensorUtilFunctions.h"


namespace itensor{

  class DMT
  {
    int presRadius_;
    bool vectorized_ = false;
    bool cacheTrace_ = true;
    bool vectorBasis_;
    bool hermitianBasis_;
    Real ctrace_;
    bool hasLinks_ = false;
    std::vector<ITensor> ctraceLeftOf_{};
    std::vector<ITensor> ctraceRightOf_{};
    std::unique_ptr<DMTSiteSet> dmtSites_;

    //With two indices rho_(i) should have site indices <Out> <In>'
    //in opposite convention to normal MPOs but to match MPS <Out>,
    //as rho_ more akin to a wavefuntion than an opeator.
    MPO rho_;
    
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
    
  public:

    
    void
    updateTraceCache()
    {
      updateTraceCacheLeft();
      updateTraceCacheRight();
    }

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


    ITensor
    getId(int siteStart, int siteEnd){
      ITensor id = siteOp("Id", siteStart);
      for (int i = siteStart+1; i<siteEnd; i++)
	id *= siteOp("Id", i); 
      return id;
    }

    ITensor
    siteOp(const char* op_name, int site_i) const{
      return op(*dmtSites_, op_name, site_i);
    }

    ITensor
    stateOp(const char* op_name, int site_i) const{
      return dmtSites_->stateOp(op_name, site_i);
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

    BondGate
    calcGate (ITensor hterm, double tDelta, int leftSite, const Args& args = Args::global()) const
    {
      int b = leftSite;
      auto & s = *dmtSites_;
      bool verbose = args.getBool("Verbose", false);
      
      if (vectorized_)
	{
	  if (verbose) printfln("Vectorizing gate.");
	  auto idterm  = s.bareOp("Id",b)*s.bareOp("Id",b+1);
	  auto siteInds = IndexSet(s.siteInd(b), s.siteInd(b+1));
	  auto vecInds = IndexSet(s.vecInd(b), s.vecInd(b+1));
	  auto hsupterm = kron(idterm, hterm, siteInds, vecInds)
	    - kron(swapPrime(hterm,0,1), idterm, siteInds, vecInds);
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
	    return BondGate(sites(),b,b+1, realPart(gate.gate()));
	  }
	  return gate;
	}
      else
	{
	  auto gp = BondGate(sites(),b,b+1,BondGate::tReal,tDelta,hterm);
	  auto gm = BondGate(sites(),b,b+1,BondGate::tReal,-tDelta,hterm);
	  return BondGate(sites(),b,b+1,
			     mapPrime(gp.gate(),1,2) * mapPrime(gm.gate(),0,3));
	}
    }
    
    
    int
    presRadius() const { return presRadius_; }

    
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

      auto useSVD = args.getBool("UseSVD",false);
      auto useSVDThird =  args.getBool("UseSVDThird",false);
      int maxDim =  args.getInt("MaxDim", MAX_DIM);
      
      
      //SVD values for "non-truncating" first and third SVDs
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
      const int sdimL = dim(csiteIndsL);
      const int sdimR = dim(csiteIndsR);

      //If the pres. dimension is greater than current, no trunc. needed.
      if(sdimL < dim(indDL) and sdimR < dim(indDR))
	{
	  //Set up a dummy index so matrix QR can be used on vector
	  Index dummyInd = hasQNs(csiteIndsL) ? dag(Index(qn(csiteIndsL,1),1)) : Index(1);
	  ITensor dummyT = ITensor(dummyInd);
	  dummyT.set(1,1.0);

	  auto idL = getId(leftmostPres, b+1)*dummyT;
	  auto idR = getId(b+1, rightmostPres+1)*dummyT;

	  Args qrArgs = Args{"Complete", true};

	  //QR the identity on pres sites to get basis with identity first el.
	  auto [QidL, RidL] = qr(CL*idL, csiteIndsL, qrArgs);
	  auto [QidR, _unused2] = qr(CR*idR, csiteIndsR, qrArgs);

	  //QR the rho left and right bases to find DMT trunc. basis.
	  auto [QbasisL, RbasisL] = qr(QidL*(dag(CL)*basisL), indDL, qrArgs);
	  auto [QbasisR, RbasisR] = qr(QidR*(dag(CR)*basisR), indDR, qrArgs);

	  auto qrLinkL = commonIndex(QbasisL, RbasisL);
	  auto qrLinkR =  commonIndex(QbasisR, RbasisR);


	  // auto ret = siteOp("Id", 1)*rho_(1); 
	  // for(int i=2; i<=length(rho_); i++)
	  //   {
	  //     ret *= siteOp("Id",i)*rho_(i); 
	  //   }

	  // ret *= D;

	  // PrintData(ret);

	  D = QbasisL * D * QbasisR;

	  

	  //PrintData(AA.inds());
	  //PrintData(D.inds());
	 
	  //PrintData(eltC(D,1,1)*eltC(RbasisL,1,1)*eltC(RbasisR,1,1));

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

	  //bool full = false;
	  int subMaxDim = maxDim;
	  subMaxDim -= sdimL + sdimR-1;
	  if (subMaxDim <= 0)
	    {
	    printfln("Warning: MaxDim <= preservation range in DMT.");
	    subD.fill(0);
	    }
	  else
	    {
	      //Second SVD: trunc. non-preserved block.
	      //PrintData(maxDim);
	      args.add("MaxDim", subMaxDim);
	      //PrintData(subMaxDim);
	      args.add("RespectDegenerate",args.getBool("RespectDegenerate",true));
	      if(useSVD)
		{
		  //Need high accuracy, use svd which calls the
		  //accurate SVD method in the MatrixRef library
		  ITensor W(subindL),S,V;
		  res = svd(subD,W,S,V,args);
		  //PrintData(S.inds());
		  subD = W*S*V;
		  // if(dim(S.inds()[0]) == subMaxDim){
		  //   full = true;
		  //   subD.fill(0);
		  // }
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

	  //if(full)
	  //  PrintData(D);

	  //Reverse the basis transformation
	  auto newAA = dag(QbasisL)*D*dag(QbasisR);

	  //Third SVD. Should only remove sing. values already removed
	  // in second SVD.
	  if(useSVDThird)
	    {
	      //PrintData(thirdSVDArgs);
	      //Need high accuracy, use svd which calls the
	      //accurate SVD method in the MatrixRef library
	      ITensor Dv, Av(indDL), Bv(indDR);
	      res = svd(newAA,Av,Dv,Bv,thirdSVDArgs);
	      //PrintData(Dv.inds());
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

    void
    finishConstruction(){
      if(not hasLinks_)
	putMPOLinks(rho_);
      // if (args.getBool("normalize"))
      // 	{
      //   rho_ *= 1/trace_();
      // 	}
      updateTraceCache();
    }  

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

    void
    unvec()
    {
      if(not vectorized_)
	Error("Cannot unvectorize a not vectorized MPO!");
      vectorized_ = false;
      for(auto i : range1(length(rho_)))
	rho_.ref(i) *= dmtSites_->vecComb(i);
    }

     void
    fromPureState(const MPS& psi){
       vectorized_ = false;
      rho_ = projector(psi);
      hasLinks_ = true;
      if(vectorBasis_)
	this->vec();
     }

    
    void
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir,
            Args args = Args::global())
    {
      this->svdBond(b,AA,dir, LocalOp(), args);
    }


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

    


    DMT(SiteSet sites, std::vector<const char*> vecBasis, Args const& args = Args::global()) {
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

    
  };







}

#endif

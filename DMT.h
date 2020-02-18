#ifndef DMT_H
#define DMT_H

#include<itensor/all.h>
#include <algorithm>
#include <cmath>

namespace itensor{

  Index
  reduceDimTop(const Index & ind, const int & reduceDim)
  {
    auto indTags =  tags(ind);
    if (reduceDim > ind.dim()){
      Error("reduceDimTop: Cannot reduce index size below zero");
    }
    if (hasQNs(ind)){
      int nb = nblock(ind);
      int cumulDim = reduceDim;
      auto qns = Index::qnstorage{};
      for (int i = 1; i <= nb; i++)
	{
	  auto d = blocksize(ind, i) - cumulDim;
	  if (d > 0)
	    {
	      qns.emplace_back(qn(ind, i), d);
	      cumulDim = 0;
	    }
	  else
	    {
	      cumulDim -= blocksize(ind,i);  
	    }
	}
      return Index(move(qns), ind.dir(), indTags);
    }
    return Index(ind.dim()-reduceDim, indTags);
  }

  ITensor
  kron(const ITensor & A , const ITensor & B,
       const IndexSet & oldInds, const IndexSet & newInds)
  {
    //Assumes oldInds in pairs of form (no prime, prime) shared by A and B.
    //Doesn't kronecker indices not in oldInds so can partial kron as well
    IndexSet indsz = findInds(oldInds, "0");
    int len = length(indsz);
    std::vector<ITensor> Combs (len), pCombs (len);
    std::vector<Index> vecInds(len), vecpInds(len);
    auto ret = replaceTags(A, "1", "2")*prime(replaceTags(B,"1", "2"));
    for (int i = 0; i < len; i++) {
      std::tie(Combs[i], vecInds[i]) = combiner(indsz[i], prime(indsz[i]));
      ret *= Combs[i];
      std::tie(pCombs[i], vecpInds[i]) = combiner(prime(indsz[i],2), prime(indsz[i],3));
      ret *= pCombs[i];
    }
    ret.replaceInds(vecInds, newInds).replaceInds(vecpInds, prime(newInds));
    ret.replaceTags("3","1").replaceTags("2","1");
    return ret;
  }

  ITensor traceSubsection(const MPO& A, int start, int end){
    auto trA_n = A(start) * delta(dag(siteInds(A, start)));
    auto L = trA_n;
    for(int n=start+1; n<end; n++)
      {
        trA_n = A(n) * delta(dag(siteInds(A,n)));
        L *= trA_n;
      }
    return L;
  }

  ITensor
  getPairedId(IndexSet pairedInds)
  {
    int order = pairedInds.r();
    if (order % 2 != 0 or order < 2)
      Error("Invalid number of indices to getId.");
    ITensor id = delta(pairedInds[0],pairedInds[1]);
    for (int i =2; i < order; i+=2){
      id *=  delta(pairedInds[i], pairedInds[i+1]);
    }
    return id;
  }


    //Adapted from pull request: https://github.com/ITensor/ITensor/pull/212 
  MPO
  projector(const MPS & psi)
  {
    const int len = length(psi);
    MPO proj (len);
    auto linkBra = commonIndex(psi.A(1),psi.A(2));
    auto linkKet = commonIndex(prime(dag(psi.A(1))),prime(dag(psi.A(2))));
    auto [CL,cindl]  = combiner({linkBra, linkKet}, {"Tags=","Link, l=1"});
    proj.ref(1) =  psi(1)*prime(dag(psi(1)))*CL;
    CL = dag(CL);
    for (int i =2; i < len; i++)
      {
	auto linkBra = commonIndex(psi.A(i),psi.A(i+1));
	auto linkKet = commonIndex(prime(dag(psi.A(i))),prime(dag(psi.A(i+1))));
	auto [CR,cindr]  = combiner({linkBra, linkKet},
				    {"Tags=","Link, l=" + std::to_string(i)});
	proj.ref(i) = psi(i)*prime(dag(psi(i)))*CL*CR;
	CL = std::move(dag(CR));
      }
    proj.ref(len) = psi(len)*prime(dag(psi(len)))*CL;
    return proj;
  }

  

  BondGate
  vecMPOBondGate(SiteSet const& sites,
			ITensor const& unit,
         int i1, 
         int i2,
	 BondGate::Type type,
         Real tau,
         ITensor bondH)
    {
    if(i1 > i2)
      std::swap(i1,i2);

    if(!(type == BondGate::tReal || type == BondGate::tImag))
        {
        Error("When providing bondH, type must be tReal or tImag");
        }
    bondH *= -tau;
    
    if(type == BondGate::tReal)
        {
        bondH *= Complex_i;
        }
    auto term = bondH;
    bondH.replaceTags("1","2");
    bondH.replaceTags("0","1");
    ITensor gate;

    // exp(x) = 1 + x +  x^2/2! + x^3/3! ..
    // = 1 + x * (1 + x/2 *(1 + x/3 * (...
    // ~ ((x/3 + 1) * x/2 + 1) * x + 1
    for(int ord = 100; ord >= 1; --ord)
        {
        term /= ord;
        gate = unit + term;
        term = gate * bondH;
        term.replaceTags("2","1");
        }
    return BondGate(sites, i1,i2,gate);
    }

  class DMT
  {
    int presRange_;
    bool vectorized_ = false;
    bool cacheTrace_;
    bool vectorBasis_;
    Real ctrace_;
    bool hasLinks_ = false;
    std::vector<ITensor> ctraceLeftOf_{};
    std::vector<ITensor> ctraceRightOf_{};
    std::vector<IndexSet> siteIndsStore_{};
    SiteSet sites_;

    //With two indices rho_(i) should have site indices <Out> <In>'
    //in opposite conventio to normal MPOs but to match MPS <Out>,
    //as rho_ more akin to a wavefuntion than an opeator.
    MPO rho_;
    
    std::vector<ITensor> vecCombs;
    IndexSet vecInds;


    ITensor
    traceLeftOf_(int presL) const
    {
      if (presL > 1) {
	if (not vectorized_ and not vectorBasis_)
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
	if (not vectorized_ and not vectorBasis_)
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
      if (not vectorized_ and not vectorBasis_)
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
      if(vectorized_)
	return (op(sites_,"Id", site_i)*vecCombs[site_i-1])*rho_(site_i);
      else if (vectorBasis_)
	return op(sites_,"Idh", site_i)*rho_(site_i);
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
      if(vectorized_)
	return op(sites_, op_name, site_i)*vecCombs[site_i-1];
      else if(vectorBasis_){
	std::string op_nameh = op_name;
	return op(sites_, op_nameh + 'h', site_i);
      }
      else
	return op(sites_, op_name, site_i);
    }

    ITensor
    stateOp(const char* op_name, int site_i) const{
      if(vectorBasis_){
	return siteOp(op_name, site_i);
      }
      else if (vectorized_)
	return swapPrime(op(sites_, op_name, site_i),0,1)*vecCombs[site_i-1];
      else
	return swapPrime(op(sites_, op_name, site_i),0,1);
    }

    ITensor
     twoSiteOpH(const char* op_name_i, int site_i,
	       const char* op_name_j, int site_j) const {
      std::string opStr_i = op_name_i;
      std::string opStr_j = op_name_j;
      Args st =  {"Super", true};
      if (vectorBasis_){
	auto braOp = op(sites_, "h"+opStr_i + "Id", site_i, st)
	  *op(sites_, "h"+opStr_j + "Id", site_j, st);
	PrintData(op_name_i);
	PrintData(op_name_j);
	PrintData(braOp);
	auto ketOp = op(sites_, "hId" + opStr_i, site_i, st)
	  *op(sites_, "hId" + opStr_j, site_j, st);
	PrintData(ketOp);
	PrintData(braOp-ketOp);
	return braOp - ketOp;	
      }
      else
	{
	return op(sites_, op_name_i, site_i)*op(sites_, op_name_j, site_j);
      }
    }

    const SiteSet &
    sites () const { return sites_;}

    void
    sites(SiteSet s) { sites_ = s;}

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

    const std::vector<ITensor> &
    vecCombiners() const{ return vecCombs;}

    const ITensor &
    vecC(int i) const{ return vecCombs[i-1];}

    const Index &
    vecInd(int i) const { return vecInds[i-1];}

    BondGate
    calcGate (ITensor hterm, double tSweep, int leftSite) const
    {
      int b = leftSite;
      if (vectorized_)
	{
	  auto idterm  = twoSiteOpH("Id",b, "Id",b+1);
	  IndexSet siteInds = idterm.inds();
	  IndexSet vecInds = IndexSet(vecInd(b), vecInd(b+1));
	  auto hsupterm = kron(idterm, hterm, siteInds, vecInds)
	    - kron(hterm, idterm, siteInds, vecInds);
	  return vecMPOBondGate(sites_,
				   kron(idterm, idterm, siteInds, vecInds),
				   b,b+1,BondGate::tReal,tSweep/2.,hsupterm);
	}
      else if(vectorBasis_)
	{
	  return BondGate(sites_,b,b+1,BondGate::tReal,tSweep/2.,hterm);
	}
      else
	{
	  auto gp = BondGate(sites_,b,b+1,BondGate::tReal,tSweep/2.,hterm);
	  auto gm = BondGate(sites_,b,b+1,BondGate::tReal,-tSweep/2.,hterm);
	  return BondGate(sites_,b,b+1,
			     mapPrime(gp.gate(),1,2) * mapPrime(gm.gate(),0,3));
	}
    }
    
   
    std::vector<IndexSet>
    siteIndsStore() const { return siteIndsStore_; }
    
    Real
    presRange() const { return presRange_; }
    void
    presRange(Real pr) { presRange_ = pr; }

    
    template <typename BigMatrixT>
    Spectrum
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir,
	    BigMatrixT const & PH,
            Args args = Args::global())
    {
	  
      auto noise = args.getReal("Noise",0.);
      auto cutoff = args.getReal("Cutoff",MIN_CUT);
      auto presCutoff = args.getReal("PresCutoff",1e-15);
      auto firstSVDCutoff = args.getReal("FirstSVDCutoff",1e-15);
      auto usesvd = args.getBool("UseSVD",false);
      int maxDim =  args.getInt("MaxDim", MAX_DIM);

      auto presArgs = Args{"Cutoff=",presCutoff};
      auto firstSVDArgs = Args{"Cutoff=",firstSVDCutoff};
	  
      //From MPS svdBond
      //rho_.setBond(b); ??? I believe this only matter for write to disk
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

      if(presRange_ < 1)
	Error("presRange must be > 0 for DMT.");

      Spectrum res;
      ITensor D;

      // Store the original tags for link b so that it can
      // be put back onto the newly introduced link index
	  
      auto original_link_tags = tags(linkIndex(rho_, b));
      //PrintData(AA);
      //PrintData(rho_(b));
      //PrintData(rho_(b+1));
      res = svd(AA,rho_.ref(b),D,rho_.ref(b+1), firstSVDArgs);
      //PrintData(D.inds());
      auto indDL = commonIndex(rho_(b),D);
      auto indDR = commonIndex(rho_(b+1),D);

      //Closest left preserved site on b
      int presL = std::max(1,   b - presRange_ + 1);
      auto basisL = rho_(presL);
	  
      //Closest right preserved site on b + 1
      int presR = std::min(length(rho_), b + presRange_);
      auto basisR = rho_(presR);

      int presLenL = b - presL + 1;
      int presLenR = presR - b;
	  
      for (int i = 1; i < presLenL; i ++)
	basisL *= rho_(presL + i);
      for (int i = 1; i < presLenR; i ++)
	basisR *= rho_(presR - i);

      //Get product of tensors for the identity on all non-preserved sites
      basisL *= this->traceLeftOf(presL);
      basisR *= this->traceRightOf(presR);

      #ifdef DEBUG
      //CHECK(abs(itensor::norm(D)), 1); 
      #endif

      //Get physical indices (i.e. not bond)
      auto siteIndsL = uniqueInds(basisL, {D});
      auto siteIndsR =  uniqueInds(basisR, {D});

      //Vectorise physical indices
      auto [CL,csiteIndsL] = combiner(dag(siteIndsL));
      auto [CR,csiteIndsR] = combiner(dag(siteIndsR));
      const int sdimL = dim(csiteIndsL);
      const int sdimR = dim(csiteIndsR);

      if(sdimL < dim(indDL) and sdimR < dim(indDR))
	{

	  //Set up a dummy index so matrix QR can be used on vector
	  Index dummyInd = hasQNs(csiteIndsL) ? dag(Index(qn(csiteIndsL,1),1)) : Index(1);
	  ITensor dummyT = ITensor(dummyInd);
	  dummyT.set(1,1.0);

	  auto idL = getId(presL, b+1)*dummyT;
	  //PrintData(idL);
	  //PrintData(CL);
	  auto idR = getId(b+1, presR+1)*dummyT;

	  Args qrArgs = Args{"Complete", true};

	  auto [QidL, RidL] = qr(CL*idL, csiteIndsL, qrArgs);
	  auto [QidR, _unused2] = qr(CR*idR, csiteIndsR, qrArgs);

	  auto [QbasisL, RbasisL] = qr(QidL*(dag(CL)*basisL), indDL, qrArgs);
	  auto [QbasisR, RbasisR] = qr(QidR*(dag(CR)*basisR), indDR, qrArgs);

	  auto qrLinkL = commonIndex(QbasisL, RbasisL);
	  auto qrLinkR =  commonIndex(QbasisR, RbasisR);
	  //PrintData(D);
	  D = QbasisL * D * QbasisR;
	  auto connectedComp = (D*setElt(qrLinkL=1).dag())*(D*setElt(qrLinkR=1).dag())/eltC(D,1,1);
	  D -= connectedComp;
	  //PrintData(D.inds());

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
	  //PrintData(subD.inds());

	  int subMaxDim = maxDim;
	  subMaxDim -= sdimL + sdimR - 1;
	  //rintData(subD);
	  if (subMaxDim <= 0)
	    {
	    printfln("Warning: MaxDim <= preservation range in DMT.");
	    subD.fill(0);
	    }
	  else
	    {
	      args.add("MaxDim", subMaxDim);
	 
	      // Truncate blocks of degenerate singular values
	      args.add("RespectDegenerate",args.getBool("RespectDegenerate",true));

	      if(usesvd || (noise == 0 && cutoff < 1E-12))
		{
		  //Need high accuracy, use svd which calls the
		  //accurate SVD method in the MatrixRef library
		  ITensor W(subindL),S,V;
		  res = svd(subD,W,S,V,args);
		  //PrintData(S.inds());
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
	  //PrintData(subD);
	  for (int i = sdimL+1; i <= dim(qrLinkL); i++)
	    for (int j = sdimR+1; j <= dim(qrLinkR); j++)
	      {
	      auto el = eltC(D,i,j);
	      if (abs(el) > 0)
		D.set(i,j, eltC(subD,i-sdimL,j-sdimR));
	      }
	  //PrintData(D.inds());
	  D += connectedComp;
	  args.add("MaxDim", maxDim);
	  presArgs.add("MaxDim", maxDim);
	  

	  //A_[b] *= QbasisL.dag()*D*QbasisR.dag();
	  //auto newAA = A_[b]*A_[b+1];
	  auto newAA = dag(QbasisL)*D*dag(QbasisR);
	  //PrintData(newAA);

	  if(usesvd || (noise == 0 && presCutoff < 1E-12))
	    {
	      //Need high accuracy, use svd which calls the
	      //accurate SVD method in the MatrixRef library
	      ITensor Dv, Av(indDL), Bv(indDR);
	      res = svd(newAA,Av,Dv,Bv,presArgs);
	      rho_.ref(b) *= Av;
	      rho_.ref(b+1) *= Bv;
	      //Normalize the ortho center if requested
	      //PrintData(Dv.inds());
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
	      res = denmatDecomp(newAA,Av,Bv,dir,PH,presArgs);
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
    finishConstruction(Args args){
      if(not hasLinks_)
	putMPOLinks(rho_);
      if (args.getBool("normalize"))
	{
        rho_ *= 1/trace_();
	}
      updateTraceCache();
    }

    

    auto
    vec()
    {
      if(vectorized_ or vectorBasis_)
	Error("Cannot vectorize already vectorized MPO!");
      siteIndsStore_ = std::vector<IndexSet>(length(rho_));
      vectorized_ = true;
      //Vectorise physical indices
      vecCombs.resize(length(rho_));
      vecInds.resize(length(rho_));
      for(auto i : range1(length(rho_))){
	auto inds = siteInds(rho_, i);
	siteIndsStore_[i-1] = inds;
	std::tie(vecCombs[i-1],vecInds[i-1]) = combiner(inds,
							{"Tags", "Site, n="+ std::to_string(i)});
	rho_.ref(i) *= vecCombs[i-1];
      }
      return std::tie(vecCombs,vecInds);
    }

    void
    unvec()
    {
      if(not vectorized_)
	Error("Cannot unvectorize a not vectorized MPO!");
      vectorized_ = false;
      for(auto i : range1(length(rho_)))
	{	  
	  rho_.ref(i) *= vecCombs[i-1];
	}
      vecCombs.clear();
      vecInds.clear();
      siteIndsStore_.clear();
    }


    /*
    void
    unvec()
    {
      if(not vectorized_)
	Error("Cannot unvectorize a not vectorized MPO!");
      vectorized_ = false;
      for(auto i : range1(length(rho_)))
	{
	  Index vecind = siteIndex(rho_,i);
	  IndexSet sinds = siteIndsStore_[i-1];
	  IndexSet linds = linkInds(rho_,i); 
	  ITensor B{IndexSet(linds, sinds)};
	  int nrows = sinds[1].dim();
	  bool twolinks = order(linds) > 1;
	  for (int l0 = 1; l0 <= linds[0].dim(); l0++){
	    int u = 1;
	    int v = 1;
	    int nels = vecind.dim();
	    for (int j =1; j <= nels; j++)
	      {
		if (twolinks)
		  {
		    for (int l1 = 1; l1 <= linds[1].dim(); l1++){
		      auto el =  eltC(rho_(i),vecind=j, linds[0] = l0, linds[1] = l1);
		      B.set(l0, l1, u,v, el);
		    }
		  }
		else
		  {
		    auto el = eltC(rho_(i),vecind=j, linds[0] = l0);
		    B.set(l0, u,v, el);
		  }
		u++;
		if (u == nrows+1){
		  u = 1;
		  v++;
		}
	      }
	  }
	  rho_.ref(i) = B;
	}
      vecCombs.clear();
      vecInds.clear();
      siteIndsStore_.clear();
    }
    */


     void
    fromPureState(const MPS& psi){
      if(vectorized_)
	this->unvec();
      rho_ = projector(psi);
      hasLinks_ = true;
      if(vectorized_)
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


    DMT(SiteSet sites, Args const& args = Args::global()): sites_(sites){
      if(cacheTrace_){
	ctraceLeftOf_.resize(length(sites_));
	ctraceRightOf_.resize(length(sites_));
      }
      presRange_ = args.getInt("presRange", 1);
      cacheTrace_ = args.getBool("cacheTrace", true);
      vectorBasis_ = args.getBool("vectorBasis", false);
      rho_ = MPO(sites);
      if(args.getBool("vectorized", false))
	this->vec();
    }

    
  };







}

#endif

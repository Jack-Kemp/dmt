#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

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

  ITensor getPairedId(IndexSet pairedInds){
    int order = pairedInds.r();
    if (order % 2 != 0 or order < 2)
      Error("Invalid number of indices to getId.");
    ITensor id = delta(pairedInds[0],pairedInds[1]);
    for (int i =2; i < order; i+=2){
      id *=  delta(pairedInds[i], pairedInds[i+1]);
    }
    return id;
  }

  class DMTDensityMatrix : public MPO
  
  {
    int presRange_ = 1;
    bool vectorized = false;
    std::vector<IndexSet> siteIndsStore_{};
    
    //using MPO::MPS::N_;
    //using MPS::A_;
    // using MPS::l_orth_lim_;
    // using MPS::r_orth_lim_;

    ITensor
    get_Aid_prodL(int presL)
    {
      return traceSubsection(*this, 1, presL);
    }

    ITensor
    get_Aid_prodR(int presR)
    {
      return traceSubsection(*this, presR+1, N_+1);
    }
    
  public:

    using MPO::MPO;
    DMTDensityMatrix (MPO&& mpo) : MPO(mpo){}
    DMTDensityMatrix (const MPO& mpo): MPO(mpo){}
    
   
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
            Args & args = Args::global())
    {
	  
      auto noise = args.getReal("Noise",0.);
      auto cutoff = args.getReal("Cutoff",MIN_CUT);
      auto presCutoff = args.getReal("PresCutoff",0.0);
      auto usesvd = args.getBool("UseSVD",false);
      int maxDim =  args.getInt("MaxDim", MAX_DIM);
      auto presArgs = Args{"Cutoff=",presCutoff};
	  
      //From MPS svdBond
      setBond(b);
      if(dir == Fromleft && b-1 > leftLim())
	{
	  printfln("b=%d, l_orth_lim_=%d",b,leftLim());
	  Error("b-1 > l_orth_lim_");
	}
      if(dir == Fromright && b+2 < rightLim())
	{
	  printfln("b=%d, r_orth_lim_=%d",b,rightLim());
	  Error("b+2 < r_orth_lim_");
	}

      if(presRange_ < 1)
	Error("presRange must be > 0 for DMT.");

      Spectrum res;
      ITensor D;

      // Store the original tags for link b so that it can
      // be put back onto the newly introduced link index
	  
      auto original_link_tags = tags(linkIndex<MPO>(*this,b));
      res = svd(AA,A_[b],D,A_[b+1], {"Truncate", false});
      //PrintData(D);
      auto indDL = commonIndex(A_[b],D);
      auto indDR = commonIndex(A_[b+1],D);

      //Closest left preserved site on b
      int presL = std::max(1,   b - presRange_ + 1);
      auto basisL = A_[presL];
	  
      //Closest right preserved site on b + 1
      int presR = std::min(N_, b + presRange_);
      auto basisR = A_[presR];

      int presLenL = b - presL + 1;
      int presLenR = presR - b;
	  
      for (int i = 1; i < presLenL; i ++)
	basisL *= A_[presL + i];
      for (int i = 1; i < presLenR; i ++)
	basisR *= A_[presR - i];

      //Get product of tensors for the identity on all non-preserved sites
      if (presL > 1){
	basisL *= this->get_Aid_prodL(presL);
	//PrintData( this->get_Aid_prodL(presL));
		}
      if (presR < N_)
	basisR *= this->get_Aid_prodR(presR);

      //Get physical indices (i.e. not bond)
      auto siteIndsL = uniqueInds(basisL, {D});
      auto siteIndsR =  uniqueInds(basisR, {D});

      //Vectorise physical indices
      auto [CL,csiteIndsL] = combiner(siteIndsL);
      auto [CR,csiteIndsR] = combiner(siteIndsR);
      const int sdimL = dim(csiteIndsL);
      const int sdimR = dim(csiteIndsR);

      if(sdimL < dim(indDL) and sdimR < dim(indDR))
	{

	  //Set up a dummy index so matrix QR can be used on vector
	  Index dummyInd = hasQNs(csiteIndsL) ? Index(qn(csiteIndsL,1),1) : Index(1);
	  ITensor dummyT = ITensor(dummyInd);
	  dummyT.set(1,1.0);

	  auto idL = getPairedId(siteIndsL)*dummyT;
	  auto idR = getPairedId(siteIndsR)*dummyT;
	 

	  // int dL = std::sqrt(dim(csiteIndsL));
	  // int dR = std::sqrt(dim(csiteIndsR));

	  // //Get vectorised identities to QR into first basis vector.
	  
	  // ITensor idL{csiteIndsL, Index()};
	  // ITensor idR{csiteIndsR};
	  // for (int i = 0; i < dL; i++)
	  //   idL.set(1+(dL+1)*i, 1.0);
	  // for (int i = 0; i < dR; i++)
	  //   idR.set(1+(dR+1)*i, 1.0);

	  Args qrArgs = Args{"Complete", true};

	  auto [QidL, RidL] = qr(CL*idL, csiteIndsL, qrArgs);
	  auto [QidR, _unused2] = qr(CR*idR, csiteIndsR, qrArgs);

	  auto [QbasisL, RbasisL] = qr(QidL.dag()*(CL*basisL), indDL, qrArgs);
	  auto [QbasisR, RbasisR] = qr(QidR.dag()*(CR*basisR), indDR, qrArgs);

	  auto qrLinkL = commonIndex(QbasisL, RbasisL);
	  auto qrLinkR =  commonIndex(QbasisR, RbasisR);
	  //PrintData(D);
	  D = QbasisL.conj() * D * QbasisR;
	  auto connectedComp = (D*setElt(qrLinkL=1).dag())*(D*setElt(qrLinkR=1).dag())/eltC(D,1,1);
	  D -= connectedComp;
	  //PrintData(D);

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
	  subMaxDim -= sdimL + sdimR;
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
		D.set(i,j, eltC(D,i-sdimL,j-sdimR));
	      }
	  D += connectedComp;
	  args.add("MaxDim", maxDim);
	  presArgs.add("MaxDim", maxDim);
	  //PrintData(D);

	  //A_[b] *= QbasisL.dag()*D*QbasisR.dag();
	  //auto newAA = A_[b]*A_[b+1];
	  auto newAA = QbasisL.dag()*D*QbasisR.dag();
	  //PrintData(newAA);

	  if(usesvd || (noise == 0 && presCutoff < 1E-12))
	    {
	      //Need high accuracy, use svd which calls the
	      //accurate SVD method in the MatrixRef library
	      ITensor Dv, Av(indDL), Bv(indDR);
	      res = svd(newAA,Av,Dv,Bv,presArgs);
	      A_[b] *= Av;
	      A_[b+1] *= Bv;
	      //Normalize the ortho center if requested
	      if(args.getBool("DoNormalize",false))
		{
		  Dv *= 1./itensor::norm(Dv);
		}

	      //Push the singular values into the appropriate site tensor
	      if(dir == Fromleft) A_[b+1] *= Dv;
	      else                A_[b]   *= Dv;
	    }
	  else
	    {
	      //If we don't need extreme accuracy
	      //or need to use noise term
	      //use density matrix approach
	      ITensor Av(indDL), Bv(indDR);
	      res = denmatDecomp(newAA,Av,Bv,dir,PH,presArgs);
	      A_[b] *= Av;
	      A_[b+1] *= Bv;
	      //Normalize the ortho center if requested
	      if(args.getBool("DoNormalize",false))
		{
		  ITensor& oc = (dir == Fromleft ? A_[b+1] : A_[b]);
		  auto nrm = itensor::norm(oc);
		  if(nrm > 1E-16) oc *= 1./nrm;
		}
	    }
	  //PrintData(A_[b]);
	  //PrintDat((A_[b+1]));
	}
      else
	{
	  //Push the singular values into the appropriate site tensor
	  if(dir == Fromleft) A_[b+1] *= D;
	  else                A_[b]   *= D;
	}

      // Put the old tags back onto the new index
      auto lb = commonIndex(A_[b],A_[b+1]);
      A_[b].setTags(original_link_tags,lb);
      A_[b+1].setTags(original_link_tags,lb);


      if(dir == Fromleft)
	{
	  l_orth_lim_ = b;
	  if(r_orth_lim_ < b+2) r_orth_lim_ = b+2;
	}
      else //dir == Fromright
	{
	  if(l_orth_lim_ > b-1) l_orth_lim_ = b-1;
	  r_orth_lim_ = b+1;
	}
      return res;
    }

    auto
    vec()
    {
      if(vectorized)
	Error("Cannot vectorize already vectorized MPO!");
      siteIndsStore_ = std::vector<IndexSet>(N_);
      vectorized = true;
      std::vector<ITensor> Cs (N_);
      IndexSet cinds (N_);
      //Vectorise physical indices
      for(auto i : range1(N_)){
	auto inds = siteInds(*this, i);
	siteIndsStore_[i-1] = inds;
	std::tie(Cs[i-1],cinds[i-1]) = combiner(inds);
	A_[i] *= Cs[i-1];
      }
      return std::tie(Cs,cinds);
    }

    void
    unvec()
    {
      if(not vectorized)
	Error("Cannot unvectorize a not vectorized MPO!");
      vectorized = false;
      for(auto i : range1(N_))
	{
	  ITensor B{siteIndsStore_[i-1]};
	  int ncols = siteIndsStore_[i-1][0].dim();
	  int u = 0;
	  int v = 0;
	  int nels = A_[i].index(1).dim();
	  for (int j =0; j < nels; j++)
	    {
	      B.set(u,v, eltC(A_[i],j));
	      v++;
	      if (v == ncols){
		v = 0;
		u++;
	      }
	    }
	  A_[i] = B;
	}
    }

    void
    svdBond(int b, 
            ITensor const& AA, 
            Direction dir,
            Args & args = Args::global())
    {
      this->svdBond(b,AA,dir, LocalOp(), args);
    }



    
    
    

    
    

  };

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





}

#endif

#ifndef ITENSORUTILFUNCTIONS_H
#define ITENSORUTILFUNCTIONS_H

#include<itensor/all.h>

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
      std::tie(Combs[i], vecInds[i]) = combiner(prime(indsz[i]), indsz[i]);
      ret *= Combs[i];
      std::tie(pCombs[i], vecpInds[i]) = combiner(prime(indsz[i],3), prime(indsz[i],2));
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
  vecBondGate(SiteSet const& sites,
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

}

#endif

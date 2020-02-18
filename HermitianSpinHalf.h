//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_HERMITIANSPINHALF_H
#define __ITENSOR_HERMITIANSPINHALF_H
#include <itensor/mps/siteset.h>
#include <itensor/util/str.h>
#include <cmath>
#include<map>

namespace itensor {

  const std::map< std::string, std::vector<std::vector<Complex>> > HermitianSpinHalfSupOps_ =
    {
     { "hIdSx", {{0, 0.5, 0, 0}, {0.5, 0, 0, 0}, {0, 0, 0, 0.5*Cplx_i}, {0, 0, 0.5*Cplx_i, 0}} },
     { "hIdSy", {{0, 0, 0.5, 0}, {0, 0, 0, 0.5*Cplx_i}, {0.5, 0, 0, 0}, {0, 0.5*Cplx_i, 0, 0}} },
     { "hIdSz", {{0, 0, 0, 0.5}, {0, 0, 0.5*Cplx_i, 0}, {0, 0.5*Cplx_i, 0, 0}, {0.5, 0, 0, 0}} },
     { "hSxId", {{0, 0.5, 0, 0}, {0.5, 0, 0, 0}, {0, 0, 0, 0.5*Cplx_i}, {0, 0, 0.5*Cplx_i, 0}} },
     { "hSyId", {{0, 0, -0.5, 0}, {0, 0, 0, 0.5*Cplx_i}, {-0.5, 0, 0, 0}, {0, 0.5*Cplx_i, 0,0}} },
     { "hSzId", {{0, 0, 0, 0.5}, {0, 0,0.5*Cplx_i, 0}, {0, 0.5*Cplx_i, 0, 0}, {0.5, 0, 0, 0}} },
     { "hIdId", {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} },
     { "hIdId", {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}} }
      };

  class HermitianSpinHalfSite;

  using HermitianSpinHalf = BasicSiteSet<HermitianSpinHalfSite>;

  class HermitianSpinHalfSite
  {
    Index s;
    static constexpr Real sqrt2 = std::sqrt(2.0);
    //static const std::map< std::string, std::vector<std::vector<Real>> > supOps;
    
  public:

    HermitianSpinHalfSite(Index const& I) : s(I) { }

    HermitianSpinHalfSite(Args const& args = Args::global())
    {
      auto ts = TagSet("Site,S=1/2");
      if( args.defined("SiteNumber") )
	ts.addTags("n="+str(args.getInt("SiteNumber")));
      auto conserveqns = args.getBool("ConserveQNs",false);
      if (conserveqns)
	Error("Cannot conserveqns with Hermitian basis.");
      s = Index(4,ts);
    }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
    {
      Error("State not supported for vector basis. Use op instead. Yes, this is a gross hack");
      return IndexVal{};
    }

    ITensor
    op(std::string const& opname,
       Args const& args = Args::global()) const
    {
      //auto sP = prime(s);
	
      //auto Up = s(1);
      //auto UpP = sP(1);
      //auto Dn = s(2);
      //auto DnP = sP(2);
      ITensor Op;
      if (not args.getBool("Super", false))
	{
      std::vector<Real> opEls;
      Op = ITensor(s); 
      if(opname == "Idh")
	  opEls = {sqrt2,0,0,0};
      else
        if(opname == "Sxh")
	  opEls = {0,sqrt2/2,0,0};
      else
        if(opname == "Syh")
	  opEls = {0,0,sqrt2/2,0};
      else
        if(opname == "Szh")
	  opEls = {0,0,0,sqrt2/2};
	else
	  {
	    Error("Operator \"" + opname + "\" name not recognized");
	  }
      for (int i=0; i<4; i++)
	  Op.set(i+1, opEls[i]);
	}
      else
	{
      std::vector<std::vector<Complex>> opEls = HermitianSpinHalfSupOps_.at(opname);
      Op = ITensor(s, prime(s));
      for (int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	  Op.set(i+1, j+1, opEls[i][j]);
	}
      return Op;
    }

    //
    // Deprecated, for backwards compatibility
    //

    HermitianSpinHalfSite(int n, Args const& args = Args::global())
    {
      *this = HermitianSpinHalfSite({args,"SiteNumber=",n});
    }

  };




} //namespace itensor
#endif

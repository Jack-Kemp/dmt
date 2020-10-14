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
#ifndef __ITENSOR_RUNGBASIS_H
#define __ITENSOR_RUNGBASIS_H
#include<vector>
#include<string>
#include "itensor/mps/siteset.h"
#include "itensor/util/str.h"

namespace itensor {

class RungSite;

using Rung = BasicSiteSet<RungSite>;

class RungSite
    {
    Index s;
    public:

    RungSite(Index const& I) : s(I) { }

    RungSite(Args const& args = Args::global())
        {
        auto ts = TagSet("Site,S=1/2");
        if( args.defined("SiteNumber") )
          ts.addTags("n="+str(args.getInt("SiteNumber")));
        


	s = Index(4,ts); // UpUp, DnUp, UpDn, DnDn


	// //Conservation currently not relevant for DMT

	// auto conserveqns = args.getBool("ConserveQNs",false);
        // auto conserveSz = args.getBool("ConserveSz",conserveqns);
        // auto conserveParity = args.getBool("ConserveParity",false);
	// auto conserveLegExchange = args.getBool("ConserveLegExchange",false);

	
	// if(conserveSz and conserveParity and conserveLegExchange)
        //     {
	//       s = Index(QN({"Sz",+2},{"Parity",0,2}),{"LegExchange",0,2},1, //UpUp
	// 		QN({"Sz", 0},{"Parity",1,2}),{"LegExchange",0,2},1, //DnUp + UpDn
	// 		QN({"Sz", 0},{"Parity",1,2}),{"LegExchange",1,2},1, //UpDn - DnUp
	// 		QN({"Sz",-2},{"Parity",0,2}),{"LegExchange",0,2},1, //DnDn
	// 	      Out,ts);
        //     }
        // else if(conserveSz)
        //     {
        //     s = Index(QN({"Sz",+1}),1,
        //               QN({"Sz",-1}),1,Out,ts);
        //     }
        // else if(conserveParity)
        //     {
        //     s = Index(QN({"Parity",1,2}),1,
        //               QN({"Parity",0,2}),1,Out,ts);
        //     }
        // else
           
        }

    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state)
        {
        if(state == "Up") 
            {
            return s(1);
            }
        else 
        if(state == "Dn") 
            {
            return s(2);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IndexVal{};
        }

	ITensor
	op(std::string const& opname,
	   Args const& args = Args::global()) const
        {

	bool firstCall = true;
        auto sP = prime(s);

        auto UpUp = s(1);
        auto UpUpP = sP(1);
        auto DnUp = s(2);
        auto DnUpP = sP(2);
	auto UpDn = s(3);
        auto UpDnP = sP(3);
	auto DnDn = s(4);
        auto DnDnP = sP(4);

        auto Op = ITensor(dag(s),sP);

	if(opname == "IdId")
            {
            Op.set(UpUp,UpUpP,+1.0);
            Op.set(DnUp,DnUpP,+1.0);
	    Op.set(UpDn,UpDnP,+1.0);
            Op.set(DnDn,DnDnP,+1.0);
            }
        else
        if(opname == "SzId")
            {
            Op.set(UpUp,UpUpP,+0.5);
            Op.set(DnUp,DnUpP,-0.5);
	    Op.set(UpDn,UpDnP,+0.5);
            Op.set(DnDn,DnDnP,-0.5);
            }
        else
        if(opname == "IdSz")
	   {
            Op.set(UpUp,UpUpP,+0.5);
            Op.set(DnUp,DnUpP,+0.5);
	    Op.set(UpDn,UpDnP,-0.5);
            Op.set(DnDn,DnDnP,-0.5);
            }
        else
        if(opname == "SxId") {
            Op.set(UpUp,DnUpP,+0.5);
            Op.set(DnUp,UpUpP,+0.5);
	    Op.set(UpDn,DnDnP,+0.5);
            Op.set(DnDn,UpDnP,+0.5);
            }
        else
        if(opname == "SyId") {
            Op.set(UpUp,DnUpP,+0.5*Cplx_i);
            Op.set(DnUp,UpUpP,-0.5*Cplx_i);
	    Op.set(UpDn,DnDnP,+0.5*Cplx_i);
            Op.set(DnDn,UpDnP,-0.5*Cplx_i);
            }
        else
        if(opname == "IdSx"){
            Op.set(UpUp,UpDnP,+0.5);
            Op.set(DnUp,DnDnP,+0.5);
	    Op.set(UpDn,UpUpP,+0.5);
            Op.set(DnDn,DnUpP,+0.5);
            }
        else
        if(opname == "IdSy"){
            Op.set(UpUp,UpDnP,+0.5*Cplx_i);
            Op.set(DnUp,DnDnP,+0.5*Cplx_i);
	    Op.set(UpDn,UpUpP,-0.5*Cplx_i);
            Op.set(DnDn,DnUpP,-0.5*Cplx_i);
            }
        else
	  if(opname.size() == 4){
	    if (firstCall)
	      {
	      Op =  multSiteOps(this->op(opname.substr(0,2) +"Id"),
			      this->op("Id"+opname.substr(2,2)));
	      firstCall = false;
	      }
	    else
	      Error("Operator \"" + opname + "\" name not recognized");
	  }
        else
        if(opname == "projUpUp")
            {
            Op.set(UpUp,UpUpP,1);
            }
        else
	  if(opname == "projDnUp")
            {
            Op.set(DnUp,DnUpP,1);
            }
	  else
	    if(opname == "projUpDn")
            {
            Op.set(UpDn,UpDnP,1);
            }
	    else
	    if(opname == "projDnDn")
            {
            Op.set(DnDn,DnDnP,1);
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }

    //
    // Deprecated, for backwards compatibility
    //

    RungSite(int n, Args const& args = Args::global())
        {
        *this = RungSite({args,"SiteNumber=",n});
        }

    };

  std::vector<std::string>
      rungVectorBasis () {
	 std::vector<std::string>  shalfBasis = {"Id", "Sx", "Sy", "Sz"};
         std::vector<std::string> rungBasis;
	 for (auto u : shalfBasis)
	   for (auto l : shalfBasis)
	     rungBasis.push_back(u+l);
	 return rungBasis;
      }

} //namespace itensor
#endif

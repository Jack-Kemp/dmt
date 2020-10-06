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
#ifndef __ITENSOR_TEBD_H
#define __ITENSOR_TEBD_H

#include "itensor/mps/mpo.h"
#include "itensor/mps/bondgate.h"
#include "TEBDObserver.h"

namespace itensor {



  class CalcTEBDGate
  {
    SiteSet sites;
  public:
    CalcTEBDGate(const SiteSet & inSites): sites(inSites){}
    BondGate calcGate(ITensor hterm, double tDelta, int leftSite, const Args& args = Args::global()) const
    {
      return BondGate(sites, leftSite, leftSite+1, BondGate::tReal,tDelta,hterm);
    }

  };

//
// Evolves an MPS in real or imaginary time by an amount ttotal in steps
// of tstep using the list of bond gates provided.
//
// Arguments recognized:
//    "Verbose": if true, print useful information to stdout
//

template <class Iterable>
Real
tebdTEvol(Iterable const& gatelist, 
          Real ttotal, 
          Real tstep, 
          MPS & psi, 
          TEBDObserver& obs,
          Args args = Args::global());

//
//
// Implementations
//

template <class Iterable>
Real
tebdTEvol(Iterable const& gatelist, 
          Real ttotal, 
          Real tstep, 
          MPS & psi, 
          TEBDObserver& obs,
          Args args)
    {
    const bool verbose = args.getBool("Verbose",false);
    const bool do_normalize = args.getBool("Normalize",true);

    const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    if(verbose) 
        {
        printfln("Taking %d steps of timestep %.5f, total time %.5f",nt,tstep,ttotal);
        }

    psi.position(gatelist.front().i1());
    Real tot_norm = norm(psi);

    Real tsofar = 0;
    Real truncerr = 0;
    Spectrum spec;

    //Initial measurement
    args.add("TimeStepNum",0);
    args.add("Time",tsofar);
    args.add("TotalTime",ttotal);
    args.add("TruncError", truncerr);
    obs.interrupt(psi, args);
    
    for(auto tt : range1(nt))
        {
        auto g = gatelist.begin();
        while(g != gatelist.end())
	  {
            auto i1 = g->i1();
            auto i2 = g->i2();
            auto AA = psi(i1)*psi(i2)*g->gate();
            AA.replaceTags("Site,1","Site,0");

            ++g;
            if(g != gatelist.end())
                {
                //Look ahead to next gate position
                auto ni1 = g->i1();
                auto ni2 = g->i2();
                //SVD AA to restore MPS form
                //before applying current gate
                if(ni1 >= i2)
                    {
                    spec = psi.svdBond(i1,AA,Fromleft,args);
                    psi.position(ni1); //does no work if position already ni1
                    }
                else
                    {
                    spec = psi.svdBond(i1,AA,Fromright,args);
                    psi.position(ni2); //does no work if position already ni2
                    }
                }
            else
                {
                //No next gate to analyze, just restore MPS form
                spec = psi.svdBond(i1,AA,Fromright,args);
                }
	    truncerr += spec.truncerr();
	  }

        if(do_normalize)
            {
            tot_norm *= psi.normalize();
            }

        tsofar += tstep;

	

        args.add("TimeStepNum",tt);
        args.add("Time",tsofar);
        args.add("TotalTime",ttotal);
	args.add("TruncError", truncerr);
        obs.interrupt(psi, args);
        }
    if(verbose) 
        {
        printfln("\nTotal time evolved = %.5f\n",tsofar);
        }

    return tot_norm;

    } // tebdTEvol

} //namespace itensor


#endif

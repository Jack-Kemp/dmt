#ifndef MPOTEBD_H
#define MPOTEBD_H

#include "DMT.h"
#include <itensor/mps/bondgate.h>
#include "DMTObserver.h"

namespace itensor {


  BondGate MPOBondGate(SiteSet const& sites, 
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
    ITensor sunit = sites.op("Id",i1)*sites.op("Id",i2);
    ITensor unit = mapPrime(sunit,1,2) * mapPrime(sunit,0,3);
    
    if(type == BondGate::tReal)
      {
        bondH *= Complex_i;
      }
    auto term = bondH;
    bondH.replaceTags("2","4");
    bondH.replaceTags("0","2");
    bondH.replaceTags("3","5");
    bondH.replaceTags("1","3");
    ITensor gate;

    // exp(x) = 1 + x +  x^2/2! + x^3/3! ..
    // = 1 + x * (1 + x/2 *(1 + x/3 * (...
    // ~ ((x/3 + 1) * x/2 + 1) * x + 1
    for(int ord = 100; ord >= 1; --ord)
      {
        term /= ord;
        gate = unit + term;
        term = gate * bondH;
        term.replaceTags("5","3");
	term.replaceTags("4","2");
      }
    return BondGate(sites, i1,i2,gate);
  }



  


  //
  // Evolves an MPO in real or imaginary time by an amount ttotal in steps
  // of tstep using the list of bond gates provided.
  //
  // Arguments recognized:
  //    "Verbose": if true, print useful information to stdout
  //

  template <class Iterable>
  void
  gateTEvol(Iterable const& gatelist, 
	    Real ttotal, 
	    Real tstep, 
	    DMT & dmt, 
	    DMTObserver& obs,
	    Args args = Args::global());

  //
  //
  // Implementations
  //

  template <class Iterable>
  void
  gateTEvol(Iterable const& gatelist, 
	    Real ttotal, 
	    Real tstep, 
	    DMT & dmt, 
	    DMTObserver& obs,
	    Args args)
  {
    const bool verbose = args.getBool("Verbose",false);
    const int maxDim = args.getInt("MaxDim", 0);
    const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    const int nSweeps = args.getInt("nSweeps", 1);
    
    if(fabs(nt*tstep-ttotal) > 1E-9)
      {
        Error("Timestep not commensurate with total time");
      }
    int siteDim = dim(siteIndex(dmt.rho(), 1));
    int vecDim = siteDim;
    if (not dmt.vectorized())
      vecDim = siteDim*siteDim;
    if(maxDim < pow(vecDim, dmt.presRadius()))
      printfln("Warning: MaxDim < Dimension of preserved range in DMT.");
    
    if(verbose) 
      {
        printfln("Taking %d steps of timestep %.5f, total time %.5f",nt,tstep,ttotal);
      }
    
    dmt.updateTraceCache();
    dmt.position(gatelist.front().i1());
    dmt.updateTraceCache();
    

    Real tsofar = 0;

    //Initial measurement
    args.add("TimeStepNum",0);
    args.add("Time",tsofar);
    args.add("TotalTime",ttotal);
    obs.interrupt(dmt, args);
    
    for(auto tt : range1(nt))
      {
	for (int sweep = 0; sweep < nSweeps; ++sweep){
	  auto g = gatelist.begin();
	  while(g != gatelist.end())
            {
	      auto i1 = g->i1();
	      auto i2 = g->i2();
	      auto AA = dmt.rho(i1)*dmt.rho(i2)*g->gate();
	      if (dmt.vectorized())
		{
		  AA.replaceTags("1","0");
		}
	      else
		{
		  AA.replaceTags("Site,2","Site,0");
		  AA.replaceTags("Site,3","Site,1");
		}

	      ++g;
	      if(g != gatelist.end())
                {
		  //Look ahead to next gate position
		  auto ni1 = g->i1();
		  //auto ni2 = g->i2();
		  //SVD AA to restore MPO form
		  //before applying current gate
		  if(ni1 >= i2)
                    {
		      dmt.svdBond(i1,AA,Fromleft,args);
		      //dmt.position(ni1); //does no work if position already ni1
                    }
		  else
                    {
		      dmt.svdBond(i1,AA,Fromright,args);
		      //dmt.position(ni2); //does no work if position already ni2
                    }
                }
	      else
                {
		  //No next gate to analyze, just restore MPO form
		  dmt.svdBond(i1,AA,Fromright,args);
                }
            }
	}
	
	dmt.updateTraceCache();
        tsofar += tstep;

        args.add("TimeStepNum",tt);
        args.add("Time",tsofar);
        args.add("TotalTime",ttotal);
        obs.interrupt(dmt, args);
      }
    if(verbose) 
      {
        printfln("\nTotal time evolved = %.5f\n",tsofar);
      }

  } // gateTEvol

} //namespace itensor


#endif


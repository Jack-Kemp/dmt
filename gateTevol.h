#ifndef MPOTEBD_H
#define MPOTEBD_H

#include "DMT.h"
#include <itensor/mps/bondgate.h>
#include "DMTObserver.h"

namespace itensor {

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
    const bool verbose = args.getBool("Verbose");
    const int maxDim = args.getInt("MaxDim");
    const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    const int nSweeps = args.getInt("nSweeps");
    
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
    dmt.position(gatelist.front().i1(), args);
    dmt.updateTraceCache();
    

    Real tsofar = args.getReal("tStart");
    Real truncerr = 0;
    Spectrum spec;
    

    //Initial measurement
    args.add("TimeStepNum",0);
    args.add("Time",tsofar);
    args.add("TotalTime",ttotal);
    args.add("TruncError", truncerr);
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
		      spec = dmt.svdBond(i1,AA,Fromleft,args);
                    }
		  else
                    {
		      spec = dmt.svdBond(i1,AA,Fromright,args);
                    }
                }
	      else
                {
		  //No next gate to analyze, just restore MPO form
		  spec = dmt.svdBond(i1,AA,Fromright,args);
                }
	      truncerr += spec.truncerr();
            }	  
	}

	
	dmt.updateTraceCache();
        tsofar += tstep;

        args.add("TimeStepNum",tt);
        args.add("Time",tsofar);
        args.add("TotalTime",ttotal);
	args.add("TruncError", truncerr);
        obs.interrupt(dmt, args);
      }
    if(verbose) 
      {
        printfln("\nTotal time evolved = %.5f\n",tsofar);
      }

  } // gateTEvol

} //namespace itensor


#endif


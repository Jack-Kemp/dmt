#ifndef MPOTEBD_H
#define MPOTEBD_H

#include "density_matrix.h"
#include <itensor/mps/bondgate.h>
#include <itensor/mps/TEvolObserver.h>

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
          DMTDensityMatrix & psi, 
          Args const& args = Args::global());

template <class Iterable>
void
gateTEvol(Iterable const& gatelist, 
          Real ttotal, 
          Real tstep, 
          DMTDensityMatrix & psi, 
          Observer& obs,
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
          DMTDensityMatrix & psi, 
          Observer& obs,
          Args args)
    {
    const bool verbose = args.getBool("Verbose",false);
    const int MaxDim = args.getInt("MaxDim", 0);
    const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }
    
    int siteDim = dim(siteIndex(psi, 1));
    if(MaxDim < pow(siteDim*siteDim, psi.presRange()))
       printfln("Warning: MaxDim < Dimension of preserved range in DMT.");
    
    if(verbose) 
        {
        printfln("Taking %d steps of timestep %.5f, total time %.5f",nt,tstep,ttotal);
        }
    if(verbose) 
        {
        printfln("Taking %d steps of timestep %.5f, total time %.5f",nt,tstep,ttotal);
        }


    psi.position(gatelist.front().i1());

    Real tsofar = 0;
    for(auto tt : range1(nt))
        {
        auto g = gatelist.begin();
        while(g != gatelist.end())
            {
            auto i1 = g->i1();
            auto i2 = g->i2();
            auto AA = psi(i1)*psi(i2)*g->gate();
            AA.replaceTags("Site,2","Site,0");
	    AA.replaceTags("Site,3","Site,1");

            ++g;
            if(g != gatelist.end())
                {
                //Look ahead to next gate position
                auto ni1 = g->i1();
                auto ni2 = g->i2();
                //SVD AA to restore MPO form
                //before applying current gate
                if(ni1 >= i2)
                    {
                    psi.svdBond(i1,AA,Fromleft,args);
                    psi.position(ni1); //does no work if position already ni1
                    }
                else
                    {
                    psi.svdBond(i1,AA,Fromright,args);
                    psi.position(ni2); //does no work if position already ni2
                    }
                }
            else
                {
                //No next gate to analyze, just restore MPO form
                psi.svdBond(i1,AA,Fromright,args);
                }
            }

        tsofar += tstep;

        args.add("TimeStepNum",tt);
        args.add("Time",tsofar);
        args.add("TotalTime",ttotal);
        obs.measure(args);
        }
    if(verbose) 
        {
        printfln("\nTotal time evolved = %.5f\n",tsofar);
        }

    } // gateTEvol

template <class Iterable>
void
gateTEvol(Iterable const& gatelist, 
          Real ttotal, 
          Real tstep, 
          DMTDensityMatrix & psi, 
          Args const& args)
    {
    TEvolObserver obs(args);
    return gateTEvol(gatelist,ttotal,tstep,psi,obs,args);
    }

} //namespace itensor


#endif


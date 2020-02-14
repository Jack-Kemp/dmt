#ifndef __ITENSOR_DMTOBSERVER_H
#define __ITENSOR_DMTOBSERVER_H
#include <itensor/util/readwrite.h>
#include "DMT.h"
#include<string>
#include<map>
#include<functional>

namespace itensor{

class DMTObserver
    {
      std::function<void(DMT & dmt, Args const& args)> measure;
    public:
    
      DMTObserver(std::function<void(DMT & dmt, Args const& args)>,
		  Args const& args = Args::global());

    virtual ~DMTObserver() { }

    void virtual
    interrupt(DMT& dmt, Args const& args = Args::global());
    
    bool virtual
    checkDone(Args const& args = Args::global());

    private:

    /////////////
    //
    // Data Members

    bool done_,
         show_percent_;

    //
    /////////////

    }; // class DMTObserver

inline DMTObserver::
DMTObserver(std::function<void(DMT & dmt, Args const& args)> mfunc,
	     const Args& args) 
    : 
    measure(mfunc),
    done_(false),
    show_percent_(args.getBool("ShowPercent",true))
    {
      
    }


void inline DMTObserver::
interrupt(DMT& dmt, const Args& args)
    {
      measure(dmt, args);
      if(show_percent_)
        {
	  const Real t = args.getReal("Time");
	  const Real ttotal = args.getReal("TotalTime");
	  Real percentdone = (100.*t)/ttotal;
	  if(percentdone < 99.5 || (std::fabs(t-ttotal) < 1E-10))
            {
	      printf("\b\b\b%2.f%%",percentdone);
	      std::cout.flush();
            }
        }
      checkDone(args);
    }
  

bool inline DMTObserver::
checkDone(const Args& args)
    {
    const Real t = args.getReal("Time");
    if(fileExists("STOP_DMT"))
        {
        println("File STOP_DMT found: stopping this time evolution run at time ",t);
        std::remove("STOP_DMT");
        return true;
        }

    //Set done_ flag to true so any outer callers using this Observer will also terminate.
    if(fileExists("STOP_DMT_ALL"))
        {
        println("File STOP_DMT_ALL found: stopping this time evolution at time ",t);
        std::remove("STOP_DMT_ALL");
        done_ = true;
        return done_;
        }
    
    return done_;
    }

} //namespace itensor


#endif // __ITENSOR_DMTOBSERVER_H

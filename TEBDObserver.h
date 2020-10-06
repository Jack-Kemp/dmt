#ifndef __ITENSOR_TEBDOBSERVER_H
#define __ITENSOR_TEBDOBSERVER_H
#include<itensor/all.h>
#include<string>
#include<map>
#include<functional>

namespace itensor{

  //A functor called during TEBD time evolution for real time control
  //of the algorithm, e.g. for adaptive stopping.

class TEBDObserver
    {
      std::function<void(MPS& psi, Args const& args)> measure;
    public:
    
      TEBDObserver(std::function<void(MPS& psi, Args const& args)>,
		  Args const& args = Args::global());

    virtual ~TEBDObserver() { }

    void virtual
    interrupt(MPS& psi, Args const& args = Args::global());
    
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

    }; // class TEBDObserver

inline TEBDObserver::
TEBDObserver(std::function<void(MPS& psi, Args const& args)> mfunc,
	     const Args& args) 
    : 
    measure(mfunc),
    done_(false),
    show_percent_(args.getBool("ShowPercent",true))
    {
      
    }


void inline TEBDObserver::
interrupt(MPS& psi, const Args& args)
    {
      measure(psi, args);
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
  

bool inline TEBDObserver::
checkDone(const Args& args)
    {
    const Real t = args.getReal("Time");
    if(fileExists("STOP_TEBD"))
        {
        println("File STOP_TEBD found: stopping this time evolution run at time ",t);
        std::remove("STOP_TEBD");
        return true;
        }

    //Set done_ flag to true so any outer callers using this Observer will also terminate.
    if(fileExists("STOP_TEBD_ALL"))
        {
        println("File STOP_TEBD_ALL found: stopping this time evolution at time ",t);
        std::remove("STOP_TEBD_ALL");
        done_ = true;
        return done_;
        }
    
    return done_;
    }

} //namespace itensor


#endif // __ITENSOR_TEBDOBSERVER_H

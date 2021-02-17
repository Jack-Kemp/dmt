#ifndef __ITENSOR_DMTOBSERVER_H
#define __ITENSOR_DMTOBSERVER_H
#include <itensor/util/readwrite.h>
#include "DMT.h"
#include<string>
#include<chrono>
#include<map>
#include<functional>

namespace itensor{

  //A functor called during DMT time evolution for real time control
  //of the algorithm, e.g. for adaptive stopping.

  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = std::chrono::time_point<Clock>;
  using Duration = std::chrono::duration<double>;
  constexpr TimePoint invalidTime_k = TimePoint::max();

class DMTObserver
    {
       
    public:

      Real getRunTime() {return totalTime_.count();}

      DMTObserver(std::function<void(DMT & dmt, Args const& args)>);
    
      DMTObserver(std::function<void(DMT & dmt, Args const& args)>,
		  Args const& args);

      DMTObserver(std::function<void(DMT & dmt, Args const& args)>,
		  std::function<void(DMT & dmt, Args const& args)>,
		  Args const& args);

    virtual ~DMTObserver() { }

    void virtual
    interrupt(DMT& dmt, Args const& args = Args::global());
    
    bool virtual
    checkDone(Args const& args = Args::global());

    private:

    /////////////
    //
    // Data Members

      bool done_, show_percent_;
      TimePoint startTime_, checkTime_;
      Duration totalTime_;
      bool checkpoint_;
      Real checkpointTime_;
      std::string checkpointName_, outputDir_;
      std::function<void(DMT & dmt, Args const& args)> measure;
      std::function<void(DMT & dmt, Args const& args)> write_checkpoint;

    //
    /////////////

    }; // class DMTObserver


inline DMTObserver::
DMTObserver(std::function<void(DMT & dmt, Args const& args)> mfunc) 
    :
  done_(false),
  show_percent_(true),
  startTime_(invalidTime_k),
  checkpoint_(false),
  measure(mfunc)  
    {
    }

inline DMTObserver::
DMTObserver(std::function<void(DMT & dmt, Args const& args)> mfunc,
	     const Args& args) 
    :
  done_(false),
  show_percent_(args.getBool("ShowPercent",true)),
  startTime_(invalidTime_k),
  checkpoint_(args.getBool("Checkpoint")),
  checkpointTime_(args.getReal("CheckpointTime")),
  checkpointName_(args.getString("CheckpointName")),
  outputDir_(args.getString("OutputDir")),
  measure(mfunc)  
    {
      write_checkpoint = [&](DMT& dmt, Args const & args){
			   dmt.writeToFile(args.getString("WriteFilename"));
			 };
    }

inline DMTObserver::
DMTObserver(std::function<void(DMT & dmt, Args const& args)> mfunc,
	    std::function<void(DMT & dmt, Args const& args)> wfunc,
	     const Args& args) 
    : 
    done_(false),
    show_percent_(args.getBool("ShowPercent",true)),
    startTime_(invalidTime_k),
    checkpoint_(args.getBool("Checkpoint")),
    checkpointTime_(args.getReal("CheckpointTime")),
    checkpointName_(args.getString("CheckpointName")),
    outputDir_(args.getString("OutputDir")),
    measure(mfunc),
    write_checkpoint(wfunc)
    {
    }


void inline DMTObserver::
interrupt(DMT& dmt, const Args & args)
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
      if (startTime_ == invalidTime_k){
	startTime_ = Clock::now();
	checkTime_ = startTime_;
      }
      auto now = Clock::now();
      auto runTime = std::chrono::duration_cast<Duration>(now - checkTime_);
      totalTime_ =  std::chrono::duration_cast<Duration>(now - startTime_);
      
      if (checkpoint_ and runTime.count()/3600 > checkpointTime_){
	Args checkArgs = args;
	std::string wfilename = outputDir_ + "/" + checkpointName_ + "_t_" + std::to_string(args.getReal("Time"));
	checkArgs.add("WriteFilename", wfilename);
	checkArgs.add("WallTime", totalTime_.count());
	write_checkpoint(dmt, checkArgs);
	checkTime_ = now;
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

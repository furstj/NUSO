#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <sys/time.h>
#include <sys/resource.h>


class Timer {
private:
  std::string name;
  double wallTime;
  double userTime;
  long   numberOfRuns;
  bool   isRunning;
  double lastWallTime;
  double lastUserTime;

  double systemUserTime() {
    struct rusage  rusage;
    getrusage(RUSAGE_SELF, &rusage);
    return (double)(rusage.ru_utime.tv_sec) + 
      (double)(rusage.ru_utime.tv_usec)*1.0e-06;
  }

  double systemWallTime() {
    struct timeval timeval;
    gettimeofday(&timeval, 0);
    return (double)(timeval.tv_sec) + (double)(timeval.tv_usec)*1.0e-06;
  }

public:


  Timer() : numberOfRuns(0), wallTime(0), userTime(0), 
	    isRunning(false), name("unnamed") {};

  Timer(std::string name) : numberOfRuns(0), wallTime(0), userTime(0), 
	    isRunning(false), name(name) {};

  ~Timer() {
    std::cout << "TIMER " << name << ":" 
      " user time: " << userTime << "s" << 
      " wall time: " << wallTime << "s" <<
      " #calls: " << numberOfRuns << std::endl;
  } 

  void start() {
    isRunning = true;
    lastWallTime = systemWallTime();
    lastUserTime = systemUserTime();
  }

  void stop() {
    if (isRunning) {
      isRunning = false;
      userTime += systemUserTime() - lastUserTime;
      wallTime += systemWallTime() - lastWallTime;
      numberOfRuns++;
    }
  }

  double getWallTime() {
    if (isRunning) 
      return systemWallTime() - lastWallTime;
    else
      return wallTime;
  }

};

#endif

#pragma once

// A C++ microsecond-accurate timer, that works like a stopwatch.
// Call start() to start counting, then stop() to get the elapsed time.
// If desired, the accuracy can be changed in the stop() method.
//
// Usage:
// 1: Instantiate.
// 2: Start counting.
// 3: Get the time elapsed.
// 
// Timer timer;
// timer.start();
// double elapsed_time = timer.stop();

#include <chrono>

class Timer 
{
public:
  void start();
  double stop();
private:
  std::chrono::high_resolution_clock::time_point t1;
  std::chrono::high_resolution_clock::time_point t2;
};

void Timer::start()
{
  t1 = std::chrono::high_resolution_clock::now();
}

double Timer::stop()
{
  t2 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> change_ = t2 - t1;
  double change 
    = std::chrono::duration_cast<std::chrono::microseconds>(change_).count()
    * 1.0e-6;
  return change;
}
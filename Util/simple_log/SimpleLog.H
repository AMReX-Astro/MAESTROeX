#ifndef _SIMPLE_LOG_H_
#define _SIMPLE_LOG_H_

#include <AMReX_Print.H>
#include <MaestroUtil.H>
#include <state_indices.H>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/*
the SimpleLog class provides a simple buffer that stores strings
and writes to the screen at the same time.

at the moment, the log size is fixed to MAX_LINES -- in the future
we can update it reallocate the buffer as needed.
*/

class SimpleLog {
   public:
    SimpleLog() {
        log_data.resize(MAX_LINES);
        log_lines = 0;
    };

    void Log(const std::string& str, const amrex::Real d);
    void Log(const std::string& str, const int i);
    void Log(const std::string& str);

    void LogBreak();

   private:
    static constexpr int MAX_LINES = 128;

    int log_lines;

    std::vector<std::string> log_data;
};

#endif
#include <SimpleLog.H>

using namespace amrex;

void 
SimpleLog::Log(const std::string str, const amrex::Real d)
{
    Log(str + std::to_string(d));
}


void 
SimpleLog::Log(const std::string str, const int i)
{
    Log(str + std::to_string(i));
}


void 
SimpleLog::Log(const std::string str)
{
    if (log_lines >= MAX_LINES) {
        Abort("ERROR: log bugger exceeded in SimpleLog");
    }

    // output to the screen
    Print() << MAESTROeX::trim(str) << std::endl;

    // and store in the log
    log_data[log_lines] = MAESTROeX::trim(str);
    log_lines++;
}

void 
SimpleLog::LogBreak()
{
    if (log_lines >= MAX_LINES) {
        Abort("ERROR: log bugger exceeded in SimpleLog");
    }

    log_data[log_lines] = "--------------------------------------------------------------------------------";
    log_lines++;
}
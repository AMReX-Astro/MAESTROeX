#include <SimpleLog.H>

using namespace amrex;

string 
trim(const string& str) {
    const char* whitespace = " \t\v\r\n";
    const auto str_begin = str.find_first_not_of(whitespace);
    const auto str_end = str.find_last_not_of(whitespace);
    return str_begin == str_end ? "" : str.substr(str_begin, str_end - str_begin - 1);
}

void 
SimpleLog::Log(const string str, const amrex::Real d)
{
    Log(str + to_string(d));
}


void 
SimpleLog:Log(const string str, const int i)
{
    Log(str + to_string(d));
}


void 
SimpleLog:Log(const string str)
{
    if (log_lines >= MAX_LINES) {
        Abort("ERROR: log bugger exceeded in SimpleLog");
    }

    // output to the screen
    Print() << trim(str) << std::endl;

    // and store in the log
    log_data[log_lines] = trim(str);
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
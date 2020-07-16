#include <MaestroUtil.H>

// this trims the whitespace from either side of a string
std::string maestro::trim(const std::string& str) {
    const char* whitespace = " \t\v\r\n";
    const auto str_begin = str.find_first_not_of(whitespace);
    const auto str_end = str.find_last_not_of(whitespace);
    return str_begin == str_end
               ? ""
               : str.substr(str_begin, str_end - str_begin + 1);
}
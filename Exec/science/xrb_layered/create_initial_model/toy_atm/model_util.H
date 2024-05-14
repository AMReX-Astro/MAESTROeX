#ifndef MODEL_UTIL_H
#define MODEL_UTIL_H

std::string num_to_unitstring(const double value) {

    std::string numstr{};
    std::stringstream ss;

    if (value > 1.e5) {
        // work in km
        ss << std::setprecision(2) << value / 1.e5;
        numstr = ss.str() + "km";
    } else {
        // work in cm
        ss << std::setprecision(3) << value;
        numstr = ss.str() + "cm";
    }

    return numstr;
}



#endif

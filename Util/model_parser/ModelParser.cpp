#include <ModelParser.H>

using namespace amrex;

string 
trim(const string& str) {
    const char* whitespace = " \t\v\r\n";
    const auto str_begin = str.find_first_not_of(whitespace);
    const auto str_end = str.find_last_not_of(whitespace);
    return str_begin == str_end ? "" : str.substr(str_begin, str_end - str_begin - 1);
}

void 
ModelParser::ReadFile(const string& model_file_name)
{
    // open the model file 
    ifstream model_file(model_file_name);
    if (!model_file.is_open()) {
        Abort("Could not open model file!")
    }

    // the first line has the number of points in the model
    string line;
    getline(model_file, line);
    int ipos = line.find('=') + 1;
    npts_model = stoi(line.substr(ipos));

    // now read in the number of variables
    getline(model_file, line);
    ipos = line.find('=') + 1;
    static const int nvars_model_file = stoi(line.substr(ipos));

    RealVector vars_stored(nvars_model_file);
    Vector<string> varnames_stored(nvars_model_file);

    // now read in the names of the variables
    for (auto i = 0; i < nvars_model_file; ++i) {
        getline(model_file, line);
        ipos = line.find('#') + 1;
        varnames_stored[i] = trim(line.substr(ipos));
    }

    // alocate storage for the model data 
    model_state.resize(npts_model);
    for (auto i = 0; i < npts_model; ++i) {
        model_state[i].resize(nvars_model);
    }
    model_r.resize(npts_model);

    Print() << "\n\nreading initial model" << std::endl;
    Print() << npts_model << " points found in the initial model" << std::endl;
    Print() << nvars_model_file << " variables found in the initial model" << std::endl;

    // start reading in the data 
    for (auto i = 0; i < npts_model; ++i) {
        getline(model_file, line);
        istringstream iss(line);
        iss >> model_r[i];
        for (auto j = 0; j < nvars_model_file; ++j) {
            iss >> vars_stored[j];
        }

        for (auto j = 0; j < nvars_model; ++j) {
            model_state[i][j] = 0.0;
        }

        // make sure that each of the variables that MAESTROeX cares about
        // are found
        bool found_dens = false;
        bool found_temp = false;
        bool found_pres = false;
        Vector<bool> found_spec(NumSpec);
        for (auto comp = 0; comp < NumSpec; ++comp) {
            found_spec[comp] = false;
        }

        for (auto j = 0; j < nvars_model_file; ++j) {
            // keep track of whether the current variable from the model 
            // file is the one that MAESTROeX cares about

            bool found_model = false;

            if (varnames_stored[j] == "density") {
                model_state[i][idens_model] = vars_stored[j];
                found_model = true;
                found_dens = true;
            } else if (varnames_stored[j] == "temperature") {
                model_state[i][itemp_model] = vars_stored[j];
                found_model = true;
                found_temp = true;
            } else if (varnames_stored[j] == "pressure") {
                model_state[i][ipres_model] = vars_stored[j];
                found_model = true;
                found_pres = true;
            } else {
                for (auto comp = 0; comp < NumSpec; ++comp) {
                    if (varnames_stored[j] == spec_names_cxx[comp]) {
                        model_state[i][ispec_model+comp] = vars_stored[j];
                        found_model = true;
                        found_spec[comp] = true;
                    }
                }
            }

            // is the current variable from the model file one that we
            // care about?
            if (!found_model && i == 0) {
                Print() << "WARNING: variable not found: " << trim(varnames_stored[j]) << std::endl;
            }
        }

        // were all the variable that we care about provided?
        if (i == 0) {
            if (!found_dens) {
                Print() << "WARNING: density not provided in inputs file" << std::endl;
            }
            if (!found_temp) {
                Print() << "WARNING: temperature not provided in inputs file" << std::endl;
            }
            if (!found_pres) {
                Print() << "WARNING: pressure not provided in inputs file" << std::endl;
            }
            for (auto comp = 0; comp < NumSpec; ++comp) {
                if (!found_spec[comp]) {
                    Print() << "WARNING: " << trim(spec_names_cxx[comp]) << " not provided in inputs file" << std::endl;
                }
            }
        }
    }

    model_initialized = true;
}

Real 
ModelParser::Interpolate(const Real r, const int ivar, 
                         bool interpolate_top)
{
    // use the module's array of model coordinates (model_r), and
    // variables (model_state), to find the value of model_var at point
    // r using linear interpolation.

    Real interpolate = 0.0;

    // find the location in the coordinate array where we want to interpolate
    int i = 0;
    while (model_r[i] < r && i <= npts_model) {
        i++;
    }
    if (i > 0 && i < npts_model) {
        if (fabs(r - model_r[i-1]) < fabs(r - model_r[i])) {
            i--;
        }
    } else if (i == npts_model) {
        i--;
    }

    if (i == 0) {
        Real slope = (model_state[i+1][ivar] - model_state[i][ivar]) / (model_r[i+1] - model_r[i]);
        interpolate = slope * (r - model_r[i]) + model_state[i][ivar];

        // safety check to make sure interpolate lies within the bounding points
        Real minvar = min(model_state[i+1][ivar], model_state[i][ivar]);
        Real maxvar = max(model_state[i+1][ivar], model_state[i][ivar]);
        interpolate = max(interpolate, minvar);
        interpolate = min(interpolate, maxvar);
    } else if (i == npts_model - 1) {
        Real slope = (model_state[i][ivar] - model_state[i-1][ivar]) / (model_r[i] - model_r[i-1]);
        interpolate = slope * (r - model_r[i]) + model_state[i][ivar];

        // safety check to make sure interpolate lies within the bounding points
        if (!interp_top) {
            Real minvar = min(model_state[i][ivar], model_state[i-1][ivar]);
            Real maxvar = max(model_state[i][ivar], model_state[i-1][ivar]);
            interpolate = max(interpolate, minvar);
            interpolate = min(interpolate, maxvar);
        }
    } else {
        if (r >= model_r[i]) {
            Real slope = (model_state[i+1][ivar] - model_state[i][ivar]) / (model_r[i+1] - model_r[i]);
            interpolate = slope * (r - model_r[i]) + model_state[i][ivar];

            // safety check to make sure interpolate lies within the bounding points
            Real minvar = min(model_state[i+1][ivar], model_state[i][ivar]);
            Real maxvar = max(model_state[i+1][ivar], model_state[i][ivar]);
            interpolate = max(interpolate, minvar);
            interpolate = min(interpolate, maxvar);
        } else {
            Real slope = (model_state[i][ivar] - model_state[i-1][ivar]) / (model_r[i] - model_r[i-1]);
            interpolate = slope * (r - model_r[i]) + model_state[i][ivar];

            // safety check to make sure interpolate lies within the bounding points
            Real minvar = min(model_state[i][ivar], model_state[i-1][ivar]);
            Real maxvar = max(model_state[i][ivar], model_state[i-1][ivar]);
            interpolate = max(interpolate, minvar);
            interpolate = min(interpolate, maxvar);
        }
    }

    return interpolate;
}
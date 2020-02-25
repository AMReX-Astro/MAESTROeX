#ifndef _network_H_
#define _network_H_

#include <string>
#include <network_properties.H>

static int network_spec_index(const std::string name) {
  int idx = -1;
  for (int n = 0; n < NumSpec; n++) {
    if (name == spec_names_cxx[n]) {
      idx = n;
      break;
    }
  }
  return idx;
}

#endif

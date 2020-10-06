#ifndef _state_indices_H_
#define _state_indices_H_
#include <network_properties.H>

  constexpr int Nscal = 4 + NumSpec + NumAux;

   // scalars
  constexpr int Rho = 0;
  constexpr int RhoH = 1;
  constexpr int Temp = 2;
  constexpr int Pi = 3;
  constexpr int FirstSpec = 4;
  constexpr int FirstAux = 4 + NumSpec;

#endif

#ifndef _MODEL_PARSER_F_H_
#define _MODEL_PARSER_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

void ca_read_model_file(int* model_file, const int* namelen);

#ifdef __cplusplus
}
#endif

#endif
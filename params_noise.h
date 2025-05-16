#include "comm.h"

#define X 0
#define Y 1
#define Z 2

const float BIAS_ACC[3] = {(0.02),(0.02),(0.02)}; // m/s2
const float NOISE_ACC[3] = {(0.012),(0.012),(0.012)}; // m/s2

const float BIAS_GYR[3] = {(0.75),(0.75),(0.75)}; // rad
const float NOISE_GYR[3] = {(0.049),(0.049),(0.049)}; // rad

const float RTK_NOISE[3] = {0.02, 0.02, 0.03};// m
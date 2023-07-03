#ifndef UserTypes_h
#define UserTypes_h
#include "Arduino.h"
#define FILE_BASE_NAME "mpubmpraw"
const size_t ADC_COUNT = 1;
struct data_t {
  unsigned long time;
  int16_t ax;
  int16_t ay;
  int16_t az;
  int16_t gx;
  int16_t gy;
  int16_t gz;
  double altitude;
};
#endif  // UserTypes_h

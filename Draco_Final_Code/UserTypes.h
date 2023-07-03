#ifndef UserTypes_h
#define UserTypes_h
#include "Arduino.h"
#define FILE_BASE_NAME "mpubmpraw"
struct data_t {
  int State;
  unsigned long time;
  int16_t ax;
  int16_t ay;
  int16_t az;
  int16_t gx;
  int16_t gy;
  int16_t gz;
  float altitude;
  double lat;
  double lng;
};
#endif  // UserTypes_h

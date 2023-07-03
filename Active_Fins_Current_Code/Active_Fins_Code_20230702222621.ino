#include <BufferedPrint.h>
#include <FreeStack.h>
#include <MinimumSerial.h>
#include <RingBuf.h>
#include <SdFat.h>
#include <SdFatConfig.h>
#include <sdios.h>
#include <MPU6050.h>  // 3 Axis Accelerometer and 3 Axis Gyroscope
#include <Adafruit_Sensor.h>
#include <TinyGPSPlus.h>
#include <Wire.h>  // allows communication with I2C/TWI devices (sensors)
#include <SPI.h>
#include <Adafruit_BMP280.h>  //Barometric Pressure and Altitude Sensor
#include <SoftwareSerial.h>
#include <EEPROM.h>
#include <string.h>
#include "IMU_Zero.h"
#include "UKF.h"
#include "UserTypes.h"
#include <BasicLinearAlgebra.h>
#include <ElementStorage.h>
#include <math.h>

// Example to demonstrate write latency for preallocated exFAT files.
// I suggest you write a PC program to convert very large bin files.
//
// The maximum data rate will depend on the quality of your SD,
// the size of the FIFO, and using dedicated SPI.

//------------------------------------------------------------------------------
// This example was designed for exFAT but will support FAT16/FAT32.
// Note: Uno will not support SD_FAT_TYPE = 3.
// SD_FAT_TYPE = 0 for SdFat/File as defined in SdFatConfig.h,
// 1 for FAT16/FAT32, 2 for exFAT, 3 for FAT16/FAT32 and exFAT.
#define SD_FAT_TYPE 2
//------------------------------------------------------------------------------
// Interval between data records in microseconds.
// Try 250 with Teensy 3.6, Due, or STM32.
// Try 2000 with AVR boards.
// Try 4000 with SAMD Zero boards.
const uint32_t LOG_INTERVAL_USEC = 150000;

// Set USE_RTC nonzero for file timestamps.
// RAM use will be marginal on Uno with RTClib.
// 0 - RTC not used
// 1 - DS1307
// 2 - DS3231
// 3 - PCF8523
#define USE_RTC 0
#if USE_RTC
#include "RTClib.h"
#endif  // USE_RTC

// LED to light if overruns occur.
#define ERROR_LED_PIN -1

/*
  Change the value of SD_CS_PIN if you are using SPI and
  your hardware does not use the default value, SS.
  Common values are:
  Arduino Ethernet shield: pin 4
  Sparkfun SD shield: pin 8
  Adafruit SD shields and modules: pin 10
*/

// SDCARD_SS_PIN is defined for the built-in SD on some boards.
#ifndef SDCARD_SS_PIN
const uint8_t SD_CS_PIN = SS;
#else   // SDCARD_SS_PIN
// Assume built-in SD is used.
const uint8_t SD_CS_PIN = SDCARD_SS_PIN;
#endif  // SDCARD_SS_PIN

// FIFO SIZE - 512 byte sectors.  Modify for your board.
#ifdef __AVR_ATmega328P__
// Use 512 bytes for 328 boards.
#define FIFO_SIZE_SECTORS 1
#elif defined(__AVR__)
// Use 2 KiB for other AVR boards.
#define FIFO_SIZE_SECTORS 4
#else  // __AVR_ATmega328P__
// Use 8 KiB for non-AVR boards.
#define FIFO_SIZE_SECTORS 16
#endif  // __AVR_ATmega328P__

// Preallocate 1GiB file.
const uint32_t PREALLOCATE_SIZE_MiB = 1024UL;

// Try max SPI clock for an SD. Reduce SPI_CLOCK if errors occur.
#define SPI_CLOCK SD_SCK_MHZ(50)

// Try to select the best SD card configuration.
#if HAS_SDIO_CLASS
#define SD_CONFIG SdioConfig(FIFO_SDIO)
#elif ENABLE_DEDICATED_SPI
#define SD_CONFIG SdSpiConfig(SD_CS_PIN, DEDICATED_SPI, SPI_CLOCK)
#else  // HAS_SDIO_CLASS
#define SD_CONFIG SdSpiConfig(SD_CS_PIN, SHARED_SPI, SPI_CLOCK)
#endif  // HAS_SDIO_CLASS

// Save SRAM if 328.
#ifdef __AVR_ATmega328P__
#include "MinimumSerial.h"
MinimumSerial MinSerial;
#define Serial MinSerial
#endif  // __AVR_ATmega328P__

const float seapress = 1014.22;  //pressure (hPa) at sea level in day in your region
const float act_alt = 43;

int16_t ax, ay, az, gx, gy, gz;

using namespace BLA;

ArrayMatrix<6, 1, double> x_aposteriori = {5, 0, 0, 0, 0, 0};
ArrayMatrix<3, 1, double> z;

int State_Machine = 1;
const int buzzer = 44;

//General
int K_P_r = -0.06;  //proportional control gain for roll
int K_D_r = -0.025;  //derivative control gain for roll

int K_P_y = 2;    //proportional control gain for yaw
int K_D_y = 2;  //derivative control gain for yaw

int K_P_p = 2;    //proportional control gain for pitch
int K_D_p = 2;  //derivative control gain for pitch

double roll_angle = 0;
double roll_angle_mvmt = 0;

double yaw_angle = 0;
double yaw_angle_mvmt = 0;

double pitch_angle = 0;
double pitch_angle_mvmt = 0;

long Servo_1_Pos = 0;
long Servo_2_Pos = 0;
long Servo_3_Pos = 0;
long Servo_4_Pos = 0;

long interval = 0;

int Servo_1_Pin = 5;
int Servo_2_Pin = 6;
int Servo_3_Pin = 7;
int Servo_4_Pin = 8;

unsigned long timePrev = 0;
double timeInt;
unsigned long timeStart = 0;
unsigned long Launch_Time = 0;
unsigned long Fin_Control_Start_Time;
long timeFreeFall;
unsigned long apogeeTime;
unsigned long unstableTime;
unsigned long parachuteTime;

float prevAltitude;
float maxAltitude = 0;


#define BMP_SCK (13)
#define BMP_MISO (12)
#define BMP_MOSI (11)
#define BMP_CS (10)  // barometic sensor chip select

Adafruit_BMP280 bmp;
Adafruit_Sensor* bmp_pressure = bmp.getPressureSensor();

MPU6050 mpu;

//==============================================================================
// Replace logRecord(), printRecord(), and ExFatLogger.h for your sensors.
void logRecord(data_t* data, uint16_t overrun) {
  if (overrun) {
    // Add one since this record has no adc data. Could add overrun field.
    overrun++;
    data->time = 0;
  } else {
      data->time = micros();
      mpu.getMotion6(&data->ax, &data->ay, &data->az,
                  &data->gx, &data->gy, &data->gz);
      data->altitude = bmp.readAltitude(seapress) + 6.25;

      timeInt = (data->time - timePrev) * 0.000001; //time interval in seconds

      if (mpu.getIntFreefallStatus() || (data->time - Launch_Time) > 30000000 && Launch_Time != 0) {
        State_Machine = 4;
      }

      if (maxAltitude > data->altitude) { //values are in pressure so the lowest pressure value will be the highest altitude
         maxAltitude = data->altitude; 
         apogeeTime = data->time;
      }

      if (((data->time - apogeeTime) > 500000 && apogeeTime != 0) || ((data->time - unstableTime) > 1000000 && unstableTime != 0)) {
        State_Machine = 4;
        parachuteTime = data->time;
        Serial.print("Apogee Time: ");
        Serial.println(apogeeTime);
        Serial.print("Parachute Time: ");
        Serial.println(parachuteTime);
      }
      

      
      switch (State_Machine) {
        case 1:                          //On Pad
          if (data->az < -2048) {  //if raw acceleration is greater than 2g
            State_Machine = 2;
            Launch_Time = data->time;
            Serial.print("Launch Time: ");
            Serial.println(Launch_Time);
          }
        break;

        case 2:  //After_Launch
          z = {int(data->gx) / 131, int(data->gy) / 131, int(data->gz) / 131 };
          UKF_full_step();
          if (data->altitude - act_alt > 30 || (data->time - Launch_Time) > 20000) {  //is our altitude greater than 30 meters different
            State_Machine = 3;
            Fin_Control_Start_Time = data->time;
            Serial.println(Fin_Control_Start_Time);
          }
        break;

        case 3:  //Fin Control
          {
            z = {data->gx / 131, data->gy / 131, data->gz / 131 };
            UKF_full_step();
            for (int i = 0; i < 6; i++) {
               Serial.print(x_aposteriori(i));
               Serial.print(", ");
            }
            roll_angle_mvmt = K_P_r * x_aposteriori(2) + K_D_r * x_aposteriori(5);  //get mvmt of servo based on roll position

            yaw_angle_mvmt = K_P_y * x_aposteriori(0) + K_D_y * x_aposteriori(3);  //get mvmt of servo based on yaw position

            pitch_angle_mvmt = K_P_p * x_aposteriori(1) + K_D_p * x_aposteriori(4);  //get mvmt of servo based on pitch position

            //servo 1 pos (Pin 5)
            OCR3A = round(16000000 * (0.0015 + 0.00001 * (roll_angle_mvmt + pitch_angle_mvmt)));  //add mvmt to determine overall position of servos
            //servo 2 pos (Pin 6)
            OCR4A = round(16000000 * (0.0015 + 0.00001 * (roll_angle_mvmt + yaw_angle_mvmt)));  //multiply mvmt by 10 microseconds (each degree on servo = 10 microseconds)
            //servo 3 pos (Pin 7)
            OCR4B = round(16000000 * (0.0015 + 0.00001 * (roll_angle_mvmt - pitch_angle_mvmt)));  //add that to 1500 microseconds (where the servo is at 90 degrees)
            //servo 4 pos (Pin 8)
            OCR4C = round(16000000 * (0.0015 + 0.00001 * (roll_angle_mvmt - yaw_angle_mvmt)));  //16000000*(0.0015 + 0.00001*(roll_angle_mvmt + yaw_angle_mvmt*cos(roll_angle) + pitch_angle_mvmt*sin(roll_angle))); //multiply that by 16 MHz to get the timer value necessary to get that servo angle
            
            if (OCR3A == 0) {
              OCR3A = 24000;
            } else if (OCR3A > 26400) {
              OCR3A = 26400;
            } else if (OCR3A < 21600) {
              OCR3A = 21600;
            }
            if (OCR4A == 0) {
              OCR4A = 24000;
            } else if (OCR4A > 26400) {
              OCR4A = 26400;
            } else if (OCR4A < 21600) {
              OCR4A = 21600;
            }
            if (OCR4B == 0) {
              OCR4B = 24000;
            } else if (OCR4B > 26400) {
              OCR4B = 26400;
            } else if (OCR4B < 21600) {
              OCR4B = 21600;
            }
            if (OCR4C == 0) {
              OCR4C = 24000;
            } else if (OCR4C > 26400) {
              OCR4C = 26400;
            } else if (OCR4C < 21600) {
              OCR4C = 21600;
            }

            
            Serial.print(OCR3A);
            Serial.print(", ");
            Serial.print(OCR4A);
            Serial.print(", ");
            Serial.print(OCR4B);
            Serial.print(", ");
            Serial.println(OCR4C);
            Serial.println(timeInt);
          }
          
        break;

        case 4:                                //Falling
          if (mpu.getIntZeroMotionStatus()) {  //
            State_Machine = 5;
          }
        break;

        case 5:  //On Ground
          {
            break;
          }
        break;
      }
      timePrev = data->time;

      prevAltitude = data->altitude;
    }
    Serial.print("State = ");
    Serial.print(State_Machine);
    Serial.println();
  }
//------------------------------------------------------------------------------
void printRecord(Print* pr, data_t* data) {
  static uint32_t nr = 0;
  if (!data) {
    pr->print(F("LOG_INTERVAL_USEC,"));
    pr->println(LOG_INTERVAL_USEC);
    pr->print(F("rec#"));
    pr->print(F("time, Acc X, Acc Y, Acc Z, Gyro X, Gyro Y, Gyro Z, Altitude"));
    pr->println();
    nr = 0;
    return;
  }
  if (data->time == 0 && data->ax == 0 && data->ay == 0 && data->az == 0 && data->gx == 0 && data->gy == 0 && data->gz == 0 && data->altitude == 0) {
     uint16_t n = data->time & 0X7FFFFFFF;
     nr += n;
     pr->print(F("-1,"));
     pr->print(n);
     pr->println(F(",overuns"));
   } else {
    pr->print(nr++);
      pr->write(',');
      pr->print(data->time);
      pr->write(',');
      pr->print(data->ax);
      pr->write(',');
      pr->print(data->ay);
      pr->write(',');
      pr->print(data->az);
      pr->write(',');
      pr->print(data->gx);
      pr->write(',');
      pr->print(data->gy);
      pr->write(',');
      pr->print(data->gz);
      pr->write(',');
      pr->println(data->altitude);
   }
}
//==============================================================================
const uint64_t PREALLOCATE_SIZE = (uint64_t)PREALLOCATE_SIZE_MiB << 20;
// Max length of file name including zero byte.
#define FILE_NAME_DIM 40
// Max number of records to buffer while SD is busy.
const size_t FIFO_DIM = 512 * FIFO_SIZE_SECTORS / sizeof(data_t);

#if SD_FAT_TYPE == 0
typedef SdFat sd_t;
typedef File file_t;
#elif SD_FAT_TYPE == 1
typedef SdFat32 sd_t;
typedef File32 file_t;
#elif SD_FAT_TYPE == 2
typedef SdExFat sd_t;
typedef ExFile file_t;
#elif SD_FAT_TYPE == 3
typedef SdFs sd_t;
typedef FsFile file_t;
#else  // SD_FAT_TYPE
#error Invalid SD_FAT_TYPE
#endif  // SD_FAT_TYPE

sd_t sd;

file_t binFile;
file_t csvFile;
// You may modify the filename.  Digits before the dot are file versions.
char binName[] = "ExFatLogger00.bin";
//------------------------------------------------------------------------------
#if USE_RTC
#if USE_RTC == 1
RTC_DS1307 rtc;
#elif USE_RTC == 2
RTC_DS3231 rtc;
#elif USE_RTC == 3
RTC_PCF8523 rtc;
#else  // USE_RTC == type
#error USE_RTC type not implemented.
#endif  // USE_RTC == type
// Call back for file timestamps.  Only called for file create and sync().
void dateTime(uint16_t* date, uint16_t* time, uint8_t* ms10) {
  DateTime now = rtc.now();

  // Return date using FS_DATE macro to format fields.
  *date = FS_DATE(now.year(), now.month(), now.day());

  // Return time using FS_TIME macro to format fields.
  *time = FS_TIME(now.hour(), now.minute(), now.second());

  // Return low time bits in units of 10 ms.
  *ms10 = now.second() & 1 ? 100 : 0;
}
#endif  // USE_RTC
//------------------------------------------------------------------------------
#define error(s) sd.errorHalt(&Serial, F(s))
#define dbgAssert(e) ((e) ? (void)0 : error("assert " #e))
//-----------------------------------------------------------------------------
// Convert binary file to csv file.
void binaryToCsv() {
  uint8_t lastPct = 0;
  uint32_t t0 = millis();
  data_t binData[FIFO_DIM];

  if (!binFile.seekSet(512)) {
    error("binFile.seek failed");
  }
  uint32_t tPct = millis();
  printRecord(&csvFile, nullptr);
  while (!Serial.available() && binFile.available()) {
    int nb = binFile.read(binData, sizeof(binData));
    if (nb <= 0) {
      error("read binFile failed");
    }
    size_t nr = nb / sizeof(data_t);
    for (size_t i = 0; i < nr; i++) {
      printRecord(&csvFile, &binData[i]);
    }

    if ((millis() - tPct) > 1000) {
      uint8_t pct = binFile.curPosition() / (binFile.fileSize() / 100);
      if (pct != lastPct) {
        tPct = millis();
        lastPct = pct;
        Serial.print(pct, DEC);
        Serial.println('%');
        csvFile.sync();
      }
    }
    if (Serial.available()) {
      break;
    }
  }
  csvFile.close();
  Serial.print(F("Done: "));
  Serial.print(0.001 * (millis() - t0));
  Serial.println(F(" Seconds"));
}
//------------------------------------------------------------------------------
void clearSerialInput() {
  uint32_t m = micros();
  do {
    if (Serial.read() >= 0) {
      m = micros();
    }
  } while (micros() - m < 10000);
}
//-------------------------------------------------------------------------------
void createBinFile() {
  binFile.close();
  while (sd.exists(binName)) {
    char* p = strchr(binName, '.');
    if (!p) {
      error("no dot in filename");
    }
    while (true) {
      p--;
      if (p < binName || *p < '0' || *p > '9') {
        error("Can't create file name");
      }
      if (p[0] != '9') {
        p[0]++;
        break;
      }
      p[0] = '0';
    }
  }
  if (!binFile.open(binName, O_RDWR | O_CREAT)) {
    error("open binName failed");
  }
  Serial.println(binName);
  if (!binFile.preAllocate(PREALLOCATE_SIZE)) {
    error("preAllocate failed");
  }

  Serial.print(F("preAllocated: "));
  Serial.print(PREALLOCATE_SIZE_MiB);
  Serial.println(F(" MiB"));
}
//-------------------------------------------------------------------------------
bool createCsvFile() {
  char csvName[FILE_NAME_DIM];
  if (!binFile.isOpen()) {
    Serial.println(F("No current binary file"));
    return false;  }

  // Create a new csvFile.
  binFile.getName(csvName, sizeof(csvName));
  char* dot = strchr(csvName, '.');
  if (!dot) {
    error("no dot in filename");
  }
  strcpy(dot + 1, "csv");
  if (!csvFile.open(csvName, O_WRONLY | O_CREAT | O_TRUNC)) {
    error("open csvFile failed");
  }
  clearSerialInput();
  Serial.print(F("Writing: "));
  Serial.print(csvName);
  Serial.println(F(" - type any character to stop"));
  return true;
}
//-------------------------------------------------------------------------------
void logData() {
  int32_t delta;  // Jitter in log time.
  int32_t maxDelta = 0;
  uint32_t maxLogMicros = 0;
  uint32_t maxWriteMicros = 0;
  size_t maxFifoUse = 0;
  size_t fifoCount = 0;
  size_t fifoHead = 0;
  size_t fifoTail = 0;
  uint16_t overrun = 0;
  uint16_t maxOverrun = 0;
  uint32_t totalOverrun = 0;
  uint32_t fifoBuf[128 * FIFO_SIZE_SECTORS];
  data_t* fifoData = (data_t*)fifoBuf;

  // Write dummy sector to start multi-block write.
  dbgAssert(sizeof(fifoBuf) >= 512);
  memset(fifoBuf, 0, sizeof(fifoBuf));
  if (binFile.write(fifoBuf, 512) != 512) {
    error("write first sector failed");
  }
  clearSerialInput();
  Serial.println(F("Type any character to stop"));

  // Wait until SD is not busy.
  while (sd.card()->isBusy()) {
  }

  // Start time for log file.
  uint32_t m = millis();

  // Time to log next record.
  uint32_t logTime = micros();
  while (true) {
    // Time for next data record.
    logTime += LOG_INTERVAL_USEC;

    // Wait until time to log data.
    delta = micros() - logTime;
    if (delta > 0) {
      Serial.print(F("delta: "));
      Serial.println(delta);
      error("Rate too fast");
    }
    while (delta < 0) {
      delta = micros() - logTime;
    }

    if (fifoCount < FIFO_DIM) {
      uint32_t m = micros();
      logRecord(fifoData + fifoHead, overrun);
      m = micros() - m;
      if (m > maxLogMicros) {
        maxLogMicros = m;
      }
      fifoHead = fifoHead < (FIFO_DIM - 1) ? fifoHead + 1 : 0;
      fifoCount++;
      if (overrun) {
        if (overrun > maxOverrun) {
          maxOverrun = overrun;
        }
        overrun = 0;
      }
    } else {
      totalOverrun++;
      overrun++;
      if (overrun > 0XFFF) {
        error("too many overruns");
      }
      if (ERROR_LED_PIN >= 0) {
        digitalWrite(ERROR_LED_PIN, HIGH);
      }
    }
    // Save max jitter.
    if (delta > maxDelta) {
      maxDelta = delta;
    }
    // Write data if SD is not busy.
    if (!sd.card()->isBusy()) {
      size_t nw = fifoHead > fifoTail ? fifoCount : FIFO_DIM - fifoTail;
      // Limit write time by not writing more than 512 bytes.
      const size_t MAX_WRITE = 512 / sizeof(data_t);
      if (nw > MAX_WRITE) nw = MAX_WRITE;
      size_t nb = nw * sizeof(data_t);
      uint32_t usec = micros();
      if (nb != binFile.write(fifoData + fifoTail, nb)) {
        error("write binFile failed");
      }
      usec = micros() - usec;
      if (usec > maxWriteMicros) {
        maxWriteMicros = usec;
      }
      fifoTail = (fifoTail + nw) < FIFO_DIM ? fifoTail + nw : 0;
      if (fifoCount > maxFifoUse) {
        maxFifoUse = fifoCount;
      }
      fifoCount -= nw;
      if (Serial.available()) {
        break;
      }
    }
  }
  Serial.print(F("\nLog time: "));
  Serial.print(0.001 * (millis() - m));
  Serial.println(F(" Seconds"));
  binFile.truncate();
  binFile.sync();
  Serial.print(("File size: "));
  // Warning cast used for print since fileSize is uint64_t.
  Serial.print((uint32_t)binFile.fileSize());
  Serial.println(F(" bytes"));
  Serial.print(F("totalOverrun: "));
  Serial.println(totalOverrun);
  Serial.print(F("FIFO_DIM: "));
  Serial.println(FIFO_DIM);
  Serial.print(F("maxFifoUse: "));
  Serial.println(maxFifoUse);
  Serial.print(F("maxLogMicros: "));
  Serial.println(maxLogMicros);
  Serial.print(F("maxWriteMicros: "));
  Serial.println(maxWriteMicros);
  Serial.print(F("Log interval: "));
  Serial.print(LOG_INTERVAL_USEC);
  Serial.print(F(" micros\nmaxDelta: "));
  Serial.print(maxDelta);
  Serial.println(F(" micros"));
}
//------------------------------------------------------------------------------
void openBinFile() {
  char name[FILE_NAME_DIM];
  clearSerialInput();
  Serial.println(F("Enter file name"));
  if (!serialReadLine(name, sizeof(name))) {
    return;
  }
  if (!sd.exists(name)) {
    Serial.println(name);
    Serial.println(F("File does not exist"));
    return;
  }
  binFile.close();
  if (!binFile.open(name, O_RDONLY)) {
    Serial.println(name);
    Serial.println(F("open failed"));
    return;
  }
  Serial.println(F("File opened"));
}
//-----------------------------------------------------------------------------
void printData() {
  if (!binFile.isOpen()) {
    Serial.println(F("No current binary file"));
    return;
  }
  // Skip first dummy sector.
  if (!binFile.seekSet(512)) {
    error("seek failed");
  }
  clearSerialInput();
  Serial.println(F("type any character to stop\n"));
  delay(1000);
  printRecord(&Serial, nullptr);
  while (binFile.available() && !Serial.available()) {
    data_t record;
    if (binFile.read(&record, sizeof(data_t)) != sizeof(data_t)) {
      error("read binFile failed");
    }
    printRecord(&Serial, &record);
  }
}
//------------------------------------------------------------------------------
void printUnusedStack() {
#if HAS_UNUSED_STACK
  Serial.print(F("\nUnused stack: "));
  Serial.println(UnusedStack());
#endif  // HAS_UNUSED_STACK
}
//------------------------------------------------------------------------------
bool serialReadLine(char* str, size_t size) {
  size_t n = 0;
  while (!Serial.available()) {
    yield();
  }
  while (true) {
    int c = Serial.read();
    if (c < ' ') break;
    str[n++] = c;
    if (n >= size) {
      Serial.println(F("input too long"));
      return false;
    }
    uint32_t m = millis();
    while (!Serial.available() && (millis() - m) < 100) {
    }
    if (!Serial.available()) break;
  }
  str[n] = 0;
  return true;
}
//------------------------------------------------------------------------------
void testSensor() {
  const uint32_t interval = 2000;
  int32_t diff;
  data_t data;
  clearSerialInput();
  Serial.println(F("\nTesting - type any character to stop\n"));
  delay(1000);
  printRecord(&Serial, nullptr);
  uint32_t m = micros();
  while (!Serial.available()) {
    m += interval;
    do {
      diff = m - micros();
    } while (diff > 0);
    logRecord(&data, 0);
    printRecord(&Serial, &data);
  }
}
//------------------------------------------------------------------------------
void setup() {
  if (ERROR_LED_PIN >= 0) {
    pinMode(ERROR_LED_PIN, OUTPUT);
    digitalWrite(ERROR_LED_PIN, HIGH);
  }
  Serial.begin(9600);

  // Wait for USB Serial
  while (!Serial) {
    yield();
  }
  delay(1000);
  FillStack();
#if !ENABLE_DEDICATED_SPI
  Serial.println(
      F("\nFor best performance edit SdFatConfig.h\n"
        "and set ENABLE_DEDICATED_SPI nonzero"));
#endif  // !ENABLE_DEDICATED_SPI

  Serial.print(FIFO_DIM);
  Serial.println(F(" FIFO entries will be used."));

  // Initialize SD.
  if (!sd.begin(SD_CONFIG)) {
    sd.initErrorHalt(&Serial);
  }
#if USE_RTC
  if (!rtc.begin()) {
    error("rtc.begin failed");
  }
  if (!rtc.isrunning()) {
    // Set RTC to sketch compile date & time.
    // rtc.adjust(DateTime(F(__DATE__), F(__TIME__)));
    error("RTC is NOT running!");
  }
  // Set callback
  FsDateTime::setCallback(dateTime);
#endif  // USE_RTC

  bmp.begin();
  //UKF Setup
    UKF_setup();

    //Buzzer Setup
    pinMode(buzzer, OUTPUT);

    tone(buzzer, 1000);
    delay(200);
    noTone(buzzer);

    delay(2000);

    mpu.initialize();

    delay(2000);
    

    IMU_Zero();  //Calibrate the sensors

    if (mpu.testConnection()) {
      //Setup MPU
      mpu.setFullScaleAccelRange(3);  //0 = +/- 2g, 1 4g, 2 8g, 3 = +/- 16g
      mpu.setFullScaleGyroRange(0);   //0 = +/- 250 deg/s, 1 500 degrees/s, 2 1000 deg/s, 3 2000 deg/s

      mpu.setFreefallDetectionThreshold(40);
      Serial.println(mpu.getFreefallDetectionThreshold());
      mpu.setFreefallDetectionDuration(100);
      Serial.println(mpu.getFreefallDetectionDuration());
      mpu.setFreefallDetectionCounterDecrement(1);
      Serial.println(mpu.getFreefallDetectionCounterDecrement());

      mpu.setZeroMotionDetectionThreshold(10);
      Serial.println(mpu.getFreefallDetectionThreshold());
      mpu.setZeroMotionDetectionDuration(255);
      Serial.println(mpu.getFreefallDetectionDuration());
      //Setup Pressure Sensor and Altimeter

      mpu.setDLPFMode(0);
      mpu.setDHPFMode(2);

      bmp.begin();

      bmp.setSampling(Adafruit_BMP280::MODE_NORMAL,      /* Operating Mode. */
                      Adafruit_BMP280::SAMPLING_X1,      /* Temp. oversampling */
                      Adafruit_BMP280::SAMPLING_X1,      /* Pressure oversampling */
                      Adafruit_BMP280::FILTER_OFF,       /* Filtering. */
                      Adafruit_BMP280::STANDBY_MS_1); /* Standby time. */

      //Send signal to buzzer to let know that the sensors and SD card are correctly connected
      tone(buzzer, 2000);
      delay(200);
      noTone(buzzer);
      tone(buzzer, 2000);
      delay(200);
      noTone(buzzer);
      delay(2000);
      Serial.println("");
    } else {
      tone(buzzer, 400);
      delay(200);
      noTone(buzzer);
      delay(2000);
      Serial.println("System initialization failed");
      //return;
    }
    mpu.getMotion6(&ax, &ay, &az, &gx, &gy, &gz);

    float off_acc_x = ax;     //if the calibration worked correctly these values should all be close to zero
    float off_acc_y = ay;     //if they aren't figure out what is wrong with the calibration process
    float off_acc_z = az - 2048;  //if they are somewhat close, but the function fails adjust values as you see fit
    float off_gyro_x = gx;
    float off_gyro_y = gy;
    float off_gyro_z = gz;
    float altitude_test = bmp.readAltitude(seapress) + 14;
    Serial.println(off_acc_x);
    Serial.println(off_acc_y);
    Serial.println(off_acc_z);
    Serial.println(off_gyro_x);
    Serial.println(off_gyro_y);
    Serial.println(off_gyro_z);
    Serial.println(altitude_test);
    if (abs(off_acc_x) < 50 && abs(off_acc_y) < 50 && abs(off_acc_z) < 50 && abs(off_gyro_x) < 50 && abs(off_gyro_y) < 50 && abs(off_gyro_z) < 50 && abs(altitude_test - act_alt) < 1) {
      Serial.println("Calibration Successful");
      tone(buzzer, 1500);
      delay(200);
      noTone(buzzer);
      delay(2000);
    } else {
      Serial.println("Calibration Failed");
      tone(buzzer, 400);
      delay(200);
      noTone(buzzer);
      delay(2000);
    }

    pinMode(Servo_1_Pin, OUTPUT);
    pinMode(Servo_2_Pin, OUTPUT);
    pinMode(Servo_3_Pin, OUTPUT);
    pinMode(Servo_4_Pin, OUTPUT);

    //Servo Setup
    /*
  //These are all bit values that determine various aspects of the Timer used

  //These bits control the output compare pins (1 0  Clear OCnA/OCnB/OCnC on compare match, set OCnA/OCnB/OCnC at BOTTOM (non-inverting mode)
  COM4A0 = 1;

  COM4B1 = 0;


/*These four bits control the wave form created (1 1 1 0 is FAST PVM with ICRn as the top)
  WGM43 = 1;

  WGM42 = 1;

  WGM41 = 1;

  WGM40 = 0;


  //These three bits control the prescaler values (0 0 1 is no prescaling)
  CS42 = 0;

  CS41 = 0;

  CS40 = 1;
  */

    TCCR4A = 1 << COM4A1 | 1 << COM4B1 | 1 << COM4C1 | 1 << WGM41 | 0 << WGM40;

    TCCR4B = 1 << WGM43 | 1 << WGM42 | 1 << CS40;

    ICR4 = 39999;
    OCR4A = 24000;  //Pin 6
    OCR4B = 24000;  //Pin 7
    OCR4C = 24000;  //Pin 8

    TCCR3A = 1<<COM3A1 | 1<<COM3B1 | 1<<COM3C1 | 1<<WGM31 | 0<<WGM30;

    TCCR3B = 1<<WGM33 | 1<<WGM32 | 1<<CS30;

    ICR3 = 39999;
    OCR3A = 24000;  //Pin 5

    tone(buzzer, 1200);  //Move Fin 1 back and forth 5 degrees
    OCR3A = 24800;
    delay(1000);
    tone(buzzer, 800);
    delay(1000);
    OCR3A = 23200;
    tone(buzzer, 1000);
    OCR3A = 24000;
    delay(1000);
    tone(buzzer, 1200);  //Move Fin 2 back and forth 5 degrees
    OCR4A = 24800;
    delay(1000);
    tone(buzzer, 800);
    delay(1000);
    OCR4A = 23200;
    tone(buzzer, 1000);
    OCR4A = 24000;
    delay(1000);
    tone(buzzer, 1200);  //Move Fin 3 back and forth 5 degrees
    OCR4B = 24800;
    delay(1000);
    tone(buzzer, 800);
    delay(1000);
    OCR4B = 23200;
    tone(buzzer, 1000);
    OCR4B = 24000;
    delay(1000);
    tone(buzzer, 1200);  //Move Fin 4 back and forth 5 degrees
    OCR4C = 24800;
    delay(1000);
    tone(buzzer, 800);
    delay(1000);
    OCR4C = 23200;
    tone(buzzer, 1000);
    OCR4C = 24000;
    delay(1000);
    noTone(buzzer);

    timePrev = micros();
    timeStart = micros();

   
}
//------------------------------------------------------------------------------
void loop() {
  printUnusedStack();
  // Read any Serial data.
  clearSerialInput();

  if (ERROR_LED_PIN >= 0) {
    digitalWrite(ERROR_LED_PIN, LOW);
  }
  Serial.println();
  Serial.println(F("type: "));
  Serial.println(F("b - open existing bin file"));
  Serial.println(F("c - convert file to csv"));
  Serial.println(F("l - list files"));
  Serial.println(F("p - print data to Serial"));
  Serial.println(F("r - record data"));
  Serial.println(F("t - test without logging"));
  // while (!Serial.available()) {
  //   yield();
  // }
  // char c = tolower(Serial.read());
  // Serial.println();

  // if (c == 'b') {
  //   openBinFile();
  // } else if (c == 'c') {
  //   if (createCsvFile()) {
  //     binaryToCsv();
  //   }
  // } else if (c == 'l') {
  //   Serial.println(F("ls:"));
  //   sd.ls(&Serial, LS_DATE | LS_SIZE);
  // } else if (c == 'p') {
  //   printData();
  // } else if (c == 'r') {
    createBinFile();
    logData();
//   } else if (c == 't') {
//     testSensor();
//   } else {
//     Serial.println(F("Invalid entry"));
//   }
}

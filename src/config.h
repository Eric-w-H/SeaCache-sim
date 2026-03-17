/*
Keep all hardware configurations.
*/

#ifndef CONFIG_H
#define CONFIG_H
#include "headers.h"

// global bandwdith & SRAM configuration
extern f64 HBMbandwidth;
extern int PEcnt, mergecnt;
extern f64 HBMbandwidthperPE;
extern int sramBank, sramReadPort, sramWritePort;

// DRAM ↔ SRAM bandwidth calculations
f64 memoryBandwidthWhole(long long ss);
f64 memoryBandwidthPE(long long ss);

// SRAM cycle
long long sramReadBandwidth(long long ss);
long long sramWriteBandwidth(long long ss);

extern bool ISCACHE;

#endif
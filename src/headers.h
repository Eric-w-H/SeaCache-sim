#ifndef HEADERS_H
#define HEADERS_H

#include <concepts>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <thread>
#include <vector>

#include <random>

typedef uint64_t    u64;
typedef uint32_t    u32;
typedef uint16_t    u16;
typedef uint8_t     u8;
typedef int64_t     i64;
typedef int32_t     i32;
typedef int16_t     i16;
typedef int8_t      i8;
typedef double      f64;
typedef float       f32;
typedef int32_t     b32;
typedef int16_t     b16;
typedef int8_t      b8;

typedef uint32_t    Coord;

template<typename A, typename B>
auto div_rup(A a, B b) -> decltype(a/b) { return (a + b - 1) / b; }

#include "config.h"
#include "data.h"
#include "estimation.h"
#include "parameters.h"
#include "simulator.h"
#include "util.h"

#endif

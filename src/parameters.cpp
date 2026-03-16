#include "parameters.h"

// smaple parameter p
double samplep = 0.01; // sampling probability
// sample parameter k
double samplek = 100;

int BSRsize = 32;

string printInterOrder[] = {"IJK", "IKJ", "JKI", "JIK", "KIJ", "KJI"};
string printDataFlow[] = {"Inner", "Outer", "Gust"};
string printFormat[] = {"RR", "RC", "CR", "CC", "BB"};

int Asize;
int Bsize;
int Csize;

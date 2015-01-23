#pragma once
#ifndef SVFINDER_H
#define SVFINDER_H


//#include <errno.h>
//#include <stdlib.h>

//#include <time.h>
//#include <math.h>
#include <iostream>
#include <stdio.h>
#include <string.h>

//#include "io.h"   // For access() in Windows.
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().

//#include "myFunc.h"

using namespace std ;
#include "ConfigFile.h"

#define MAX_ALPHABET 4
#define NA -1
#define FileNameLength 281 //maximal length of file name
#define MaxWordLength 281 // LLIB
#define Infinity 1000000000000
#define GOODSD 0.5
#define TELO_CENTRO_FLANCS 50000
//#define GOODMAPPABILITY 0.85

#define BALANCED 1
#define UNBALANCED 0
#define DIRECT 0
#define INVERTED 1
#define SAMECHR 1
#define DIFFCHR 0

#define MAXUncertainty 0.5

#define false 0
#define true 1


typedef int SYMB;
typedef long double ldouble;
extern int verbose;
extern char *ULetters[90];
extern int nProcesses;
extern double minMappabilityPerWindow;
extern bool uniqueMatch;

//extern double PI = 3.141592;

#include "ConfigFile.h"
#include "Chameleon.h"
#include "Genome.h"

#endif //SVFINDER_H


#ifndef HEADER_F465D40157E459E
#define HEADER_F465D40157E459E

#pragma once
#ifndef __MY_FUNC_H__
#define __MY_FUNC_H__

#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <algorithm>
//#include <cstdint> //for windows
#include <stdint.h>  //for Linux

#ifndef NA
#define NA -1
#endif

#define MAXDEGREE 10

#define SIMPLERIGHT 0
#define LARGECLOSE 1
#define NORMALLEVEL 2
#define HALFLENGTH 3
#define NOCALL 4

#define FLANKSEQ 50
#define MINCHRLENFORSV 1000

#define FASTASTRINGLEN 50
#define LEFT 0
#define RIGHT 1

#define HOMOZYG_MEAN 0.11
#define middleComponentMinWeight 0.1

 #define pi 3.14159
 #define ZERO 0.0000000001
 #define MAX_BUFFER 2048


#ifndef INFINITY
# define INFINITY 1000000000
#endif

class myFunc
{
public:
	myFunc(void);
	~myFunc(void);
};

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
float get_median(const std::vector<float>& data) ;
float get_median(const std::vector<float>& data, int start, int end) ;
float get_mean(const std::vector<float>& data) ;
float get_weighted_mean(const std::vector<float>& data, const std::vector<float>& weights) ;
float get_sum(const std::vector<float>& data) ;
float get_iqr(const std::vector<float>& data);
void readFileWithGenomeInfo(const std::string &chrLenFileName, std::vector<std::string>& names, std::vector<int>& lengths);
unsigned long sum(const std::vector<int>& data);
long getLineNumber(std::string const& file);
int factorial (int num);
int get_max_index(const std::vector<float>& data);
int get_min_index(const std::vector<float>& data);
void vector_sub (std::vector<float>& a, const std::vector<float>& b) ;
void findnextbreakpoint( const std::vector<float>&, int & besti, int & bestsign,float & bestlambda, int segLength, float a,float b);
int calculateBreakpoints_general(double threshold, int length, const std::vector<float>& ratio, std::vector<int>& bpfinal, int normal_chrom_length, int breakPointType);
void vector_scale(std::vector<float>& a, float const b);
void vector_add_constant(std::vector<float>& a, float const b) ;
float sd(std::vector<float>& a, float const b) ;
float round_by_ploidy(float value, int ploidy);
int argmin(const std::vector<double> & myvector);
std::string pathAppend(const std::string& p1, const std::string& p2);
int isCG (const char & a) ;
int isSpaceCharacter (const char & a) ;
int isN (const char & a) ;
float polynomial(const float x, const double a, const double b, const double c); //ax^2+bx+c
float polynomial(const float x, const double a0, const double a1, const double a2, const double a3);//a0x^3+a1x^2+a2x+a3
float polynomial(const float x, const double * a, double ratio, int degree); //any degree
int round(float r);
void processChrName(std::string & chr) ;
size_t strccnt(const char *s, int c);

void getBAFinfo(std::string BAFValuesInTheSegment,float copyNumber,float &estimatedBAF,float &fittedBAF, std::string &medianBAFSym,float &uncertainty, float normalContamination);
void getCopyNumbers (float copyNumber, std::vector <int> & copyNumbers);
double calculateLogLikelyHoodNormalMixtureForBAFs(std::vector <float> BAFs,std::vector <float> mu,float middleComponentMinW, bool isMuFixed);

double NormalDistributionDensity (double x, double mu, double sigma) ;

char complement(char nucleotide);
void myReplace(std::string& str, const std::string& oldStr, const std::string& newStr);

std::string findSmallestSuffix (std::string s1, std::string s2);
std::string getNormalBAFforPloidy(int ploidy);
std::string getXYBAFforPloidy(int ploidy);

bool getSAMinfo(std::string line,std::string &chr1,std::string &chr2,std::string &orient1,std::string &orient2,int &left,int &right);
bool getELANDinfo(std::string line,std::string &chr1,std::string &chr2,std::string &orient1,std::string &orient2,int &left,int &right,int &insertSize);

void fillSizes(std::vector <int> & insertSize,int number,int minLength,int maxLength);

std::string complementary(std::string & line);

int calculateTotalLength(std::vector <int> lefts,std::vector <int> rights);

std::string randomSeq(int insLength);

#ifdef _WIN32
double expm1(double x);
double log1p(double x);
#endif

//template <class T>
//T merge_no_dups(const T& v1, const T& v2);

template <class T>
T merge_no_dups(const T& v1, const T& v2) {
    T ret(v1);
    for (typename T::const_iterator i = v2.begin(); i != v2.end(); ++i) {
        bool notfound = true;
        for (typename T::const_iterator z = ret.begin(); z != ret.end(); ++z) {
            if (*z == *i) {
                notfound = false;
                break;
            }
        }
        if (notfound) {
            ret.push_back(*i);
        }
    }
    return ret;
}
//int round(double r);
#endif

#endif // header guard

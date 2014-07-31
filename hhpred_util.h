#pragma once
#include <iostream>   // cin, cout, cerr
#include <fstream>    // ofstream, ifstream
#include <cstdio>     // printf
#include <stdlib.h>   // exit
#include <time.h>     // clock
#include <math.h>     // sqrt, pow
#include <limits.h>   // INT_MIN
#include <float.h>    // FLT_MIN
#include <string.h>   // strcmp, strstr
#include <cassert>
#ifdef HH_SSE3
#include <emmintrin.h>
#include <pmmintrin.h>
#endif
using namespace std;


#define NAA 20
#define NTRANS 7
#define MAXRES 15002
#define LINELEN 524288
#define HMMSCALE 1000 
#define SELFEXCL 3
#define LAMDA 0.388
//------- ws_util -------//
enum transitions {M2M,M2I,M2D,I2M,I2I,D2M,D2D};
enum pair_states {STOP=0,SAME=1,GD=2,IM=3,DG=4,MI=5,MS=6,ML=7,SM=8,LM=9,MM=10};


/////////////////////////////////////////////////////////////////////////////////////
// Arithmetics
/////////////////////////////////////////////////////////////////////////////////////

//// max and min
inline double dmax(double x, double y) { return (x>y? x : y);}
inline double dmin(double x, double y) { return (x<y? x : y);}
inline int imax(int x, int y) { return (x>y? x : y);}
inline int imin(int x, int y) { return (x<y? x : y);}
inline int iabs(int x) { return (x>=0? x : -x);}

// Rounding up, rounding down and rounding to nearest integer
inline int iceil(double x)  {return int(ceil(x));}
inline int ifloor(double x) {return int(floor(x));}
inline int iround(double x) {return int(floor(x+0.5));}

//// Generalized mean: d=0: sqrt(x*y)  d=1: (x+y)/2  d->-inf: min(x,y)  d->+inf: max(x,y)
inline double fmean(double x, double y, double d) { return pow( (pow(x,d)+pow(y,d))/2 ,1./d);}
inline float frand() { return rand()/(RAND_MAX+1.0); }

// log base 2
inline float log2(float x)  {return (x<=0? (float)(-100000):1.442695041*log(x));}
inline float log10_(float x) {return (x<=0? (float)(-100000):0.434294481*log(x));}

//------- fast function -------//
extern float flog2(float x);
extern float fast_log2(float x);
extern float fast_log_gamma(float x);
extern float fpow2(float x);

//------- normalization -------//
extern double normalize_to_one(double* array, size_t length, const float* def_array=NULL);
extern float normalize_to_one(float* array, size_t length, const float* def_array=NULL);
extern float NormalizeTo1(float* array, int length, float* def_array=NULL);
extern float NormalizeToX(float* array, int length, float x, float* def_array=NULL);
extern char* sprintg(float val, int w);

//------- string function -----//
extern int strmcpy(char* dest, const char* source, size_t maxlen);
extern int strtoi(const char*& ptr);
extern int strtoi_(const char*& ptr, int deflt=INT_MAX);
extern int strint(char*& ptr);
extern int strinta(char*& ptr, int deflt=99999);
extern float strflt(char*& ptr);
extern float strflta(char*& ptr, float deflt=99999);
extern int chomp(char str[]);
extern char* fgetline(char str[], const int maxlen, FILE* file);
extern char *substr(char* substr, char* str, int a, int b);
extern char* strscn(char* str);
extern char* strscn_ws(char* str);
extern const char* strscn_c(const char* str);
extern char* strscn_(char* str);
extern char* strscn(char* str, const char c);
extern char* strscn_(char* str, const char c);
extern char* strcut(char* str);
extern char* strcut_(char* str);
extern char* strcut(char* str, const char c);
extern char* strcut_(char* str, const char c);
extern char* strcut(char* str, const char* substr);
extern char* strcut_(char* str, const char* substr);
extern char* strwrd(char* str, char* ptr);
extern char* strwrd(char* str, char* ptr, int maxlen);
extern char* strwrd(char* str, char* ptr, const char c);
extern int strtr(char* str, const char oldchars[], const char newchars[]);
extern int strtrd(char* str, const char chars[]);
extern int strtrd(char* str, char char1, char char2);
extern int strcount(char* str, char char1, char char2);
extern char* uprstr(char* str);
extern char* lwrstr(char* str);
extern char uprchr(char chr);
extern char lwrchr(char chr);
extern char* strsubst(char* str, const char str1[], const char str2[]);

//------ path related ------//
//extern void ElapsedTimeSinceFirstCall(const char str[]);
//extern void ElapsedTimeSinceLastCall(const char str[]);
extern char* RemovePath(char outname[], char filename[]);
extern char* RemoveExtension(char outname[], char filename[]);
extern char* RemovePathAndExtension(char outname[], char filename[]);
extern char* Extension(char extension[], char filename[]);
extern char* Pathname(char pathname[], char filename[]);

//------- qsort -------//
extern void swapi(int k[], int i, int j);
extern void QSortInt(int v[], int k[], int left, int right, int up=+1);
extern void QSortFloat(float v[], int k[], int left, int right, int up=+1);
extern void QSortDouble(double v[], int k[], int left, int right, int up=+1);

//------ others -----//
extern float fast_dot_product_single2(float* qi, float* tj,int num);
extern float ScalarProd20(float* qi, float* tj);





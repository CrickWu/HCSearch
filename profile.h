#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;


//---- definitions ----//
//secondary structure
#ifndef HELIX
#define HELIX	0
#endif
#ifndef SHEET
#define SHEET	1
#endif
#ifndef LOOP
#define LOOP	2
#endif
//solvent accessibility
#ifndef BURIED
#define BURIED		0
#endif
#ifndef INTERMEDIATE
#define INTERMEDIATE	1
#endif
#ifndef EXPOSED
#define EXPOSED		2
#endif

//---- HHpred related ----//
const int M_M = 0;
const int M_I = 1;
const int M_D = 2;
const int I_M = 3;
const int I_I = 4;
const int D_M = 5;
const int D_D = 6;
const int _NEFF = 7;
const int I_NEFF = 8;
const int D_NEFF = 9;



//============== New+Delete+Equal Array ==============//
//new
template <class A>
inline void NewArray2D_(A *** warray,int Narray1,int Narray2)
{
	*warray=new A*[Narray1];
	for(int i=0;i<Narray1;i++) *(*warray+i)=new A[Narray2];
};
template <class A>
inline void NewArray3D_(A **** warray,int Narray1,int Narray2,int Narray3)
{
	*warray=new A**[Narray1];
	for(int i=0;i<Narray1;i++) *(*warray+i)=new A*[Narray2];
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++)  *(*(*warray+i)+j)=new A[Narray3];                  
};
template <class A>
inline void NewArray4D_(A ***** warray,int Narray1,int Narray2,int Narray3,int Narray4)
{
	*warray=new A***[Narray1];
	for(int i=0;i<Narray1;i++) *(*warray+i)=new A**[Narray2];
	for(int i=0;i<Narray1;i++)for(int j=0;j<Narray2;j++) *(*(*warray+i)+j)=new A*[Narray3];
	for(int i=0;i<Narray1;i++)for(int j=0;j<Narray2;j++)for(int k=0;k<Narray3;k++) *(*(*(*warray+i)+j)+k)=new A[Narray4];
};
//delete
template <class A>
inline void DeleteArray2D_(A *** warray,int Narray)
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray) delete [] (*warray);
	(*warray)=NULL;
};
template <class A>
inline void DeleteArray3D_(A **** warray,int Narray1,int Narray2)
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++)  if(*(*(*warray+i)+j)) delete [] *(*(*warray+i)+j);
	for(int i=0;i<Narray1;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray1) delete [] (*warray);
	(*warray)=NULL;
};
template <class A>
inline void DeleteArray4D_(A ***** warray, int Narray1,int Narray2,int Narray3 )
{
	if((*warray)==NULL)return;
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++) for(int k=0;k<Narray3;k++) if(*(*(*(*warray+i)+j)+k)) delete []*(*(*(*warray+i)+j)+k);
	for(int i=0;i<Narray1;i++) for(int j=0;j<Narray2;j++) if(*(*(*warray+i)+j)) delete [] *(*(*warray+i)+j);
	for(int i=0;i<Narray1;i++) if(*(*warray+i)) delete [] *(*warray+i);
	if(Narray1) delete [] (*warray);
	(*warray)=NULL;
};
//equal
template <class A> 
void EqualArray_(A *out,A *in,int len)
{
	for(int i=0;i<len;i++) *(out+i)=*(in+i);
};
template <class A> 
void EqualArray2D_(A **out,A **in,int len1,int len2)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++) *(*(out+i)+j)=*(*(in+i)+j);
};
template <class A> 
void EqualArray3D_(A ***out,A ***in,int len1,int len2,int len3)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++)for(int k=0;k<len3;k++) *(*(*(out+i)+j)+k)=*(*(*(in+i)+j)+k);
};
template <class A> 
void EqualArray4D_(A ****out,A ****in,int len1,int len2,int len3,int len4)
{
	for(int i=0;i<len1;i++)for(int j=0;j<len2;j++)for(int k=0;k<len3;k++)for(int l=0;l<len4;l++) *(*(*(*(out+i)+j)+k)+l)=*(*(*(*(in+i)+j)+k)+l);
};


//------ required data -------//
//the conversion of AA to subScript
const int AA2SUB[26]={0,20,1,2,3,4,5,6,7,20,8,9,10,11,20,12,13,14,15,16,20,17,18,20,19,20};
//AA1Coding
const int AA1Coding[21]={0,4,3,6,13,7,8,9,11,10,12,2,14,5,1,15,16,19,17,18,20};
//AA3Coding
//char* const AA3Coding[26]={"ALA","XXX","CYS","ASP","GLU","PHE","GLY","HIS","ILE","XXX","LYS","LEU","MET","ASN","XXX","PRO","GLN","ARG","SER","THR","XXX","VAL","TRP","XXX","TYR","XXX"};


//========= class PROFILE ========//
//[note]: read in PSM,PSP and HMM
class PROFILE
{
    public:
	int length;
	short *residue;
	float **PSP;
	float **PSM;
	float **ProfHMM;
	float **ProfHMM_original;
	float **EmissionScore;
	float **EmissionScore_original;
	float **EmissionProb;
	float **EmissionProb_original;
	//misc
	int WantOriginal;
	int failure;
	//init
	void Profile_Create_Matrix(int length);
	void Profile_Delete_Matrix(int length);
	//function related
	void WS_Process_PSM(ifstream &fin,string &filename,int StartPos=0);
	void WS_Process_PSP(ifstream &fin,string &filename,int StartPos=0);
	void WS_Process_HMM(ifstream &fin,string &filename,int StartPos=0);
	//constructor
	PROFILE(void);
	~PROFILE(void);
};

#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#ifdef HH_SSE3
#include <emmintrin.h>
#include <pmmintrin.h>
#endif
#include "template.h"
#include "seq.h"
#include "profile.h"
using namespace std;


//====== class: CNFalign_Basic =====//
class CNFalign_Basic
{
public:
	//---- constructor and destructor -----//
	CNFalign_Basic(SEQUENCE* t, SEQUENCE* s,float DISO_THRES=0.0);  //we set default to 0.0
	CNFalign_Basic(TEMPLATE* t, SEQUENCE* s,float DISO_THRES=0.0);  //we set default to 0.0
	CNFalign_Basic();
	~CNFalign_Basic();
	void CNFalign_Basic_Pure(float DISO_THRES=0.0);
	void CNFalign_Basic_Norm(float DISO_THRES=0.0);
	

	//---- main data structure ----//
	int NORM_or_PURE;  // temp<->seq (1) or seq<->seq (0)
	TEMPLATE* temp;    // template
	SEQUENCE* tseq;    // tempsequ
	SEQUENCE* seq;     // sequence
	PROFILE* temp_p;   // temp_profile
	PROFILE* seq_p;    // seq_profile
	string tnam;       // temp name
	string snam;       // targ name
	int length_x;      // the length of template
	int length_y;      // the length of sequence
	
	//======== related to alignment score calculate ======//
	float fast_dot_product_single(float* qi, float* tj,int num);
	//-> sequence
	int AApair(int tPos, int sPos);
	int BLOSUM62_Calc(int tPos, int sPos);
	float MutationOf2Pos4(int tPos,int sPos);
	float MutationOf2Pos4_(int tPos,int sPos);
	float MutationOf2Pos4__(int tPos,int sPos);
	float MutationOf2Pos5(int tPos,int sPos);
	float MutationOf2Pos6(int tPos,int sPos);
	float MutationOf2Pos7(int tPos,int sPos);
	//-> temp_seq structure
	float SingleOf2Pos_ACC_our_10_42_norm(int tPos,int sPos);
	float SSOf2Pos2_norm(int tPos,int sPos);
	float SSOf2Pos3_norm(int tPos,int sPos);
	float SSOf2Pos4_norm(int tPos,int sPos);
	//-> seq_seq structure
	float SingleOf2Pos_ACC_our_10_42_pure(int tPos,int sPos);
	float SSOf2Pos2_pure(int tPos,int sPos);
	float SSOf2Pos3_pure(int tPos,int sPos);
	float SSOf2Pos4_pure(int tPos,int sPos);
};


//========= extern basic data ==========//
//-> blosum matrix related 
extern const char ThreeLetterOrder[21];
extern const int Ori_BLOSUM_62[21][21];
extern const int Blo_AA_Map[21];
extern const int Ori_AA_Map[26];
//-> other substitution matrix
extern const float HDSM[400];
extern const float CC50[400];
extern const float subm[3][3];
extern const float HMMNull_f[21];
extern const float PSM_aver[21];
//-> for Mutation4 
extern const float log_20;
extern const float inv_log_20;
extern const float WS_Pow2_Null[20];
extern const float WS_Pow4_Null[20];
//-> ACC cutoff
extern const float BuryCore1;
extern const float BuryCore2;
//-> structure background
extern const float SS8_freq[8];
extern const float SS3_freq[3];
extern const float ACC_freq[3];

//======== extern temporary data =======//
extern float mut4_vector1[200];
extern float mut4_vector2[200];

#include "CNFalign_Basic.h"



//------------- BLOSUM matrix --------------//
const char ThreeLetterOrder[21] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','Z'};

const int Ori_BLOSUM_62[21][21]={
{  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -5 },  //A
{ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -5 },  //R
{ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -5 },  //N
{ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -5 },  //D
{  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -5 },  //C
{ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -5 },  //Q
{ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2, -5 },  //E
{  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -5 },  //G
{ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -5 },  //H
{ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -5 },  //I
{ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -5 },  //L
{ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -5 },  //K
{ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -5 },  //M
{ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -5 },  //F
{ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -5 },  //P
{  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -5 },  //S
{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -5 },  //T
{ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -5 },  //W
{ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -5 },  //Y
{  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -5 },  //V
{ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 }}; //Z
// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z

const int Blo_AA_Map[21]=
{ 0,19, 4, 3, 6, 13,7, 8, 9, 17,11,10,12,2, 18,14,5, 1, 15,16,20};
//A  V  C  D  E  F  G  H  I  W  K  L  M  N  Y  P  Q  R   S  T  Z
//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17  18 19 20
//Ori_Mapping//-----------------AVCDEFGHIWKLMNYPQRSTZ
const int Ori_AA_Map[26]=
{ 0,20,2,3,4,5,6,7,8,20,10,11,12,13,20,15,16,17,18,19,20, 1, 9,20,14,20};




//------------- two important matrix ----------------//
const float HDSM[400]={
  2.09,  -0.50,  -0.57,  -0.73,   0.33,  -0.75,  -0.12,   0.27,  -1.42,  -0.97,  -0.39,  -0.38,  -0.04,  -0.76,  -0.53,   0.34,   0.13,  -0.66,  -1.25,   0.02, //A
 -0.50,   2.87,   0.60,   0.13,  -1.30,   0.13,   0.99,  -0.96,   0.54,  -1.40,  -1.19,   1.42,  -0.63,  -1.40,   0.21,  -0.06,  -0.15,  -0.04,  -0.75,  -1.52, //R
 -0.57,   0.60,   3.60,   1.78,  -2.08,   0.33,  -0.16,   0.79,   0.76,  -2.43,  -2.10,   0.83,  -2.01,  -2.25,  -1.10,   0.40,   0.30,  -2.89,  -0.36,  -2.17, //N
 -0.73,   0.13,   1.78,   4.02,  -2.51,   0.34,   1.20,  -1.20,  -0.01,  -2.77,  -2.65,   0.66,  -2.58,  -2.19,   0.72,   0.71,  -0.75,  -1.91,  -1.21,  -2.02, //D
  0.33,  -1.30,  -2.08,  -2.51,   6.99,  -0.83,  -1.97,  -2.11,  -1.50,   0.13,  -0.31,  -2.19,   1.04,   1.13,  -2.19,   0.31,  -0.59,  -0.76,   0.13,   0.34, //C
 -0.75,   0.13,   0.33,   0.34,  -0.83,   2.60,   1.23,  -0.12,  -0.46,  -1.47,  -1.49,   0.92,  -0.13,  -2.31,   0.24,   1.04,   0.60,  -0.81,  -0.61,  -1.38, //Q
 -0.12,   0.99,  -0.16,   1.20,  -1.97,   1.23,   2.97,  -0.41,  -0.62,  -1.81,  -2.11,   1.11,  -1.86,  -1.61,  -0.26,   0.31,  -0.21,  -2.70,  -1.64,  -1.84, //E
  0.27,  -0.96,   0.79,  -1.20,  -2.11,  -0.12,  -0.41,   4.36,  -0.40,  -2.93,  -1.98,  -0.71,  -1.86,  -2.67,  -0.04,   0.29,  -0.81,  -1.21,  -1.62,  -1.96, //G
 -1.42,   0.54,   0.76,  -0.01,  -1.50,  -0.46,  -0.62,  -0.40,   5.89,  -1.76,  -0.93,   0.31,  -1.04,  -0.22,  -1.44,  -0.74,  -0.52,  -1.48,  -0.12,  -0.35, //H
 -0.97,  -1.40,  -2.43,  -2.77,   0.13,  -1.47,  -1.81,  -2.93,  -1.76,   2.76,   1.56,  -1.81,   0.99,   0.76,  -2.00,  -1.75,  -0.96,   0.25,   0.08,   1.94, //I
 -0.39,  -1.19,  -2.10,  -2.65,  -0.31,  -1.49,  -2.11,  -1.98,  -0.93,   1.56,   2.43,  -1.96,   1.61,   1.23,  -1.56,  -2.30,  -0.86,  -0.14,   0.70,   0.81, //L
 -0.38,   1.42,   0.83,   0.66,  -2.19,   0.92,   1.11,  -0.71,   0.31,  -1.81,  -1.96,   2.91,  -1.62,  -2.41,  -0.19,  -0.06,  -0.10,  -1.94,  -1.72,  -1.27, //K
 -0.04,  -0.63,  -2.01,  -2.58,   1.04,  -0.13,  -1.86,  -1.86,  -1.04,   0.99,   1.61,  -1.62,   3.75,   0.80,  -1.09,  -1.34,  -1.58,   0.87,  -0.41,   0.61, //M
 -0.76,  -1.40,  -2.25,  -2.19,   1.13,  -2.31,  -1.61,  -2.67,  -0.22,   0.76,   1.23,  -2.41,   0.80,   3.28,  -0.91,  -1.11,  -0.69,   2.29,   1.96,   0.51, //F
 -0.53,   0.21,  -1.10,   0.72,  -2.19,   0.24,  -0.26,  -0.04,  -1.44,  -2.00,  -1.56,  -0.19,  -1.09,  -0.91,   5.45,  -0.29,   0.93,  -5.34,  -1.98,  -1.11, //P
  0.34,  -0.06,   0.40,   0.71,   0.31,   1.04,   0.31,   0.29,  -0.74,  -1.75,  -2.30,  -0.06,  -1.34,  -1.11,  -0.29,   2.36,   1.20,  -1.18,  -1.56,  -1.11, //S
  0.13,  -0.15,   0.30,  -0.75,  -0.59,   0.60,  -0.21,  -0.81,  -0.52,  -0.96,  -0.86,  -0.10,  -1.58,  -0.69,   0.93,   1.20,   2.04,  -0.57,  -0.41,   0.05, //T
 -0.66,  -0.04,  -2.89,  -1.91,  -0.76,  -0.81,  -2.70,  -1.21,  -1.48,   0.25,  -0.14,  -1.94,   0.87,   2.29,  -5.34,  -1.18,  -0.57,   6.96,   2.15,  -1.09, //W
 -1.25,  -0.75,  -0.36,  -1.21,   0.13,  -0.61,  -1.64,  -1.62,  -0.12,   0.08,   0.70,  -1.72,  -0.41,   1.96,  -1.98,  -1.56,  -0.41,   2.15,   3.95,   0.21, //Y
  0.02,  -1.52,  -2.17,  -2.02,   0.34,  -1.38,  -1.84,  -1.96,  -0.35,   1.94,   0.81,  -1.27,   0.61,   0.51,  -1.11,  -1.11,   0.05,  -1.09,   0.21,   2.05};//V
//  A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W       Y       V

const float CC50[400]={
1.000,0.620,0.257,0.133,0.411,0.684,0.681,0.368,0.361,0.481,0.728,0.586,0.735,0.488,0.431,0.355,0.440,0.475,0.479,0.484, //A
0.620,1.000,0.407,0.283,0.303,0.862,0.630,0.410,0.679,0.203,0.525,0.929,0.739,0.683,0.399,0.689,0.655,0.659,0.594,0.208, //R
0.257,0.407,1.000,0.872,0.325,0.465,0.404,0.463,0.692,0.146,0.134,0.410,0.284,0.197,0.411,0.790,0.403,0.320,0.315,0.155, //N
0.133,0.283,0.872,1.000,0.215,0.383,0.466,0.362,0.674,0.034,0.026,0.535,0.155,0.076,0.302,0.707,0.300,0.202,0.201,0.046, //D
0.411,0.303,0.325,0.215,1.000,0.389,0.255,0.414,0.861,0.365,0.356,0.275,0.428,0.651,0.426,0.825,0.742,0.696,0.673,0.368, //C
0.684,0.862,0.465,0.383,0.389,1.000,0.740,0.477,0.711,0.270,0.710,0.859,0.626,0.640,0.465,0.475,0.711,0.657,0.655,0.278, //Q
0.681,0.630,0.404,0.466,0.255,0.740,1.000,0.370,0.362,0.157,0.401,0.743,0.245,0.184,0.363,0.396,0.359,0.288,0.297,0.168, //E
0.368,0.410,0.463,0.362,0.414,0.477,0.370,1.000,0.468,0.282,0.273,0.394,0.396,0.333,0.465,0.470,0.455,0.423,0.412,0.288, //G
0.361,0.679,0.692,0.674,0.861,0.711,0.362,0.468,1.000,0.276,0.265,0.416,0.387,0.567,0.454,0.930,0.704,0.661,0.738,0.285, //H
0.481,0.203,0.146,0.034,0.365,0.270,0.157,0.282,0.276,1.000,0.499,0.168,0.465,0.491,0.364,0.263,0.380,0.443,0.449,0.948, //I
0.728,0.525,0.134,0.026,0.356,0.710,0.401,0.273,0.265,0.499,1.000,0.606,0.791,0.741,0.358,0.252,0.371,0.688,0.447,0.498, //L
0.586,0.929,0.410,0.535,0.275,0.859,0.743,0.394,0.416,0.168,0.606,1.000,0.512,0.444,0.363,0.687,0.643,0.292,0.307,0.174, //K
0.735,0.739,0.284,0.155,0.428,0.626,0.245,0.396,0.387,0.465,0.791,0.512,1.000,0.815,0.445,0.366,0.438,0.821,0.487,0.466, //M
0.488,0.683,0.197,0.076,0.651,0.640,0.184,0.333,0.567,0.491,0.741,0.444,0.815,1.000,0.400,0.550,0.650,0.918,0.919,0.492, //F
0.431,0.399,0.411,0.302,0.426,0.465,0.363,0.465,0.454,0.364,0.358,0.363,0.445,0.400,1.000,0.449,0.463,0.470,0.470,0.370, //P
0.355,0.689,0.790,0.707,0.825,0.475,0.396,0.470,0.930,0.263,0.252,0.687,0.366,0.550,0.449,1.000,0.791,0.646,0.724,0.271, //S
0.440,0.655,0.403,0.300,0.742,0.711,0.359,0.455,0.704,0.380,0.371,0.643,0.438,0.650,0.463,0.791,1.000,0.699,0.904,0.387, //T
0.475,0.659,0.320,0.202,0.696,0.657,0.288,0.423,0.661,0.443,0.688,0.292,0.821,0.918,0.470,0.646,0.699,1.000,0.742,0.444, //W
0.479,0.594,0.315,0.201,0.673,0.655,0.297,0.412,0.738,0.449,0.447,0.307,0.487,0.919,0.470,0.724,0.904,0.742,1.000,0.452, //Y
0.484,0.208,0.155,0.046,0.368,0.278,0.168,0.288,0.285,0.948,0.498,0.174,0.466,0.492,0.370,0.271,0.387,0.444,0.452,1.000};//V
//A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V

const float subm[3][3]={{3.23,-22.57,-5.28}, {-22.57,   4.56,  -2.52}, { -5.28,  -2.52,   4.24}};

//---- constant for calculate log score ----//
const float HMMNull_f[21]={3.706,5.728,4.211,4.064,4.839,3.729,4.763,4.308,4.069,3.323,5.509,4.64,4.464,4.937,4.285,4.423,3.815,3.783,6.325,4.665,0};
const float PSM_aver[21]={6.4661,13.0051,11.824,13.9727,25.1246,11.04,11.2611,16.757,14.2443,13.6784,14.1105,10.9099,12.1047,17.5259,19.3949,7.7492,7.6928,27.4168,15.3683,10.9346,0};

//---- for Mutation4 ---//__added by WS__//
const float log_20 = log(20.0);	
const float inv_log_20 = 1.0/log(20.0);
const float WS_Pow2_Null[20]={13.050200,53.002922,18.519843,16.725762,28.620957,13.259918,27.152253,19.807845,16.783829,10.007433,45.538031,24.933267,22.069775,30.632687,19.494564,21.451401,14.074385,13.765642,80.170524,25.369092};
const float WS_Pow4_Null[20]={170.307717,2809.309775,342.984602,279.751101,819.159171,175.825438,737.244835,392.350710,281.696924,100.148708,2073.712245,621.667781,487.074957,938.361492,380.038034,460.162625,198.088319,189.492902,6427.312930,643.590848};

//---- constant for solvent accessibility ----//
const float BuryCore1=10.25; //the first cutoff for the solvent access
const float BuryCore2=42.9; //the second cuttoff for the solvent access

//------- background distribution --------//
const float SS8_freq[8]={0.30,0.045,0.005,0.19,0.01,0.125,0.125,0.20}; //H,G,I,E,B,T,S,L
const float SS3_freq[3]={0.35,0.20,0.45};     //H,E,C
const float ACC_freq[3]={0.333,0.333,0.333};  //B,M,E


//------- temp structure -----//
float mut4_vector1[200];
float mut4_vector2[200];



//======================== constructor and destructor ========================//
//-> default value for DISO_THRES is 0.0
CNFalign_Basic::CNFalign_Basic(SEQUENCE* t, SEQUENCE* s,float DISO_THRES) 
{
	//--- data init ---//
	tseq=t;
	seq=s;
	temp_p=(PROFILE*)t;
	seq_p=(PROFILE*)s;
	length_x = t->length;
	length_y = s->length;
	tnam=t->seq_name;
	snam=s->seq_name;
	//--- basic set ---//
	CNFalign_Basic_Pure(DISO_THRES);
}
CNFalign_Basic::CNFalign_Basic(TEMPLATE* t, SEQUENCE* s,float DISO_THRES) 
{
	//--- data init ---//
	temp=t;
	seq=s;
	temp_p=(PROFILE*)t;
	seq_p=(PROFILE*)s;
	length_x = t->length;
	length_y = s->length;
	tnam=t->temp_name;
	snam=s->seq_name;
	//--- basic set ---//
	CNFalign_Basic_Norm(DISO_THRES);
}
CNFalign_Basic::CNFalign_Basic()
{
}
CNFalign_Basic::~CNFalign_Basic()
{
}

//----- basic set -----//
void CNFalign_Basic::CNFalign_Basic_Pure(float DISO_THRES)
{
	NORM_or_PURE=0;   //pure
	SEQUENCE* t=tseq;
	SEQUENCE* s=seq;
	//--- recalculate probability ---//
	{
		//we set background for ss8, ss3 and acc
		//for template, we mark those missing and bad residue with background probability
		for(int i=0;i<t->length;i++)
		{
			float alpha=t->DISO[i];
			if(alpha>DISO_THRES) //default: 0.2
			{
				for(int k=0;k<8;k++)t->SS8_[i][k]=0;
				for(int k=0;k<3;k++)t->SS2_[i][k]=0;
				for(int k=0;k<3;k++)t->acc_our_10_42_[i][k]=0;
			}
			else
			{
				for(int k=0;k<8;k++)t->SS8_[i][k]=t->SS8[i][k];
				for(int k=0;k<3;k++)t->SS2_[i][k]=t->SS2[i][k];
				for(int k=0;k<3;k++)t->acc_our_10_42_[i][k]=t->acc_our_10_42[i][k];
			}
		}
		//for target, we mark each residue according to the disorder region
		for(int i=0;i<s->length;i++)
		{
			float alpha=s->DISO[i];
			if(alpha>DISO_THRES) //default: 0.2
			{
				for(int k=0;k<8;k++)s->SS8_[i][k]=0;
				for(int k=0;k<3;k++)s->SS2_[i][k]=0;
				for(int k=0;k<3;k++)s->acc_our_10_42_[i][k]=0;
			}
			else
			{
				for(int k=0;k<8;k++)s->SS8_[i][k]=s->SS8[i][k];
				for(int k=0;k<3;k++)s->SS2_[i][k]=s->SS2[i][k];
				for(int k=0;k<3;k++)s->acc_our_10_42_[i][k]=s->acc_our_10_42[i][k];
			}
		}
	}
}
void CNFalign_Basic::CNFalign_Basic_Norm(float DISO_THRES)
{
	NORM_or_PURE=1;   //norm
	TEMPLATE* t=temp;
	SEQUENCE* s=seq;
	//--- recalculate probability ---//
	{
		//we set background for ss8, ss3 and acc
		//for template, we mark those missing and bad residue with background probability
		for(int i=0;i<t->length;i++)
		{
			if(t->isMissing[i]==1 || t->residue[i]==20)
			{
				for(int k=0;k<8;k++)t->SS8_[i][k]=0;
				for(int k=0;k<3;k++)t->SS2_[i][k]=0;
				for(int k=0;k<3;k++)t->acc_[i][k]=0;
			}
			else
			{
				for(int k=0;k<8;k++)t->SS8_[i][k]=t->SS8[i][k];
				for(int k=0;k<3;k++)t->SS2_[i][k]=t->SS2[i][k];
				for(int k=0;k<3;k++)t->acc_[i][k]=t->acc[i][k];
			}
		}
		//for target, we mark each residue according to the disorder region
		for(int i=0;i<s->length;i++)
		{
			float alpha=s->DISO[i];
			if(alpha>DISO_THRES) //default: 0.2
			{
				for(int k=0;k<8;k++)s->SS8_[i][k]=0;
				for(int k=0;k<3;k++)s->SS2_[i][k]=0;
				for(int k=0;k<3;k++)s->acc_our_10_42_[i][k]=0;
			}
			else
			{
				for(int k=0;k<8;k++)s->SS8_[i][k]=s->SS8[i][k];
				for(int k=0;k<3;k++)s->SS2_[i][k]=s->SS2[i][k];
				for(int k=0;k<3;k++)s->acc_our_10_42_[i][k]=s->acc_our_10_42[i][k];
			}
		}
	}
}


/////===============   vice calculate part =========////
//---- fast dot ----//
float CNFalign_Basic::fast_dot_product_single(float* qi, float* tj,int num)
{
#ifdef HH_SSE3

	if(num>=8)
	{
		float __attribute__((aligned(16))) res;
		__m128 P; // query 128bit SSE2 register holding 4 floats
		__m128 R; // result
		__m128* Qi = (__m128*) qi;
		__m128* Tj = (__m128*) tj;
		int iter=(int)(num/4)-2;
		//head
		R = _mm_mul_ps(*(Qi++),*(Tj++));  //first 4
		//iter
		for(int i=0;i<iter;i++)   //mid 4*
		{
			P = _mm_mul_ps(*(Qi++),*(Tj++)); //1
			R = _mm_add_ps(R,P);  //1
		}
		//tail
		P = _mm_mul_ps(*Qi,*Tj);  //last 4
		R = _mm_add_ps(R,P);
		R = _mm_hadd_ps(R,R);
		R = _mm_hadd_ps(R,R);
		_mm_store_ss(&res, R);
	
		//--- final remain ---//
		return res;
	}
	else
	{
		float ws2=0;
		for(int i=0;i<num;i++)ws2+=qi[i]*tj[i];
		return ws2;
	}
	
#else

	float ws2=0;
	for(int i=0;i<num;i++)ws2+=qi[i]*tj[i];
	return ws2;
	
#endif
}





//======================= alignment score related =====================//
//------ pair calculation -----//
int CNFalign_Basic::AApair(int tPos, int sPos)
{
	int ta = temp_p->residue[tPos];
	int sa = seq_p->residue[sPos];
	if(ta==20 || sa==20)return -1;
	return ta*20+sa;
}
int CNFalign_Basic::BLOSUM62_Calc(int tPos, int sPos)
{
	int ta=temp_p->residue[tPos];
	int sa=seq_p->residue[sPos];
	if(ta==20 || sa==20)return 0;
	char a=ThreeLetterOrder[ta];
	char b=ThreeLetterOrder[sa];
	int ii,jj;
	if(a<'A' || a>'Z') {
		return 0; //a non-standard amino acid, neither reward nor penalize
	}
	ii=Blo_AA_Map[Ori_AA_Map[a-'A']];
	if(b<'A' || b>'Z') {
		return 0; //a non-standard amino acid
	}
	jj=Blo_AA_Map[Ori_AA_Map[b-'A']];
	return Ori_BLOSUM_62[ii][jj];
}

//--------- mutation score --------//
float CNFalign_Basic::MutationOf2Pos4(int tPos,int sPos)
{
	if(tPos<0 || sPos<0) return 0;
	if(tPos>=length_x || sPos >=length_y ) return 0;
	for(int i=0;i<20;i++)
	{
		mut4_vector1[i]=seq_p->EmissionProb[sPos][i];
		mut4_vector2[i]=temp_p->EmissionProb[tPos][i]*WS_Pow2_Null[i];
	}
	float m = fast_dot_product_single(mut4_vector1, mut4_vector2, 20);
	if(m>0) return log(m); 
	return -10;
}
float CNFalign_Basic::MutationOf2Pos4_(int tPos,int sPos)
{
	for(int i=0;i<20;i++)
	{
		mut4_vector1[i]=seq_p->EmissionProb[sPos][i];
		mut4_vector2[i]=temp_p->EmissionProb[tPos][i]*WS_Pow4_Null[i];
	}
	float m = fast_dot_product_single(mut4_vector1, mut4_vector2,20);
	if(m>0) return log(m);
	return 0;
}
float CNFalign_Basic::MutationOf2Pos4__(int tPos,int sPos)
{
	for (int i = 0; i < 20; i++) 
	{
		double sco=0.5 * (temp_p->EmissionProb[tPos][i] + seq_p->EmissionProb[sPos][i]);
		mut4_vector1[i]= sco;
		if(sco>0)mut4_vector2[i]= log(mut4_vector1[i]);
		else mut4_vector2[i]=0;
	}	
	float C = fast_dot_product_single(mut4_vector1, mut4_vector2,20);
	return (inv_log_20 * (C + log_20));
}
float CNFalign_Basic::MutationOf2Pos5(int tPos,int sPos)
{
	//--- seq_profile --//
	for(int i=0;i<20;i++)
	{
		mut4_vector1[i]=seq_p->EmissionProb[sPos][i];
		mut4_vector2[i]=temp_p->EmissionScore[tPos][i]+HMMNull_f[i];
	}
	float m1 = fast_dot_product_single(mut4_vector1, mut4_vector2,20);
	//--- profile_seq --//
	for(int i=0;i<20;i++)
	{
		mut4_vector1[i]=temp_p->EmissionProb[tPos][i];
		mut4_vector2[i]=seq_p->EmissionScore[sPos][i]+HMMNull_f[i];
	}
	float m2 = fast_dot_product_single(mut4_vector1, mut4_vector2,20);
	return (m1+m2);
}
float CNFalign_Basic::MutationOf2Pos6(int tPos,int sPos)
{
	int sAA=seq_p->residue[sPos];  //the sequence residue at sPos
	int y=AA2SUB[ThreeLetterOrder[sAA]-'A'];
	double m=0;
	if (y>=0 && y<20)
	{
		m+= temp_p->EmissionScore[tPos][y]+HMMNull_f[y]; // score between the template profile and sequence residue
	}
	return m;
}
float CNFalign_Basic::MutationOf2Pos7(int tPos,int sPos)
{
	int tAA=temp_p->residue[tPos]; //the template residue at tPos
	int x=AA2SUB[ThreeLetterOrder[tAA]-'A'];
	double m=0;
	if (x<20 && x>=0)
	{
		m+=seq_p->EmissionScore[sPos][x]+HMMNull_f[x]; //score between the sequence profile and template residue
	}
	return m;
}

//------- structure related ---------//norm
float CNFalign_Basic::SingleOf2Pos_ACC_our_10_42_norm(int tPos,int sPos)
{
	int t_acc=temp->ACC[tPos];
	int s_acc=seq->ACC[sPos];
	short diff = t_acc-s_acc>0?t_acc-s_acc:s_acc-t_acc;
	switch(diff)
	{
		case 0: return 3;
		case 1: return 1;
		case 2: return -2;
		default:
		{
			fprintf(stderr,"bad switch acc_10_42 here !!\n");
			exit(-1);
		}
	}
}
float CNFalign_Basic::SSOf2Pos2_norm(int tPos,int sPos)
{
	return (temp->SS[tPos]==seq->SS[sPos]);
}
float CNFalign_Basic::SSOf2Pos3_norm(int tPos,int sPos)
{
	int pred_type = seq->SS[sPos];
	float score = 0;
	for(int i=0;i<3;i++)
		score += subm[pred_type][i]*seq->SS2_[sPos][i];	
	return score;
}
float CNFalign_Basic::SSOf2Pos4_norm(int tPos,int sPos)
{
	float score = 0;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			score += temp->SS2[tPos][i] * seq->SS2[sPos][j] * subm[i][j];
	return score;
}

//------- structure related ---------//pure
float CNFalign_Basic::SingleOf2Pos_ACC_our_10_42_pure(int tPos,int sPos)
{
	int t_acc=tseq->ACC[tPos];
	int s_acc=seq->ACC[sPos];
	short diff = t_acc-s_acc>0?t_acc-s_acc:s_acc-t_acc;
	switch(diff)
	{
		case 0: return 3;
		case 1: return 1;
		case 2: return -2;
		default:
		{
			fprintf(stderr,"bad switch acc_10_42 here !!\n");
			exit(-1);
		}
	}
}
float CNFalign_Basic::SSOf2Pos2_pure(int tPos,int sPos)
{
	return (tseq->SS[tPos]==seq->SS[sPos]);
}
float CNFalign_Basic::SSOf2Pos3_pure(int tPos,int sPos)
{
	int pred_type = seq->SS[sPos];
	float score = 0;
	for(int i=0;i<3;i++)
		score += subm[pred_type][i]*seq->SS2_[sPos][i];	
	return score;
}
float CNFalign_Basic::SSOf2Pos4_pure(int tPos,int sPos)
{
	float score = 0;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			score += tseq->SS2[tPos][i] * seq->SS2[sPos][j] * subm[i][j];
	return score;
}


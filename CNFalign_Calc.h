#pragma once
#include <time.h>
#include <algorithm>
#include "Score.h"
#include "ScoreMatrix.h"
#include "mt19937ar.h"       //-> random number generator
#include "CNFalign_Feat.h"
using namespace std;


//========= extern states data ========//
//for Viterbi2 constraint
extern const int CB2CB; // = 1; //bond between two CB atoms
extern const int CA2CA; // = 2; //bond between two CA atoms
extern const int CB2CA; // = 4; //bond between CB and CA
extern const int CA2CB; // = 4; //bond between CA and CB
//other states transition
extern const int DUMMY_STATE;   // = -1;     //dummy state for dynamic programming
extern const int NO_ASSIGNMENT; // = 100;  //no_assignment state for dynamic programming

//====== class: CNFalign_Calc =====//
class CNFalign_Calc : public CNFalign_Feat
{
//in this class, X always correspond to the template and Y to the sequence
public:
	//---- constructor and destructor -----//
	CNFalign_Calc(SEQUENCE* t, SEQUENCE* s,float DISO_THRES=0.0);  //we set default to 0.0
	CNFalign_Calc(TEMPLATE* t, SEQUENCE* s,float DISO_THRES=0.0);  //we set default to 0.0
	CNFalign_Calc();
	~CNFalign_Calc();
	void CNFalign_Calc_basic(void); //-> init all macros and parameters

	//---- input parameter ---//__130331__//
	int TRACE_THRES;   // this parameter is for same-state-path record (default: 20)
	int BOUND_THRES;   // for bound
	int start_x,start_y,end_x,end_y; //for bound

	//---- penalty parameter ---//__140226__//
	//-> during alignment
	float MM_penalty; // -7   , -4.9 
	float XX_penalty; // -11  , -11
	float MI_penalty; // -3   , -2.6
	float II_penalty; // -3.7 , -2.6
	float MD_penalty; // -3   , -2.6
	float ID_penalty; // 0    , 0
	float DD_penalty; // -3.7 , -2.6
	//-> after alignment
	float MATCH_FIN_PENALTY;
	float GAP_FIN_PENALTY;
	//-> match bound
	float M_THRES;           // disabled
	float MATCH_BOUND_UPPER; // 15 , 
	float MATCH_BOUND_LOWER; // -5 ,
	//---- normalize factor ----//__140226__//
	float VITERBI2_FACTOR;  //4.285 , 1.0


	//======== related to Viterbi2 alignment ========//
	int SSMatch(int pos_x);
	Score ComputeViterbi(vector<pair<int, int> > & alignment);
	Score ComputeViterbi2(vector<pair<int, int> >& alignment_path,vector<double> & align_score_vector,int &lali_out,int NORM_or_PURE);

	//======== related to MEA alignment ========//
	//forward(s, i, j) is the sum of the probabilities of all the alignment paths ending at position (i, j) with state s
	//the occurring probability of state s at position (i, j) is included in forward(s, i, j)
	//it is stored in the log form to reduce numerical error
	ScoreMatrix* forward; //forward probability matrix
	//backward(s, i, j) is the sum of the probabilities of all the alignment paths starting at position (i, j) with state s
	//the occurring probability of state s at position (i, j) is not included in backward(s, i, j)
	//stored in the log form to reduce numerical error	
	ScoreMatrix* backward; //backward probability matrix
	//posterior matrix
	ScoreMatrix* Post; 
	//----- MEA alignment related ------//
	// delete fb matrices
	void ClearFB();
	//calculate the forward probability matrix, stored in the log form
	void ComputeForward();
	//calculate the backward probability matrix, stored in the log form
	void ComputeBackward();
	//calculate log of the partition function
	LogScore ComputePartitionCoefficient();
	//print posterior matrix
	void PrintPosterior(string &file);
	//get bound from alignment
	void Bound_Detect_Given_Alignment(vector<pair<int,int> > &alignment);
	void Alignment_Get_Bound_New(int moln1,int moln2,vector<pair<int,int> > &align,vector<pair<int,int> > &bound,int neib);
	//---- MEA alignment ---//
	Score MEA_Alignment(vector<pair<int,int> > & alignment_path);
	Score MEA_Alignment2(vector<pair<int, int> > & alignment_path);

	//======== Sample alignment ========//
	//set seeds
	void SetSeed48();
	mt19937 rand_gen;
	//sample an alignment
	Score SampleAnAlignment(vector<pair<int,int> >& alignment);

	//======== alignment fix functions ====//
	int AA26_to_AA20(int amino);
	void Get_Alignment(vector<pair<int, int> > &in,int *wali1,int *wali2,int moln1,int moln2);
	int Ali_To_Cor(int *AFP_Cor,int thres,int *ali1,int *ali2,int moln1,int moln2);
	void Ali_To_Vect(vector<pair<int, int> > &alignment,int *ali1,int *ali2,int moln1,int moln2);
	void Vect_To_Str(vector<pair<int, int> > &alignment, string &str1,string &str2,string &out1,string &out2);
	void Fix_Minor_Error(vector<pair<int, int> > &in,string &str1,string &str2);
	
};

//--alignment related --//
extern int *outer_ali1;                 //lenx+leny
extern int *outer_ali2;                 //lenx+leny
extern int *outer_AFP;                  //4*(lenx+leny)
extern float* MM_LogProb;

//--- feature create function ----//
extern void CNFalign_Calc_Total_Volumn_New(int length_x_,int length_y_);  //create all the needed data structure


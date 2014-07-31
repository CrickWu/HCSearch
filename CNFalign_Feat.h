#pragma once
#include "Score.h"
#include "CNFalign_Basic.h"
using namespace std;


//----- neuron network state ----//__140226__//
extern const int D_; // = 398;
extern const int D1; // = 398;
extern const int D2; // = 176;
extern const int D3; // = 176;
extern const int MM_;// = 8;
extern const int G_; // = 12;
extern const int _G; // = 7;
extern const int WD; // = MM_*G + MM_*G*D;//134720
extern const int MM_G; // = MM_ * G;
extern const int G_D;  // = G * D;

//=========== state and transition ============//
#define STATE_INSERT_X  "INSERT_X"
#define STATE_INSERT_Y  "INSERT_Y"
#define STATE_MATCH     "MATCH"
#define tM2M 0
#define tM2IX 1
#define tM2IY 2
#define tIX2M 3
#define tIX2IX 4
#define tIX2IY 5
#define tIY2M 6
#define tIY2IY 7
enum StateType {
  StateType_NONE,
  StateType_MATCH,
  StateType_INSERT_X,
  StateType_INSERT_Y
};
struct State {
  string name;
  int type;
  int mtype,stype;
  int emit_x;
  int emit_y;
  string param;
  LogScore val;
  LogScore counts;
  State (){}
  State (const State &rhs) :
    name(rhs.name), type(rhs.type), emit_x(rhs.emit_x), emit_y(rhs.emit_y),
    param(rhs.param), val(rhs.val), counts(rhs.counts){}
};
struct Transition {
  int from;
  int to;
  int ttype;
  string param;
  LogScore val;
  LogScore counts;
  Transition(){}
  Transition(const Transition &rhs) :
    from(rhs.from), to(rhs.to), param(rhs.param),
    val(rhs.val), counts(rhs.counts){}
};


//====== class: CNFalign_Feat =====//
class CNFalign_Feat : public CNFalign_Basic
{
//in this class, X always correspond to the template and Y to the sequence
public:
	//---- constructor and destructor -----//
	CNFalign_Feat(SEQUENCE* t, SEQUENCE* s,float DISO_THRES=0.0);  //we set default to 0.0
	CNFalign_Feat(TEMPLATE* t, SEQUENCE* s,float DISO_THRES=0.0);  //we set default to 0.0
	CNFalign_Feat();
	~CNFalign_Feat();
	void CNFalign_Feat_basic(void);  //-> init all macros and parameters
	
	//=========== weight related =============//
	//---- weight data structure ----//
	float* weights;                  // weights to calculate CNF
	vector<State> states;            // the states describing an alignment
	vector<Transition> transitions;  // the state transition matrix 
	//---- set functions ----//
	void SetWeights(float* x);       // set weights
	void SetStatesAndTransitions(    // set states and transitions
		vector<State>& s, vector<Transition>& trans); 

	//========== for training purpose ===========//
	//------ inner data -----//
	int Inner_Data_Creat;            // default=0, no Data_Creat (just link outer data)
	short *inner_0_feature_data2;    // lenx*leny*177
	float *inner_0_feature_data3;    // lenx*leny*177
	short *inner_0_non_zero_length;  // lenx*leny
	short *inner_1_feature_data2;    // lenx*177
	float *inner_1_feature_data3;    // lenx*177
	short *inner_1_non_zero_length;  // lenx
	short *inner_2_feature_data2;    // leny*177
	float *inner_2_feature_data3;    // leny*177
	short *inner_2_non_zero_length;  // leny
	void Set_Inner_Data(int METHOD=0);  //-> METHOD: 0, link outer data; 1, create inner data

	/////===============   main calculate part =========////
	//======== generate features ========//
	void GenerateFeatures();
	//-> temp_seq
	void GenerateFeatures_norm();
	void GenerateOneSample_norm(int currState, int leftState,int pos_x, int pos_y, 
		float *mutation_matrix, float *stru_matrix, float* features);
	void GenerateIniSample_Match_norm(int pos_x,int pos_y, float *mutation_features, float *stru_features);
	short GenerateIniSample_IX_norm(int currState, int leftState, int pos_x, short* data2, float* data3);
	short GenerateIniSample_IY_norm(int currState, int leftState, int pos_x, short* data2, float* data3);
	//-> seq_seq
	void GenerateFeatures_pure();
	void GenerateOneSample_pure(int currState, int leftState,int pos_x, int pos_y, 
		float *mutation_matrix, float *stru_matrix, float* features);
	void GenerateIniSample_Match_pure(int pos_x,int pos_y, float *mutation_features, float *stru_features);
	short GenerateIniSample_IX_pure(int currState, int leftState, int pos_x, short* data2, float* data3);
	short GenerateIniSample_IY_pure(int currState, int leftState, int pos_x, short* data2, float* data3);
	//======= calculate gap and match prob ====//
	void ComputeGap_Gate_Output();
	void Compute_Match_Protential();


	/////===============   vice calculate part =========////
	//-> basic function
	Score gate(Score x);
	float powFastLookup(const float val,const float ilog2,const unsigned int* pTable,const unsigned int  precision);
	//======== calculate log probability =======//
	int index_trans(int s, int t);
	Score getGateOutput(int x, int y, int m, int g);                                  //-> general purpose
	LogScore CalculateLogProb(int currState, int pos_x, int pos_y, int leftState);    //-> general purpose
	Score getGateOutput_0(int x, int y, int m, int g,int idx);                        //-> for match
	LogScore CalculateLogProb_0(int currState, int pos_x, int pos_y, int leftState);  //-> for match
	Score getGateOutput_1(int x, int y, int m, int g,int idx);                        //-> for IX
	LogScore CalculateLogProb_1(int currState, int pos_x, int pos_y, int leftState);  //-> for IX
	Score getGateOutput_2(int x, int y, int m, int g,int idx);                        //-> for IY
	LogScore CalculateLogProb_2(int currState, int pos_x, int pos_y, int leftState);  //-> for IY

};


//======== extern temporary data =======//
extern float feat_vector[480];
extern float weight_vector[480];
extern float gate_vector[120];
extern float gate_weight_vector[120];

//======= extern feature data ==========//
extern const int outer_feat_dim; //180
extern float outer_features[480];
//--main matrix --//
extern float *outer_mutation_matrix;    //lenx*leny*8
extern float *outer_stru_matrix;        //lenx*leny*26
extern short *outer_0_feature_data2;    //lenx*leny*177
extern float *outer_0_feature_data3;    //lenx*leny*177
extern short *outer_0_non_zero_length;  //lenx*leny
extern short *outer_1_feature_data2;    //lenx*177
extern float *outer_1_feature_data3;    //lenx*177
extern short *outer_1_non_zero_length;  //lenx
extern short *outer_2_feature_data2;    //leny*177
extern float *outer_2_feature_data3;    //leny*177
extern short *outer_2_non_zero_length;  //leny
//--vice matrix--//
extern float *MM_LogProb;               //lenx*leny*8
extern float *MI_LogProb;               //lenx+leny
extern float *II_LogProb;               //lenx+leny
extern float *MD_LogProb;               //lenx+leny
extern float *DD_LogProb;               //lenx+leny
extern float *ID_LogProb;               //lenx+leny

//--- feature create function ----//
extern void CNFalign_Feat_Total_Volumn_New(int length_x_,int length_y_);  //create all the needed data structure


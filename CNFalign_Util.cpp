#include "CNFalign_Util.h"
#include "CNF_constants.h"

namespace CNFalign_Util {
	void initialize_states_and_transitions(vector<State>& states,
			vector<Transition>& transitions) {

		//currently the states and transitions define a global alignment model
		//to do global-local alignment model, needs to define two more states at least.

		//define states
		states.clear();
		transitions.clear();
		State alignState1;
		alignState1.name = STATE_MATCH;
		alignState1.type = StateType_MATCH;
		alignState1.emit_x = 1;
		alignState1.emit_y = 1;

		State insertXState;
		insertXState.name = STATE_INSERT_X;
		insertXState.type = StateType_INSERT_X;
		insertXState.mtype = 3;
		insertXState.stype = 1;
		insertXState.emit_x = 1;
		insertXState.emit_y = 0;

		State insertYState;
		insertYState.name = STATE_INSERT_Y;
		insertYState.type = StateType_INSERT_Y;
		insertYState.mtype = 0;
		insertYState.stype = 2;
		insertYState.emit_x = 0;
		insertYState.emit_y = 1;

		states.push_back(alignState1);
		states.push_back(insertXState);
		states.push_back(insertYState);

		//define transitions
		for (int i = 0; i < (int)states.size(); i++)
			for (int j = 0; j < (int)states.size(); j++) {
				if (states[i].type == StateType_INSERT_Y && states[j].type
						== StateType_INSERT_X)
					continue;
				Transition trans;

				if (states[i].type == StateType_MATCH) {
					if (states[j].type == StateType_MATCH)
						trans.ttype = tM2M;
					if (states[j].type == StateType_INSERT_X)
						trans.ttype = tM2IX;
					if (states[j].type == StateType_INSERT_Y)
						trans.ttype = tM2IY;

				}
				if (states[i].type == StateType_INSERT_X) {
					if (states[j].type == StateType_MATCH)
						trans.ttype = tIX2M;
					if (states[j].type == StateType_INSERT_X)
						trans.ttype = tIX2IX;
					if (states[j].type == StateType_INSERT_Y)
						trans.ttype = tIX2IY;

				}
				if (states[i].type == StateType_INSERT_Y) {
					if (states[j].type == StateType_MATCH)
						trans.ttype = tIY2M;
					if (states[j].type == StateType_INSERT_Y)
						trans.ttype = tIY2IY;
				}

				trans.from = i;
				trans.to = j;
				transitions.push_back(trans);
			}
	}

	void generate_matrix(TEMPLATE* t,SEQUENCE* s,
			vector<State> &states, vector<Transition> &transitions,
			vector<pair<int, int> > &alignment) {
		//--- init ---
		CNFalign_Calc * align = new CNFalign_Calc(t, s, CNF_constants::DISO_THRES);
		align->SetStatesAndTransitions(states, transitions);
		align->GenerateFeatures();
		align->Set_Inner_Data();
		align->SetWeights(CNF_constants::CNFpred_Model_rel);
		//--- alignment ---
		align->ComputeGap_Gate_Output();
		align->Compute_Match_Protential();
		//align->TRACE_THRES=CNF_constants::TRACE_THRES;
		align->BOUND_THRES=CNF_constants::BOUND_THRES;
		//align->ComputeViterbi(alignment);
		int lali;
		vector<double> align_score_vector;
		double score = align->ComputeViterbi2(alignment, align_score_vector, lali, 1);

		/*
		string result_root = ".";
		//--- output ---//
		if(result_root.length()!=0){
			string hhpred_style_output="";
			//if(OUT_LEVEL>0)HHpred_Style_Output(alignmentResult,align_score_vector,t,s,hhpred_style_output,80);
			HHpred_Style_Output(alignmentResult,align_score_vector,t,s,hhpred_style_output,80);
			SaveToFile(result_root,t->temp_name,s->seq_name,t->sequence,s->sequence,alignmentResult,hhpred_style_output,score); //output SEQRES in TPL
			//else SaveToFile(result_root,t->temp_name,s->seq_name,t->dssp_sequence,s->sequence,alignmentResult,hhpred_style_output,score);       //output DSSP in TPL
		}
		*/
		/*
		 *if( post_root.length()!=0 ){
		 *    //-> get bound //
		 *    align->Bound_Detect_Given_Alignment(alignmentResult);  //__130408__//
		 *    //-> output post //
		 *    string post_file = post_root + "/" + t->temp_name+"-"+s->seq_name+".post";
		 *    align->ComputeForward();
		 *    align->ComputeBackward();
		 *    align->PrintPosterior(post_file);
		 *}
		 */

		// here MM_LogProb is modified by ``extern'' in CNFalign_Calc.h
		// and MM_LogProb is a ``global'' variable
		// MM_LogProb is contained in , thus cannot be destructed inside the function??? not really
		delete align;
	}
}

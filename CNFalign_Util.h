#pragma once
#include "template.h"
#include "seq.h"
#include "CNFalign_Calc.h"

namespace CNFalign_Util {
	/** initilize the states and transitions needed **/
	void initialize_states_and_transitions(vector<State>& states,
			vector<Transition>& transitions);
	/** generete the 8 scoring matrices **/
	void generate_matrix(TEMPLATE* t,SEQUENCE* s,
			vector<State> &states, vector<Transition> &transitions, vector<pair<int, int> > &alignment);
		
}

#pragma once
#include <iostream>
#include <vector>
#include "template.h"
#include "seq.h"
#include "CNFalign_Calc.h"
#include "CNFalign_Util.h"

typedef double (*lossfunc) (vector<pair<int, int> >);

class HCSearch {
	public:
		// this function utilizes the visible MM_LogProb[][] matrix as the scoring matrix
		HCSearch(int length_x, int length_y, int out_search_range = 4, int out_max_depth = 5) : st_map(3, vector<int>(3, -1)), search_range(out_search_range), max_depth(out_max_depth) {
			// length_x is the length of the template
			this->length_x = length_x;
			// length_y is the length of the target/sequence
			this->length_y = length_y;
		}
		void put_single_sample(lossfunc lf);
		void train_hc();
		void set_states_and_transitions(vector<State>& s, vector<Transition>& trans);
		// here, we implement the binary-modification as modifying the middle residue of the total alignment
		void divide_alignment(int start2, int end2, vector<pair<int, int> > &prev_alignment, int depth);
		// recompute the alignment template[start1, end2) and target[start2, end2)
		void compute_alignment(int start1, int start2, int end1, int end2, vector<pair<int, int> > &alignment);
		int loss_func(vector<pair<int, int> > &prev_alignment);
		// this function actually can be similar as with put_single_sample, since it also has to exploit the MM_LogProb[][] matrix at last
		vector<pair<int, int> > predict(SEQUENCE *s, TEMPLATE *t);
		vector<vector<int> > st_map;
		vector<pair<int, int> > correct_alignment, cur_alignment;

		// perform binary search one sequence and the search is for ``index'' snum, and it seems that the orginal non-pair appointed dp results can be different from any of the sub-searches (so we need to compare with the original one ?)
		int search_seq_pair(int snum, vector<pair<int, int> > &alignment);
		int common_pairs(vector<pair<int, int> > &align1, vector<pair<int, int> > &align2);
	private:
		// denotes the maximu number of depth the hcsearch can do
		int max_depth;
		// the number denotes the +- range will consider
		int search_range;
		// remaining .. (first_half) + (one pair) + (second_half) + remaining ..
		// include start - 1] [end = -[start, end)
		void concatenate(int start, int end, vector<pair<int, int> > &prev_alignment, vector<pair<int, int> > &first_half, pair<int, int> p, vector<pair<int, int> > &second_half, vector<pair<int, int> > &return_alignment);
		void put_constrains();
		int length_x, length_y;
		double h(vector<pair<int, int> > align);
		double c(vector<pair<int, int> > align);
		vector<State> states;
		vector<Transition> transitions;

		// constrains_*[i] represents the i-th training data
		// constrains_*[i][j] represents the j-th feature of the training data
		// constrains_*[i][n] represents the loss score (which implies the ranking)
		// in the training, only pairwise constrains_*[i]'s relations are generated
		vector<vector<double> > constrains_h, constrains_c;

		// if the rank function is of pairwise form, we could use 1(>) -1(<) 0(=) to represents the corresponding pairwise relations
		vector<vector<int> > rank_h, rank_c;
};

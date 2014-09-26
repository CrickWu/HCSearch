#pragma once
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "template.h"
#include "seq.h"
#include "CNFalign_Calc.h"
#include "CNFalign_Util.h"
#include "uGDT_Calc.h"

// sofia-ml includes
#include "sofia-ml/sofia-ml-methods.h"
#include "floatfann.h"

typedef double (*lossfunc) (vector<pair<int, int> >);
enum task_type {TRAIN = 1, PREDICT = 2};

class HCSearch {
	public:
		// this function utilizes the visible MM_LogProb[][] matrix as the scoring matrix
		HCSearch(TEMPLATE* t, SEQUENCE* s, int out_search_range, int out_seq_search_range, int out_max_depth, string &outfile) : st_map(3, vector<int>(3, -1)), search_range(out_search_range), seq_search_range(out_seq_search_range), max_depth(out_max_depth), qid_counter(1), out_train(outfile.c_str()), task(TRAIN), verbose(1), w(NULL) {
			this->t = t;
			this->s = s;
			// length_x is the length of the template
			this->length_x = t->length;
			// length_y is the length of the target/sequence
			this->length_y = s->length;
		}
		~HCSearch() {out_train.close();}
		void set_correct_alignment(const vector<pair<int, int > > &alignment);
		void put_single_sample(lossfunc lf);
		void train_hc();
		void set_states_and_transitions(vector<State>& s, vector<Transition>& trans);
		double get_percent(const vector<pair<int, int> > &alignment);
		void seq_divide_alignment(int start2, int end2, vector<pair<int, int> >&prev_alignment, int dpeth);
		// predict using the binary division method
		void predict_divide(int start2, int end2, vector<pair<int, int> > &prev_alignment, int depth);
		// here, we implement the binary-modification as modifying the middle residue of the total alignment
		void divide_alignment(int start2, int end2, vector<pair<int, int> > &prev_alignment, int depth);
		// recompute the alignment template[start1, end2) and target[start2, end2)
		void compute_alignment(int start1, int start2, int end1, int end2, vector<pair<int, int> > &alignment);
		double loss_func(vector<pair<int, int> > &prev_alignment);
		double loaded_loss_func(vector<double> &features);
		// this function actually can be similar as with put_single_sample, since it also has to exploit the MM_LogProb[][] matrix at last
		vector<pair<int, int> > predict(SEQUENCE *s, TEMPLATE *t);
		vector<vector<int> > st_map;
		vector<pair<int, int> > correct_alignment, cur_alignment;

		// perform binary search one sequence and the search is for ``index'' snum, and it seems that the orginal non-pair appointed dp results can be different from any of the sub-searches (so we need to compare with the original one ?)
		int search_seq_pair(int snum, vector<pair<int, int> > &alignment);
		int common_pairs(const vector<pair<int, int> > &align1, const vector<pair<int, int> > &align2);
		void output_map_inf();
		void output_layer_training(double priority, const vector<double> features, int qid, int depth, ostream &out);
		void load_model(const string& file_name);
		void set_seqs_template(TEMPLATE * s_t) {this->s_t = s_t;}
		// need to be moved to private
		// denotes the maximu number of depth the hcsearch can do
		int max_depth;
		// the number denotes the +- range will consider
		int search_range;
		int seq_search_range;
		ofstream out_predict;
		// reduced number of features
		bool reduced;
	private:
		SfWeightVector* w;
        fann* net;
		ofstream out_train, out_test;
		task_type task;
		// 0 - nothing, 1 - intermediate candidate, 2 - layer info
		int verbose;
		int non_zero_pairs_in_structure;
		int qid_counter;
		void generate_features(vector<pair<int, int> >&alignment, vector<double> &features);
		// remaining .. (first_half) + (one pair) + (second_half) + remaining ..
		// include start - 1] [end = -[start, end)
		void concatenate(int start, int end, vector<pair<int, int> > &prev_alignment, vector<pair<int, int> > &first_half, pair<int, int> p, vector<pair<int, int> > &second_half, vector<pair<int, int> > &return_alignment);
		void put_constrains();
		SEQUENCE* s;
		TEMPLATE* t;
		// the real structure for sequence
		TEMPLATE* s_t;
		int length_x, length_y;
		double h(vector<pair<int, int> > align);
		double c(vector<pair<int, int> > align);
		map<int, int, greater<int> > common_pair_count;
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

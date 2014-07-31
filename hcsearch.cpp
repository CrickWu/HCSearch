#include "hcsearch.h"
using namespace std;

void printPair(vector<pair<int, int> > &alignment) {
	for (int i = 0; i < alignment.size(); i++)
		cout << "[" << alignment[i].first << ", " << alignment[i].second << "] ";
	cout << endl;
}

void HCSearch::train_hc() {}
vector<pair<int, int> > HCSearch::predict(SEQUENCE *s, TEMPLATE *t) {}
void HCSearch::set_states_and_transitions(vector<State>& s, vector<Transition>& trans)
{
	states = s;
	transitions = trans;
	// generate st_map based on given transitions
	for (int m = 0; m < transitions.size(); m++) {
		int s = transitions[m].from;
		int t = transitions[m].to;
		st_map[s][t] = m;
	}
}

void HCSearch::put_constrains() {}
double HCSearch::h(vector<pair<int, int> >) {}
double HCSearch::c(vector<pair<int, int> >) {}
void HCSearch::put_single_sample(lossfunc lf) {
}

void HCSearch::compute_alignment(int start1, int start2, int end1, int end2, vector<pair<int, int> > &alignment) {
	if (start1 > end1 || start2 > end2) {
		cerr << "Unexpected bound values found in \"compute_alignment\"." << endl;
		exit(-1);
	}
	alignment.clear();

	//best stores the log of the alignment score probability
	ScoreMatrix best((int) states.size(), length_x + 1, length_y + 1);
	best.Fill(0);

	//traceback matrix is used for trace back
	ScoreMatrix traceback((int) states.size(), length_x + 1, length_y + 1);
	traceback.Fill(DUMMY_STATE);

	//Local alignment
	for (int i = 0; i <= length_x; i++)
		best(0, i, 0) = 0;
	for (int j = 0; j <= length_y; j++)
		best(0, 0, j) = 0;

	//calculate the whole alignment score matrix
	best(0, 0, 0) = 0;
	for (int i = start1; i <= end1; i++) {
		for (int j = start2; j <= end2; j++) {
			for (int m = 0; m < (int) transitions.size(); m++) {
				int s = transitions[m].from;
				int t = transitions[m].to;

				//Check the secondary structures
				//if (s == 0 && t != 0 && SSMatch(i - 1))
				//	continue;
				if (i - start1 >= states[t].emit_x && j - start2 >= states[t].emit_y) {
					Score new_score = best(s, i - states[t].emit_x, j
							- states[t].emit_y)// + CalculateLogProb(t, i, j, s);
					+ MM_LogProb[m*(length_x + 1)*(length_y+1) + i*(length_y+1) + j];
					if (best(t, i, j) < new_score) {
						best(t, i, j) = new_score;
						traceback(t, i, j) = s;
					}
				}
				if (best(t, i, j) < 0)
					best(t, i, j) = 0, traceback(t, i, j) = DUMMY_STATE;
			}
		}
	}

	// Do traceback //
	double maxv = 0;
	int ts = 0, lx = end1, ly = end2;
	for (int s = 0; s < (int) states.size(); s++)
		for (int x = start1; x <= end1; x++)
			for (int y = start2; y <= end2; y++) {
				if (maxv < best(s, x, y))
					maxv = best(s, x, y), ts = s, lx = x, ly = y;
			}
	Score score = best(ts, lx, ly);

	//string alignment_path;
	for (int i = end1; i != lx; i--)
		alignment.push_back(pair<int, int> (i, -ly));
	for (int j = end2; j != ly; j--)
		alignment.push_back(pair<int, int> (-end1, j));
	int ti = lx, tj = ly;
	while (ts != DUMMY_STATE && ti > start1 && tj > start2)
	{
		int new_ts = (int) traceback(ts, ti, tj);
		switch (states[ts].type)
		{
			case StateType_INSERT_X:
				alignment.push_back(pair<int, int> (ti, -tj));
				ti--;
				break;
			case StateType_INSERT_Y:
				alignment.push_back(pair<int, int> (-ti, tj));
				tj--;
				break;
			case StateType_MATCH:
				alignment.push_back(pair<int, int> (ti, tj));
				ti--;
				tj--;
				break;
			default:
				fprintf(stderr,"Unexpected value found in traceback matrix.");
				exit(-1);
		}
		ts = new_ts;
	}
	while (ti > start1)
		alignment.push_back(pair<int, int> (ti, -start2)), ti--;
	while (tj > start2)
		alignment.push_back(pair<int, int> (-start1, tj)), tj--;
	reverse(alignment.begin(), alignment.end());
}

int HCSearch::common_pairs(vector<pair<int, int> > &align1, vector<pair<int, int> > &align2) {
	int size1 = align1.size(), size2 = align2.size();
	int p1 = 0, p2 = 0;
	int res = 0;
	while (p1 < size1 && p2 < size2) {
		if (align1[p1].first <= 0 || align1[p1].second <= 0)
			p1++;
		else if (align2[p2].first <= 0 || align2[p2].second <= 0)
			p2++;
		else if (align1[p1].first > align2[p2].first || align1[p1].second > align2[p2].second)
			p2++;
		else if (align1[p1].first < align2[p2].first || align1[p1].second < align2[p2].second)
			p1++;
		// any necessity to check so ?
		else if (align1[p1] == align2[p2]) {
			res++;
			/*
			 *cout << "(" << align1[p1].first << "," << align1[p1].second << ") ";
			 */
			p1++;
			p2++;
		}
	}
	/*
	 *cout << endl;
	 *if (res == 3)
	 *    printPair(align1), printPair(align2);
	 */
	return res;
}

int HCSearch::loss_func(vector<pair<int, int> > &prev_alignment) {
	return common_pairs(prev_alignment, correct_alignment);
}

// range: [start2, end2)
// [start1, end1)
// start2 is the regular form of resiudes
void HCSearch::divide_alignment(int start2, int end2, vector<pair<int, int> > &prev_alignment, int depth) {
	cout << "start from :" << start2 << ":" << end2 << endl;
	if (start2 >= end2)
		return;
	// range: prev_alignment[start_pos, end_pos]
	int start_pos = search_seq_pair(start2 + 1, prev_alignment);
	int end_pos = search_seq_pair(end2, prev_alignment);

	int start1 = prev_alignment[start_pos].first;
	// value == 0, then must align from 0
	// value > 0, align from the normal starting one (-1)
	// value < 0, adjust the start1 to preceding reisidues until we find the first residue of the template which is not aligned with a gap
	if (start1 <= 0) {
		int new_start_pos = start_pos - 1;
		while (new_start_pos >= 0 && prev_alignment[new_start_pos].second <= 0 )
			new_start_pos--;
		if (prev_alignment[new_start_pos + 1].first > 0) {
			start1 = prev_alignment[new_start_pos + 1].first - 1;
			start_pos = new_start_pos + 1;
		}
		else
			start1 = -start1;
	} else
		start1 = start1 - 1;

	int end1 = prev_alignment[end_pos].first;
	// value > 0, align from the normal starting one ()
	// value < 0, align from the normal starting one (abs)
	if (end1 < 0)
		end1 = -end1;

	// consider the situation when there is no residue in template
	// just fill it with gaps
	if (end1 < start1) {
		vector<pair<int, int> > tmp;
		for (int i = start2; i < end2; i++)
			tmp.push_back(pair<int, int>(-start1, i + 1));
		return;
	}

	// will not have == since there is possibility that we can change it from the gap to an alignment
	if (end1 - start1 < 1)
		return;
	//cout << "pos:" << start_pos << " " << end_pos << endl;
	int mid2 = (start2 + end2) / 2;
	int mid2_pos = search_seq_pair(mid2 + 1, prev_alignment);
	int old_mid1 = prev_alignment[mid2_pos].first;
	// including this residue
	// value == 0, then must align from 0
	// value > 0, align from the normal starting one (-1)
	// value < 0, align from the normal starting one (abs)
	if (old_mid1 <= 0)
		old_mid1 = -old_mid1;
	else
		old_mid1 = old_mid1 - 1;
	// this is a rough range: consider [0, 1]<->[3,7]
	int range = ((end1 - start1 + 1) / 2) < search_range ? ((end1 - start1 + 1) / 2) : search_range;
	/*
	 *cout << "old_mid1 " << old_mid1 << " old_mid2 " << mid2 << endl;
	 *cout << "seq 1: " << start1 <<" " << end1;
	 *cout << "seq 2: " << start2 <<" " << end2 << endl;
	 */

	// record the next parameter to expand (greedy), note that we will use a DFS method
	// which means that the later one (reference parameter to prev_alignment)
	// is influenced by the first method
	// should not initialize it with 0 since the modification set might not include the pre_alignment
	vector<pair<int, int> > opt_alignment = prev_alignment;
	int next_start = 0, next_mid = 0, next_end = 0;
	int opt_common_pair = loss_func(opt_alignment);

	for (int i = -range; i <= range; i++) {
		vector<pair<int, int> > first_half, second_half;
		int mid1 = old_mid1 + i;
		if (mid1 < start1 || mid1 >= end1)
			continue;
		//cout << "cur_mid: " << mid1 << " " << mid2 << ":";

		vector<pair<int, int> > tmp;
		// consider the forced aligned situation
		compute_alignment(start1, start2, mid1, mid2, first_half);
		compute_alignment(mid1 + 1, mid2 + 1, end1, end2, second_half);
		concatenate(start_pos, end_pos + 1, prev_alignment, first_half, pair<int, int>(mid1 + 1, mid2 + 1), second_half, tmp);
		// can be recored in some priority array in the object
		// together with the range
		/*
		printPair(tmp), printPair(prev_alignment);
		cout << "common_pair: " <<common_pairs(tmp, prev_alignment) << endl;
		divide_alignment(start2, mid2, tmp);
		divide_alignment(mid2 + 1, end2, tmp);
		*/
		int tmp_loss = loss_func(tmp);
		if (opt_common_pair < tmp_loss) {
			opt_alignment = tmp;
			opt_common_pair = tmp_loss;
		}

		// then contain a gap
		compute_alignment(start1, start2, mid1, mid2, first_half);
		compute_alignment(mid1, mid2 + 1, end1, end2, second_half);
		concatenate(start_pos, end_pos + 1, prev_alignment, first_half, pair<int, int>(-mid1, mid2 + 1), second_half, tmp);

		tmp_loss = loss_func(tmp);
		if (opt_common_pair < tmp_loss) {
			opt_alignment = tmp;
			opt_common_pair = tmp_loss;
		}
	}

	// acutually we can later modify it to be a class member to speed up
	prev_alignment = opt_alignment;
	cout << opt_common_pair << endl;
	//cout << "seq 1: " << start1 <<" " << end1;
	//cout << " seq 2: " << start2 <<" " << end2 << " ";
	//printPair(opt_alignment);
	if (depth + 1 < max_depth) {
		divide_alignment(start2, mid2, prev_alignment, depth+1);
		divide_alignment(mid2 + 1, end2, prev_alignment, depth+1);
	}
}

// include start - 1] [end = -[start, end)
void HCSearch::concatenate(int start, int end, vector<pair<int, int> > &prev_alignment, vector<pair<int, int> > &first_half, pair<int, int> p, vector<pair<int, int> > &second_half, vector<pair<int, int> > &return_alignment) {
	return_alignment.clear();
	if (start > 0)
		// note the contruction function will not include the last element
		return_alignment = vector<pair<int, int> >(prev_alignment.begin(), prev_alignment.begin() + start);
	return_alignment.insert(return_alignment.end(), first_half.begin(), first_half.end());
	return_alignment.push_back(p);
	return_alignment.insert(return_alignment.end(), second_half.begin(), second_half.end());
	if (end < prev_alignment.size()) {
		vector<pair<int, int> > last_part(prev_alignment.begin() + end, prev_alignment.end());
		return_alignment.insert(return_alignment.end(), last_part.begin(), last_part.end());
	}
}

int HCSearch::search_seq_pair(int snum, vector<pair<int, int> > &alignment) {
	//length_y is the number we need
	int start = 0;
	int end = alignment.size();
	int res = (start + end) / 2;
	// +- 1???
	while (alignment[res].second != snum && start < end) {
		if (abs(alignment[res].second) < snum) {
			start = res + 1;
		} else {
			end = res - 1;
		}
		res = (start + end) / 2;
	}
	return res;
}

void write_to_file(string str, vector<pair<int, int> > &alignment, TEMPLATE *t, SEQUENCE *s) {
	ofstream out;
	if (str.length() > 0)
		out.open(str.c_str());
	else
		out.open("tmp.fasta");
	int p1 = 0, p2 = 0;
	out << ">" << t->temp_name << endl;
	for (int i = 0; i < alignment.size(); i++) {
		if (alignment[i].first > 0)
			out << t->sequence[ alignment[i].first - 1 ];
		else
			out << "-";
	}
	out << "\n>" << s->seq_name << endl;
	for (int i = 0; i < alignment.size(); i++) {
		if (alignment[i].second> 0)
			out << s->sequence[ alignment[i].second - 1 ];
		else
			out << "-";
	}
	out.close();
}

void convert_to_pair(string s, vector<pair<int, int> > &alignment) {
	alignment.clear();
	ifstream in(s.c_str());
	string tmp_s, str1, str2;
	if (! in.is_open()) {
		cerr << "File not exists\n";
		exit(-1);
	}
	in >> tmp_s >> str1 >> tmp_s >> str2;

	int p1 = 0, p2 = 0;
	int char1, char2;
	for (int i = 0; i < str1.length(); i++) {
        if (str1[i] != '-') {
            char1 = p1 + 1;
            p1 = p1 + 1;
		} else
            char1 = -p1;

        if (str2[i] != '-') {
            char2 = p2 + 1;
            p2 = p2 + 1;
		} else
            char2 = -p2;
		alignment.push_back(pair<int, int>(char1, char2));
	}
	in.close();
}

int ide(vector<pair<int, int> >&alig1, vector<pair<int, int> > &alig2) {
	int len = alig1.size();
	for (int i = 0; i < len; i++)
		if (alig1[i]!= alig2[i])
			return i;

	return -1;
}

int non_zero_pairs(vector<pair<int, int> > &alignment) {
	int res = 0;
	for (int i = 0; i < alignment.size(); i++)
		res += alignment[i].first > 0 && alignment[i].second > 0;
	return res;
}

int main(int argc, char **argv) {

	// ./** tpl tgt
	char tmp_char[30];
	strcat(tmp_char, argv[1]);
	strcat(tmp_char, "-");
	strcat(tmp_char, argv[2]);
	strcat(tmp_char, ".str.fasta");
	//string str = "1cjmA-1texA.str.fasta";
	string str(tmp_char);
	vector<pair<int, int> > alig, gen;
	cout << "start\n";
	convert_to_pair(str, alig);
	cout << "bn 1\n";

	cout << "Hello World!" << endl;
	string temp_name(argv[1]), seq_name(argv[2]);
	temp_name = "./required_tpl/" + temp_name;
	seq_name = "./required_tgt/" + seq_name;
	TEMPLATE *t = new TEMPLATE(temp_name, ".", 1);
	SEQUENCE *s = new SEQUENCE( seq_name, ".", 1);
	cout << "bn 2\n";

	vector<State> states;
	vector<Transition> transitions;

	/* Compute the dynamic matrix */
	CNFalign_Util::initialize_states_and_transitions(states, transitions);
	CNFalign_Feat_Total_Volumn_New(t->length+100,s->length+100);
	CNFalign_Calc_Total_Volumn_New(t->length+100,s->length+100);
	cout << "bn 3\n";

	int window_size = 4;
	if (argc > 3)
		window_size  = atoi(argv[3]);
	HCSearch hc = HCSearch(t->length, s->length, window_size);
	hc.set_states_and_transitions(states, transitions);
	hc.correct_alignment = alig;
	CNFalign_Util::generate_matrix(t, s, states, transitions, gen);
	// convert_to_pair("filename", hc.correct_alignment);

	//length is the # of residues while MM_LogProb[i][j] represents the score that [i-1][j-1] needs

	/*
	int real[7][2] = { {2,1}, {-2,2}, {-2,3}, {3,4}, {-3,5}, {4,6}, {5,7} };
	for (int i = 0; i < 7; i++)
		hc.correct_alignment.push_back(pair<int,int>(real[i][0], real[i][1]));

	int data[9][2] = {{2, 0}, {-2, 1}, {3, -1}, {-3, 2}, {-3, 3}, {-3, 4}, {4, 5}, {-4, 6}, {5, 7}};
	vector<pair<int, int> > test_ali;
	for (int i = 0; i < 9; i++)
		test_ali.push_back(pair<int, int>(data[i][0], data[i][1]));
	cout << hc.loss_func(test_ali) << endl;

	hc.compute_alignment(0, 10, 10, 20, test_ali);
	// the second parameter is whether the index of the last ``negative'' number (of the same value) or a ``positive'' number's index
	hc.divide_alignment(0, 7, test_ali, 0);
	printPair(test_ali);
	*/

	cout << "before: " << hc.common_pairs(alig, gen) << " " << non_zero_pairs(alig) << " " << non_zero_pairs(gen);
	cout << "before size: " << hc.common_pairs(alig, gen) << " " << alig.size() << " " << gen.size();
	write_to_file("bf.fasta", gen, t, s);
	hc.divide_alignment(0, s->sequence.length(), gen, 0);
	cout << "after: " << hc.common_pairs(alig, gen) << " " << non_zero_pairs(alig) << " " << non_zero_pairs(gen);
	cout << "after size: " << hc.common_pairs(alig, gen) << " " << alig.size() << " " << gen.size();
	write_to_file("af.fasta", gen, t, s);

	/*
	 *cout << hc.search_seq_pair(11, test_ali) << " " << hc.search_seq_pair(14, test_ali) << " " << hc.search_seq_pair(16, test_ali) << endl;
	 */

	/*
	 *for (int i = 0; i < t->length; i++) {
	 *    cout << MM_LogProb[s->length + i] << " ";
	 *    cout << II_LogProb[i] << " ";
	 *}
	 */

	/* Move residues */

	/* Put constraints in the set */

	/* After iteratively operating on the previous rounds, put it all in a training mechanism and generate heuristic function ``H'' and cost function ``C'' */
	return 0;
}

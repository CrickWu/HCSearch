#include "CNFalign_Calc.h"

//for Viterbi2 constraint
const int CB2CB = 1; //bond between two CB atoms
const int CA2CA = 2; //bond between two CA atoms
const int CB2CA	= 4; //bond between CB and CA
const int CA2CB = 4; //bond between CA and CB
//other states transition
const int DUMMY_STATE = -1;     //dummy state for dynamic programming
const int NO_ASSIGNMENT = 100;  //no_assignment state for dynamic programming


//--alignment related --//
int *outer_ali1;                 //lenx+leny
int *outer_ali2;                 //lenx+leny
int *outer_AFP;                  //4*(lenx+leny)
//----------- new array --------//
void CNFalign_Calc_Total_Volumn_New(int length_x_,int length_y_)
{
	//--- alignment related --//
	outer_ali1=new int[length_x_+length_y_];
	outer_ali2=new int[length_x_+length_y_];
	outer_AFP=new int[4*(length_x_+length_y_)];
}

//======================== constructor and destructor ========================//
//-> default value for DISO_THRES is 0.0
CNFalign_Calc::CNFalign_Calc(SEQUENCE* t, SEQUENCE* s,float DISO_THRES) //-> pure
{
	//--- penalty parameter --//__140226__//
	//-> during alignment
	MM_penalty = -4.9;
	XX_penalty = -3.0;
	MI_penalty = -2.6;
	II_penalty = -2.6;
	MD_penalty = -2.6;
	ID_penalty = -2.6;
	DD_penalty = -2.6;
	//-> after alignment
	MATCH_FIN_PENALTY=0;
	GAP_FIN_PENALTY=0;
	//-> match bound
	M_THRES=0.0;
	MATCH_BOUND_UPPER=100000;   // should be calculated from statistics
	MATCH_BOUND_LOWER=-100000;  // should be calculated from statistics
	//-> final factor
	VITERBI2_FACTOR=1.0;

	//--- member init ---//
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
	CNFalign_Feat_basic();
	CNFalign_Calc_basic();
}
CNFalign_Calc::CNFalign_Calc(TEMPLATE* t, SEQUENCE* s,float DISO_THRES) //-> norm
{
	//--- penalty parameter --//__140226__//
	//-> during alignment
	MM_penalty = -7;
	XX_penalty = -11;
	MI_penalty = -3;
	II_penalty = -3.7;
	MD_penalty = -3;
	ID_penalty = 0;
	DD_penalty = -3.7;
	//-> after alignment
	MATCH_FIN_PENALTY=0;
	GAP_FIN_PENALTY=0;
	//-> match bound
	M_THRES=0.0;
	MATCH_BOUND_UPPER=15;
	MATCH_BOUND_LOWER=-5;
	//-> final factor
	VITERBI2_FACTOR=4.285;

	//--- member init ---//
	temp=t;
	seq=s;
	temp_p=(PROFILE*)t;
	seq_p=(PROFILE*)s;
	tnam=t->temp_name;
	snam=s->seq_name;
	length_x = t->length;
	length_y = s->length;
	//--- basic set ---//
	CNFalign_Basic_Norm(DISO_THRES);
	CNFalign_Feat_basic();
	CNFalign_Calc_basic();
}
CNFalign_Calc::CNFalign_Calc()
{
}
CNFalign_Calc::~CNFalign_Calc()
{
	//MEA data
	if (forward) delete forward;
	if (backward) delete backward;
	if (Post) delete Post;
}

//================= basic setup ==================//
void CNFalign_Calc::CNFalign_Calc_basic(void)
{
	//--- macro parameters ---//__130331__//
	TRACE_THRES=30000;    //disabled  //-> for short-cut alignment
	BOUND_THRES=-1;       //disabled  //-> for bound calculation for posterior
	//--- data init ---//
	//for MEA
	forward = 0;
	backward = 0;
	Post = 0;
	//--- bound init ---//__130408__/
	start_x=start_y=end_x=end_y=-1;
}


//============= compute viterbi2 =============//
//---- pre calculate for viterbi2 -----//
int CNFalign_Calc::SSMatch(int pos_x)
{
	int thisType = 0, lastType = 0;
	if (pos_x >= 0)
		thisType = temp->SS[pos_x];
	if (pos_x > 0)
		lastType = temp->SS[pos_x - 1];
	if (thisType == lastType && lastType != 2)
		return 1;
	if (pos_x == 0)
		return 0;
	return 0;
}

//----------------- Compute Viterbi path, that is, find the best alignment -------------//
Score CNFalign_Calc::ComputeViterbi(vector<pair<int, int> > & alignment)
{
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
	for (int i = 0; i <= length_x; i++) {
		for (int j = 0; j <= length_y; j++) {
			for (int m = 0; m < (int) transitions.size(); m++) {
				int s = transitions[m].from;
				int t = transitions[m].to;

				//Check the secondary structures
				//if (s == 0 && t != 0 && SSMatch(i - 1))
				//	continue;
				if (i >= states[t].emit_x && j >= states[t].emit_y) {

					Score new_score = best(s, i - states[t].emit_x, j
							- states[t].emit_y) + CalculateLogProb(t, i, j, s);
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
	int ts = 0, lx = length_x, ly = length_y;
	for (int s = 0; s < (int) states.size(); s++)
		for (int x = 0; x <= length_x; x++)
			for (int y = 0; y <= length_y; y++) {
				if (maxv < best(s, x, y))
					maxv = best(s, x, y), ts = s, lx = x, ly = y;
			}
	Score score = best(ts, lx, ly);

	//string alignment_path;
	for (int i = length_x; i != lx; i--)
		alignment.push_back(pair<int, int> (i, -ly));
	for (int j = length_y; j != ly; j--)
		alignment.push_back(pair<int, int> (-length_x, j));
	int ti = lx, tj = ly;
	while (ts != DUMMY_STATE && ti > 0 && tj > 0) 
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
	while (ti > 0)
		alignment.push_back(pair<int, int> (ti, 0)), ti--;
	while (tj > 0)
		alignment.push_back(pair<int, int> (0, tj)), tj--;
	reverse(alignment.begin(), alignment.end());

	//return
	return score;
}

//--------- main function for viterbi2 ---------//
Score CNFalign_Calc::ComputeViterbi2(vector<pair<int, int> > & alignment,vector<double> & align_score_vector,int &lali_out,int NORM_or_PURE) 
{
	alignment.clear();
	//best stores the log of the alignment score probability
	ScoreMatrix best((int) states.size(), length_x + 1, length_y + 1);
	ScoreMatrix cons((int) states.size(), length_x + 1, length_y + 1);
	ScoreMatrix dist(1, length_x + 1, length_x + 1);
	best.Fill(0);
	dist.Fill(-1);
	cons.Fill(-1);
	if(NORM_or_PURE==1) //norm only
	{
		//-- constrains --//start
		for (int i = 1; i < length_x; i++)
		{
			for (int j = i + 1; j <= length_x; j++)
			{
				if (temp->DistOf2AAs(i - 1, j - 1, CA2CA) < 8 && temp->DistOf2AAs(i - 1, j - 1, CA2CA) > 0)
				{
					dist(0, i, j) = 1;
				}
			}
		}
		//-- constrains --//over
	}

	//-- alignment --//
	ScoreMatrix traceback((int) states.size(), length_x + 1, length_y + 1);
	ScoreMatrix tracereco((int) states.size(), length_x + 1, length_y + 1);
	traceback.Fill(DUMMY_STATE);
	tracereco.Fill(0);
	for (int i = 0; i <= length_x; i++) best(0, i, 0) = 0;
	for (int j = 0; j <= length_y; j++) best(0, 0, j) = 0;

	//calculate the whole alignment score matrix
	best(0, 0, 0) = 0;
	best(1, 0, 0) = 0;
	best(2, 0, 0) = 0;
	for (int i = 0; i <= length_x; i++) 
	{
		for (int j = 0; j <= length_y; j++) 
		{
			for (int m = 0; m < (int) transitions.size(); m++) 
			{
				int s = transitions[m].from;
				int t = transitions[m].to;
				if (t != 0) 
				{
					if (i >= states[t].emit_x && j >= states[t].emit_y) 
					{
						Score new_score =  best(s, i - states[t].emit_x, j - states[t].emit_y);
						if(t==1)  // IX
						{
							if(s==0)
							{
								new_score += MI_LogProb[i-1] + MI_penalty;
							}
							else if(s==1)
							{
								//--- check for previous same-state-path
								int trace_reco = (int)tracereco(s, i - states[t].emit_x, j - states[t].emit_y);
								if(trace_reco < TRACE_THRES)  //TRACE_THRES: default is disabled
								{
									new_score += II_LogProb[i-1] + II_penalty;
								}
								else
								{
									double cursco = II_LogProb[i-1] + II_penalty;
									if(cursco>0.5)new_score+=cursco;
								}
							}
						}
						else if(t==2) // IY
						{
							if(s==0)
							{
								new_score += MD_LogProb[j-1] + MD_penalty;
							}
							else if(s==2)
							{
								//--- check for previous same-state-path
								int trace_reco = (int)tracereco(s, i - states[t].emit_x, j - states[t].emit_y);
								if(trace_reco < TRACE_THRES)  //TRACE_THRES: default is disabled
								{
									new_score += DD_LogProb[j-1] + DD_penalty;
								}
								else
								{
									double cursco = DD_LogProb[j-1] + DD_penalty;
									if(cursco>0.5)new_score+=cursco;
								}
							}
							else
							{
								new_score += ID_LogProb[j-1] + ID_penalty;
							}
						}
						//record best score
						if (best(t, i, j) < new_score) 
						{
							best(t, i, j) = new_score;
							traceback(t, i, j) = s;
							//-- record constrains --//
							if(NORM_or_PURE==1) //norm only
							{
								if (t == 2)cons(t, i, j) = -1;
								else cons(t, i, j) = cons(s, i - 1, j);
							}
							//-- record trace_record --//
							if (s != 0) tracereco(t, i, j) = tracereco(s, i - states[t].emit_x, j - states[t].emit_y)+1;
							else tracereco(t, i, j) = 0;
						}
					}
				} 
				else 
				{
					if (i >= states[t].emit_x && j >= states[t].emit_y )
					{
						//calc cur_bool
						int cur_bool=1;
						if(NORM_or_PURE==1) //norm only
						{
							int wscons = (int) cons(s, i - 1, j - 1);
							int valid=0;
							if(wscons!=-1)
							{
								if(dist(0, wscons, i)> 0.1) valid=1;
							}
							cur_bool = ( (s == 0 && SSMatch(i - 1) ) || wscons  == -1 || valid==1 );
						}

						//judge
						if ( cur_bool ) 
						{
							//check bound
							double aaa = MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j];
							if(aaa>=MATCH_BOUND_UPPER)aaa = MATCH_BOUND_UPPER;
							else if(aaa<=MATCH_BOUND_LOWER)aaa = MATCH_BOUND_LOWER;

							//---- solve head tail tail problem ----//start
							if(i==1 || j==1 || i==length_x || j==length_y) aaa = aaa - MM_penalty;
							//---- solve head tail tail problem ----//over

							//---- solve bad and missing problem ----//start
							Score new_score;
							if(NORM_or_PURE==1) //norm only
							{
								if( temp->residue[i-1]==20 || seq->residue[j-1]==20 || temp->isMissing[i-1]==1 )
								{
									new_score = best(s, i - states[t].emit_x, j - states[t].emit_y) + aaa + XX_penalty;
//									new_score = best(s, i - states[t].emit_x, j - states[t].emit_y) + aaa + MM_penalty;
								}
								else
								{
									new_score = best(s, i - states[t].emit_x, j - states[t].emit_y) + aaa + MM_penalty;
								}
							}
							else
							{
								new_score = best(s, i - states[t].emit_x, j - states[t].emit_y) + aaa + MM_penalty;
							}	
							//---- solve bad and missing problem ----//over

							//record best
							if (best(t, i, j) < new_score) 
							{
								best(t, i, j) = new_score;
								traceback(t, i, j) = s;
								if(NORM_or_PURE==1) //norm only
								{
									cons(t, i, j) = i;
								}
							}
						}
					}
				}
				if (best(t, i, j) < 0) 
				{
					best(t, i, j) = 0, traceback(t, i, j) = DUMMY_STATE;
					if(NORM_or_PURE==1) //norm only
					{
						cons(t, i, j) = -1;
					}
				}
			}
		}
	}

	//--- get max ---//
	Score maxv = 0;
	int ts = 0, lx = length_x, ly = length_y;
	for (int s = (int) states.size() - 1; s >= 0; s--)
		for (int x = 0; x <= length_x; x++)
			for (int y = 0; y <= length_y; y++) 
			{
				if (maxv < best(s, x, y))
					maxv = best(s, x, y), ts = s, lx = x, ly = y;
			}
	Score score_ret = best(ts, lx, ly);

	//--- get alignment ---//
	for (int i = length_x; i != lx; i--)
		alignment.push_back(pair<int, int> (i, -ly));
	for (int j = length_y; j != ly; j--)
		alignment.push_back(pair<int, int> (-length_x, j));
	int ti = lx, tj = ly;
	//--- final extract ---//
	while (ts != DUMMY_STATE && best(ts, ti, tj) > 0 && ti > 0 && tj > 0) 
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
	while (ti > 0)
		alignment.push_back(pair<int, int> (ti, 0)), ti--;
	while (tj > 0)
		alignment.push_back(pair<int, int> (0, tj)), tj--;

	reverse(alignment.begin(), alignment.end());

//----- alignment_fix ----//__130408__//
{
	if(NORM_or_PURE==1) //norm
	{
		Fix_Minor_Error(alignment,temp->sequence,seq->sequence);
	}
	else                //pure
	{
		Fix_Minor_Error(alignment,tseq->sequence,seq->sequence);
	}
}
//----- alignment_fix ----//__130408__//over

int lali;
double ranking_match_sco_ful;
double ranking_gap_score_ful;
double ws_new_score;
//----  deal with bad and missing -----//__130230__//
{
	int k;
	int *ali1=outer_ali1;
	int *ali2=outer_ali2;
	for(k=0;k<length_x;k++)ali1[k]=-1;
	for(k=0;k<length_y;k++)ali2[k]=-1;
	//get ali
	lali=0;
	int size=(int)alignment.size();
	for(k=0;k<size;k++)
	{
		int ii=alignment[k].first;
		int jj=alignment[k].second;
		if(ii>0 && jj>0)
		{
			lali++;
			//check bad (we only consider disorder > 0.5 in tgt as bad)
			if(NORM_or_PURE==1) //norm
			{
//				if(temp->residue[ii-1]==20 || seq->residue[jj-1]==20 || temp->isMissing[ii-1]==1 || seq->DISO[jj-1]>0.5 )
				if(temp->residue[ii-1]==20 || seq->residue[jj-1]==20 || temp->isMissing[ii-1]==1 )
				{
				}
				else
				{
					ali1[ii-1]=jj-1;
					ali2[jj-1]=ii-1;
				}
			}
			else                //pure
			{
//				if(tseq->residue[ii-1]==20 || seq->residue[jj-1]==20 || tseq->DISO[jj-1]>0.5 || seq->DISO[jj-1]>0.5 )
				if(tseq->residue[ii-1]==20 || seq->residue[jj-1]==20 )
				{
				}
				else
				{
					ali1[ii-1]=jj-1;
					ali2[jj-1]=ii-1;
				}
			}
		}
	}
	//get vector
	{
		//ali1->alignment
		//init
		alignment.clear();
		//start
		int i;
		int j;
		int wlen;
		int ii,jj;
		int pre_ii=0;
		int pre_jj=0;
		for(i=1;i<=length_x;i++)
		{
			ii=i;
			jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
			if(jj==-1)
			{
				continue;
			}
			else
			{
				jj++;
				//previous_path
				wlen=ii-pre_ii;
				for(j=1;j<wlen;j++)
				{
					pre_ii++;
					alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
				}
				wlen=jj-pre_jj;
				for(j=1;j<wlen;j++)
				{
					pre_jj++;
					alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
				}
				//current_path
				alignment.push_back (pair<int,int>(ii, jj)); //Match
				//update
				pre_ii=ii;
				pre_jj=jj;
			}
		}
		//termi
		pre_ii++;
		for(i=pre_ii;i<=length_x;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
		pre_jj++;
		for(i=pre_jj;i<=length_y;i++)alignment.push_back (pair<int,int>(-length_x, i));  //Iy
	}
	
	//----- recalculate score -----//
	lali=0;  //-> targ length
	ranking_match_sco_ful = 0;
	ranking_gap_score_ful = 0;
	ws_new_score = 0;
	align_score_vector.clear();
	{
		int i;
		int start,end;
		int size=(int)alignment.size();
		//get start and end
		start=-1;
		for(i=0;i<size;i++)
		{
			int ii=alignment[i].first;
			int jj=alignment[i].second;
			if(ii>0 && jj>0)
			{
				start=i;
				break;
			}
		}
		end=-1;
		for(i=size-1;i>=0;i--)
		{
			int ii=alignment[i].first;
			int jj=alignment[i].second;
			if(ii>0 && jj>0)
			{
				end=i;
				break;
			}
		}
		//valid check
		if(start!=-1 && end!=-1)
		{
			//re-calculate score
			int pre,cur;
			pre=0;
			int gap_len=0;
			double ranking_match_thres = M_THRES;
			double ranking_match_sco = MATCH_FIN_PENALTY;
			double ranking_gap_score = GAP_FIN_PENALTY;
			for(i=start;i<=end;i++)
			{
				cur=-1;
				int ii=alignment[i].first;
				int jj=alignment[i].second;
				if(ii>0 && jj>0) //match
				{
					cur=0;
					int ws_m = index_trans(pre, cur);
					//check bound
					double aaa=MM_LogProb[ws_m*(length_x+1)*(length_y+1) + ii*(length_y+1)+jj];
					if(aaa>=MATCH_BOUND_UPPER)aaa = MATCH_BOUND_UPPER;
					else if(aaa<=MATCH_BOUND_LOWER)aaa = MATCH_BOUND_LOWER;
					align_score_vector.push_back(aaa);
					aaa=aaa+MM_penalty;
					//--- modify for ranking ---//
					if(aaa>0 && aaa<ranking_match_thres)aaa=0;
					ws_new_score+=(aaa);
					ranking_match_sco_ful+=ranking_match_sco;
					//--- modify for ranking ---//over
					gap_len=0;
					//-- lali ---//
					lali++;
					//-- lali ---//over
				}
				else
				{
					if(ii>0)   //Ix
					{
						cur=1;
						if(pre==0)
						{
							ws_new_score += MI_LogProb[ii-1] + MI_penalty;
						}
						else if(pre==1)
						{
							if(gap_len<TRACE_THRES)
							{
								ws_new_score += II_LogProb[ii-1] + II_penalty;
							}
							else
							{
								double cursco = II_LogProb[ii-1] + II_penalty;
								if(cursco>0.5)ws_new_score += cursco;
							}
						}
						gap_len++;
						//-- final penalty for gap --//
						ranking_gap_score_ful+=(ranking_gap_score);
					}
					if(jj>0)   //Iy
					{
						//-- lali ---//
						lali++;
						//-- lali ---//over
						cur=2;
						if(pre==0)
						{
							ws_new_score += MD_LogProb[jj-1] + MD_penalty;
						}
						else if(pre==2)
						{
							if(gap_len<TRACE_THRES)
							{
								ws_new_score += DD_LogProb[jj-1] + DD_penalty;
							}
							else
							{
								double cursco = DD_LogProb[jj-1] + DD_penalty;
								if(cursco>0.5)ws_new_score += cursco;
							}
						}
						else
						{
							if(gap_len<TRACE_THRES)
							{
								ws_new_score += ID_LogProb[jj-1] + ID_penalty;
							}
							else
							{
								double cursco = ID_LogProb[jj-1] + ID_penalty;
								if(cursco>0.5)ws_new_score += cursco;
							}
						}
						gap_len++;
						//-- final penalty for gap --//
						ranking_gap_score_ful+=(ranking_gap_score);
					}
				}
				//next stage
				if(cur!=-1)pre=cur;
			}
		}
	}
}

	//---- final return ----//
	ws_new_score=ws_new_score+ranking_match_sco_ful+ranking_gap_score_ful;
	score_ret=score_ret+ranking_match_sco_ful+ranking_gap_score_ful;
	if(ws_new_score<0)ws_new_score=0;
	if(score_ret<0)score_ret=0;
	lali_out=lali;
	if(NORM_or_PURE==1) //norm
	{
		return ws_new_score/VITERBI2_FACTOR ;    //default: 4.285 or 1.0
	}
	else                //pure
	{
		return score_ret/VITERBI2_FACTOR ;    //default: 4.285 or 1.0
	}
}


//================ MEA alignment related functions ===============//

void CNFalign_Calc::ComputeForward()
{

	//initialize forward
	if (forward)
	{
		delete forward;
		forward = 0;
	}

//---- add bound ----//
if(start_x<0)start_x=0;
if(start_y<0)start_y=0;
if(end_x<0 || end_x >length_x)end_x=length_x;
if(end_y<0 || end_y >length_y)end_y=length_y;
//---- add bound ----//over



	// forward(s,i,j) = sum of the probs of all the paths from the beginning which end in state s and which end at grid point (i,j) //
	// the occurring probability of (s, i, j) is included in forward(s, i, j) //

	//the template has length_x residues, indexed by 1, 2,..., length_x in this matrix
	//the sequence has length_y residues, indexed by 1, 2,..., length_y in this matrix

	forward = new ScoreMatrix(states.size(), length_x + 1, length_y + 1);
	for (int i = 0; i <= length_x; i++)
		for (int j = 0; j <= length_y; j++)
			for (int s = 0; s < (int) states.size(); s++)
				(*forward)(s, i, j) = LogScore_ZERO;

	(*forward)(0, start_x, start_y) = 0; //log of 1

	//process
	for (int i = start_x; i <= end_x; i++)
	{
		for (int j = start_y; j <= end_y; j++)
		{
			for (int m = 0; m < (int) transitions.size(); m++)
			{
				int s = transitions[m].from;
				int t = transitions[m].to;
				if (i - states[t].emit_x >= start_x && j - states[t].emit_y >=start_y)
				{
					double new_score=0;
					if (t != 0)  //--- gap state ---//
					{
						if(t==1)
						{
							if(s==0)
							{
								if(i - states[t].emit_x >= start_x) new_score = MI_LogProb[i-1];
							}
							else if(s==1)
							{
								if(i - states[t].emit_x >= start_x) new_score = II_LogProb[i-1];
							}
						}
						else if(t==2)
						{
							if(s==0)
							{
								if(j - states[t].emit_y >= start_y) new_score = MD_LogProb[j-1];
							}
							else if(s==2)
							{
								if(j - states[t].emit_y >= start_y) new_score = DD_LogProb[j-1];
							}
							else
							{
								if(j - states[t].emit_y >= start_y) new_score = ID_LogProb[j-1];
							}
						}
					}
					else  //--- match state ---//
					{
						new_score = MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j];
					}
					//--- calculate score ---//
					LogScore_PLUS_EQUALS((*forward)(t, i, j), (*forward)(s, i
							- states[t].emit_x, j - states[t].emit_y)
							+ new_score);
							
				}//end of if(i >= states[t].emit_x && j >= states[t].emit_y)
			}
		}
	}
}
//this function computes the alignment probability matrix from back to front
void CNFalign_Calc::ComputeBackward()
{

	if (backward)
	{
		delete backward;
		backward = 0;
	}

//---- add bound ----//
if(start_x<0)start_x=0;
if(start_y<0)start_y=0;
if(end_x<0 || end_x >length_x)end_x=length_x;
if(end_y<0 || end_y >length_y)end_y=length_y;
//---- add bound ----//over

	// backward(s,i,j) = sum of the probabilities of all the alignment paths which begin at state s and grid point (i,j) //
	// however, the occurring probability of (s,i,j) is not included in backward(s,i,j) //

	//the template has length_x residues, indexed by 1, 2,..., length_x in this matrix
	//the sequence has length_y residues, indexed by 1, 2,..., length_y in this matrix

	backward = new ScoreMatrix(states.size(), length_x + 2, length_y + 2);
	for (int i = 0; i <= length_x + 1; i++)
		for (int j = 0; j <= length_y + 1; j++)
			for (int s = 0; s < (int) states.size(); s++)
				(*backward)(s, i, j) = LogScore_ZERO; //LogScore of 0

	for (int s = 0; s < (int) states.size(); s++) {
		(*backward)(s, end_x, end_y) = 0; //LogScore of 1
	}

	//process
	for (int i = end_x; i >= start_x; i--) 
	{
		for (int j = end_y; j >= start_y; j--) 
		{
			for (int m = 0; m < (int) transitions.size(); m++) 
			{
				int s = transitions[m].from;
				int t = transitions[m].to;
				if ( i + states[t].emit_x <= end_x && j + states[t].emit_y <= end_y) 
				{
					double new_score=0;
					if (t != 0)  //--- gap state ---//
					{
						if(t==1)
						{
							if(s==0)
							{
								if(i + states[t].emit_x <= end_x ) new_score = MI_LogProb[i + states[t].emit_x-1];
							}
							else if(s==1)
							{
								if(i + states[t].emit_x <= end_x ) new_score = II_LogProb[i + states[t].emit_x-1];
							}
						}
						else if(t==2)
						{
							if(s==0)
							{
								if(j + states[t].emit_y <= end_y ) new_score = MD_LogProb[j + states[t].emit_y-1];
							}
							else if(s==2)
							{
								if(j + states[t].emit_y <= end_y ) new_score = DD_LogProb[j + states[t].emit_y-1];
							}
							else
							{
								if(j + states[t].emit_y <= end_y ) new_score = ID_LogProb[j + states[t].emit_y-1];
							}
						}
					}
					else  //--- match state ---//
					{
						new_score = MM_LogProb[m*(length_x+1)*(length_y+1) + (i + states[t].emit_x)*(length_y+1)+(j + states[t].emit_y)];
					}
					//--- calculate score ---//
					LogScore_PLUS_EQUALS((*backward)(s, i, j), (*backward)(t, i
							+ states[t].emit_x, j + states[t].emit_y)
							+ new_score);

				}//end of if ((i + states[t].emit_x) <= length_x && (j + states[t].emit_y) <= length_y) 
			}
		}
	}
}



//----- clear forward, backward and Post -----//
void CNFalign_Calc::ClearFB()
{
	if (forward)delete forward;
	if (backward)delete backward;
	if (Post)delete Post;
	forward = 0;
	backward = 0;
	Post = 0;
}

//forward and backward matrices should be calculated before calling this function
//return log of the partition function
LogScore CNFalign_Calc::ComputePartitionCoefficient()
{
	if(forward==0 || backward==0)
	{
		fprintf(stderr,"Forward/backward matrix not yet computed.");
		exit(-1);	
	}


//---- add bound ----//
if(start_x<0)start_x=0;
if(start_y<0)start_y=0;
if(end_x<0 || end_x >length_x)end_x=length_x;
if(end_y<0 || end_y >length_y)end_y=length_y;
//---- add bound ----//over


	//the partition coefficient of the forward and backward matries may not exactly same
	//so we use sqrt(Z(forward)*Z(backward)) as the final partition coefficient if both matrices are calculated

	LogScore coeff_forward = LogScore_ZERO;
	LogScore coeff_backward = LogScore_ZERO;

	for (int s = 0; s < (int) states.size(); s++) {
		//forward probability
		if (forward)
			LogScore_PLUS_EQUALS(coeff_forward, (*forward)(s, end_x,end_y));
	}
	//backward probability
	LogScore_PLUS_EQUALS(coeff_backward, (*backward)(0, start_x, start_y));

	//return the average of forward and backward probability
	return (coeff_forward + coeff_backward) / 2;
}

void CNFalign_Calc::PrintPosterior(string &file)
{
//	ComputeForward();
//	ComputeBackward();
	LogScore partCoeff = ComputePartitionCoefficient();

	//output match
	ofstream fout(file.c_str());
	for (int i = 1; i <= length_y; i++)
		for (int j = 1; j <= length_x; j++) {
			if (temp->isMissing[j - 1])
				continue;
			double post = exp(((*forward)(0, j, i) + (*backward)(0, j, i)
					- partCoeff));
			fout << post << " ";
		}

	//output i<->j_g  //insert_y
	for (int i = 1; i <= length_y; i++)
		for (int j = 1; j <= length_x; j++) {
			if (temp->isMissing[j - 1])
				continue;
			double post = exp(((*forward)(2, j, i) + (*backward)(2, j, i)
					- partCoeff));
			fout << post << " ";
		}

	//output j<->i_g // insert_x
	for (int j = 1; j <= length_x; j++)
		for (int i = 1; i <= length_y; i++) {
			if (temp->isMissing[j - 1])
				continue;
			double post = exp(((*forward)(1, j, i) + (*backward)(1, j, i)
					- partCoeff));
			fout << post << " ";
		}

}


//============= bound detect given alignment ==========//-> for Posterior Matrix
void CNFalign_Calc::Bound_Detect_Given_Alignment(vector<pair<int,int> > &alignment)
{
	if(BOUND_THRES<0)return;
	int i;
	int moln1=length_x;
	int moln2=length_y;
	int *ali1=outer_ali1;
	int *ali2=outer_ali2;
	for(i=0;i<moln1;i++)ali1[i]=-1;
	for(i=0;i<moln2;i++)ali2[i]=-1;
	int size=(int)alignment.size();
	for(i=0;i<size;i++)
	{
		int ii=alignment[i].first;
		int jj=alignment[i].second;
		if(ii>0 && jj>0)
		{
			ali1[ii-1]=jj-1;
			ali2[jj-1]=ii-1;
		}
	}
	//--- detect bound ---//x
	for(i=0;i<moln1;i++)
	{
		if(ali1[i]!=-1)
		{
			start_x=i-BOUND_THRES;
			break;
		}
	}
	for(i=moln1-1;i>=0;i--)
	{
		if(ali1[i]!=-1)
		{
			end_x=i+1+BOUND_THRES;
			break;
		}
	}
	//--- detect bound ---//y
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]!=-1)
		{
			start_y=i-BOUND_THRES;
			break;
		}
	}
	for(i=moln2-1;i>=0;i--)
	{
		if(ali2[i]!=-1)
		{
			end_y=i+1+BOUND_THRES;
			break;
		}
	}
}

//======== WS_Written Alignment_Bound =========//-> for Posterior Matrix
//to use this result, just set following for j-loop
//for (int j = bound[i].first; j <= bound[i].second; j++)
void CNFalign_Calc::Alignment_Get_Bound_New(int moln1,int moln2,vector<pair<int,int> > &align,vector<pair<int,int> > &bound,int neib)
{
	int k;
	int ii,jj;
	int pre_ii,pre_jj;
	int size=(int)align.size();
	//[1]get real alignment
	bound.resize(moln1+1);
	bound[0].first=0;
	bound[0].second=0;
	pre_jj=0;
	for(k=0;k<size;k++)
	{
		ii=align[k].first;
		jj=align[k].second;
		if(ii>0&&jj>0)
		{
			bound[ii].first=jj;  //positive
			bound[ii].second=jj;
			pre_ii=ii+1;
			pre_jj=jj+1;
		}
		else
		{
			if(ii>0)
			{
				bound[ii].first=pre_jj;
				bound[ii].second=pre_jj;
			}
			else
			{
				bound[abs(ii)].second++;
			}
		}
	}

	//[2]get horizontal expand
	for(k=0;k<=moln1;k++)
	{
		ii=bound[k].first;
		jj=bound[k].second;
		//if(ii!=jj)continue;
		if(ii-neib<0)ii=0;
		else ii=ii-neib;
		if(jj+neib>moln2)jj=moln2;
		else jj=jj+neib;
		bound[k].first=ii;
		bound[k].second=jj;
	}
	for(int i=0;i<(int)bound.size();i++)
	{
		if(bound[i].first==0)
			bound[i].first=1;	
	}
}



//============== MEA_Alignment =============//
Score CNFalign_Calc::MEA_Alignment(vector<pair<int, int> > & alignment_path)
{
	alignment_path.clear();
//	ComputeForward();
//	ComputeBackward();
	LogScore partCoeff = ComputePartitionCoefficient();

	//initialize Post
	if (Post)
	{
		delete Post;
		Post = 0;
	}
	Post = new ScoreMatrix(3, length_x + 1, length_y + 1);

	//calculate Post
	for (int i = 0; i <= length_x; i++) {
		for (int j = 0; j <= length_y; j++) {
			for (int s = 0; s < (int) states.size(); s++)
				(*Post)(s, i, j) = exp(((*forward)(s, i, j) + (*backward)(
						s, i, j) - partCoeff));
		}
	}

	float T = 0.00;
	//best stores the log of the alignment score probability
	ScoreMatrix best((int) states.size(), length_x + 1, length_y + 1);

	best.Fill((Score) LogScore_ZERO);
	best(0, 0, 0) = 0;
	//traceback matrix is used for trace back
	ScoreMatrix traceback(3, length_x + 1, length_y + 1);
	traceback.Fill(DUMMY_STATE);

	//calculate the whole alignment score matrix
	for (int i = 0; i <= length_x; i++) {
		for (int j = 0; j <= length_y; j++) {
			for (int m = 0; m < (int) transitions.size(); m++) {
				int s = transitions[m].from;
				int t = transitions[m].to;
				if (i >= states[t].emit_x && j >= states[t].emit_y) {
					Score new_score = best(s, i - states[t].emit_x, j
							- states[t].emit_y) + (*Post)(t, i, j);
					float offset = 0;
					
					if (t == 0)
						offset -= T;
					else
						offset -= T / 2;
					
					new_score += offset;
					if (best(t, i, j) < new_score) {
						best(t, i, j) = new_score;
						traceback(t, i, j) = s;
					}
				}
			}
		}
	}

	Score maxv = 0;
	int ts = 0, lx = length_x, ly = length_y;
	for (int s = 0; s < (int) states.size(); s++)
		if (maxv < best(s, lx, ly))
			maxv = best(s, lx, ly), ts = s;
	ts = 0;
	Score score = best(ts, lx, ly);
	//string alignment_path;
	int ti = lx, tj = ly;
	while (ts != DUMMY_STATE && best(ts, ti, tj) > 0 && ti > 0 && tj > 0)
	{
		int new_ts = (int) traceback(ts, ti, tj);
		switch (states[ts].type) 
		{
			case StateType_INSERT_X:
				//ti--;
				alignment_path.push_back(pair<int, int> (ti, -tj));
				ti--;
				break;
			case StateType_INSERT_Y:
				//tj--;
				alignment_path.push_back(pair<int, int> (-ti, tj));
				tj--;
				break;
			case StateType_MATCH:
				//ti--; tj--;
				alignment_path.push_back(pair<int, int> (ti, tj));
				ti--;
				tj--;
				break;
			default:
				fprintf(stderr,"Unexpected value found in traceback matrix.");
				exit(-1);
		}
		ts = new_ts;
	}
	//ti++; tj++;
	while (ti > 0)
		alignment_path.push_back(pair<int, int> (ti, 0)), ti--;
	while (tj > 0)
		alignment_path.push_back(pair<int, int> (0, tj)), tj--;
	reverse(alignment_path.begin(), alignment_path.end());

	//return
	return score;
}

// With Distance Constraint
Score CNFalign_Calc::MEA_Alignment2 (vector<pair<int, int> > & alignment)
{
	alignment.clear();
//	ComputeForward();
//	ComputeBackward();
	LogScore partCoeff = ComputePartitionCoefficient();

	//initialize Post
	if (Post)
	{
		delete Post;
		Post = 0;
	}
	Post = new ScoreMatrix(3, length_x + 1, length_y + 1);

	//calculate Post
	for (int i = 0; i <= length_x; i++) {
		for (int j = 0; j <= length_y; j++) {
			(*Post)(0, i, j) = exp(((*forward)(0, i, j) + (*backward)(0,
					i, j) - partCoeff));
		}
	}

	//best stores the log of the alignment score probability
	ScoreMatrix best((int) states.size() + 1, length_x + 1, length_y + 1);
	ScoreMatrix dist(1, length_x + 1, length_x + 1);
	ScoreMatrix cons((int) states.size(), length_x + 1, length_y + 1);
	best.Fill((Score) LogScore_ZERO);

	cons.Fill(-1);
	for (int i = 1; i < length_x; i++)
		for (int j = i + 1; j <= length_x; j++)
			if (temp->DistOf2AAs(i - 1, j - 1, CA2CA) < 7 && temp->DistOf2AAs(i
					- 1, j - 1, CA2CA) > 0)
				dist(0, i, j) = 1;


	//traceback matrix is used for trace back
	ScoreMatrix traceback((int) states.size(), length_x + 1, length_y + 1);
	traceback.Fill(DUMMY_STATE);

	best(0, 0, 0) = 0;

	for (int i = 0; i <= length_x; i++) {
		for (int j = 0; j <= length_y; j++) {
		//for (int j =  temp_label_number[i].first; j <= temp_label_number[i].second; j++) {
			for (int m = 0; m < (int) transitions.size(); m++) {
				int s = transitions[m].from;
				int t = transitions[m].to;
				if (t != 0) {
					if (i >= states[t].emit_x && j >= states[t].emit_y) {
						Score new_score = best(s, i - states[t].emit_x, j
								- states[t].emit_y) + (*Post)(t, i, j);
						if (best(t, i, j) < new_score) {
							best(t, i, j) = new_score;
							traceback(t, i, j) = s;
							if (t == 2)
								cons(t, i, j) = -1;
							else
								cons(t, i, j) = cons(s, i - 1, j);
						}
					}
				} else {
					if (i >= states[t].emit_x && j >= states[t].emit_y && ((s
							== 0 && SSMatch(i - 1)) || cons(s, i - 1, j - 1)
							== -1 || dist(0, (int) cons(s, i - 1, j - 1), i)
							> 0.1)) {
						if ( (s == 0 && SSMatch(i - 1) ) || cons(s, i - 1, j - 1)
								== -1 || cons(s, i - 1, j - 1) >= 0) {

							Score new_score = best(s, i - states[t].emit_x, j
									- states[t].emit_y) + (*Post)(t, i, j);
							if (best(t, i, j) < new_score) {
								best(t, i, j) = new_score;
								traceback(t, i, j) = s;
								cons(t, i, j) = i;
							}
						}
					}
				}

				if (best(t, i, j) < 0) {
					best(t, i, j) = 0, traceback(t, i, j) = DUMMY_STATE;
					cons(t, i, j) = -1;
				}
			}
		}
	}

	// Do traceback //
	double maxv = 0;

	int ts = 0, lx = length_x, ly = length_y;
	for (int s = 0; s < (int) states.size(); s++)
		for (int x = 0; x <= length_x; x++)
			for (int y = 0; y <= length_y; y++) {
				if (maxv < best(s, x, y))
					maxv = best(s, x, y), ts = s, lx = x, ly = y;
			}

	for (int i = length_x; i != lx; i--)
		alignment.push_back(pair<int, int> (i, -ly));
	for (int j = length_y; j != ly; j--)
		alignment.push_back(pair<int, int> (-length_x, j));

	int ti = lx, tj = ly;
	while (ts != DUMMY_STATE && best(ts, ti, tj) > 0 && ti > 0 && tj > 0) {

		int new_ts = (int) traceback(ts, ti, tj);
		switch (states[ts].type) 
		{
			case StateType_INSERT_X:
				//ti--;
				alignment.push_back(pair<int, int> (ti, -tj));
				ti--;
				break;
			case StateType_INSERT_Y:
				//tj--;
				alignment.push_back(pair<int, int> (-ti, tj));
				tj--;
				break;
			case StateType_MATCH:
				//ti--; tj--;
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
	while (ti > 0)
		alignment.push_back(pair<int, int> (ti, 0)), ti--;
	while (tj > 0)
		alignment.push_back(pair<int, int> (0, tj)), tj--;
	reverse(alignment.begin(), alignment.end());

	//return
	return maxv;
}


//================= sample alignment part =============//
//-------- set seed ------//
void CNFalign_Calc::SetSeed48()
{
	srand(time(NULL));
	unsigned long randomSeed = rand();
	rand_gen.init_genrand(randomSeed);
}
//----------- sample alignment ----------//
//sample an alignment from the Forward probability matrix
//the Forward matrix should be calculated before calling this function
// With Geometric Constraints
Score CNFalign_Calc::SampleAnAlignment(vector<pair<int, int> >& alignment)
{
	//init alignment
	alignment.clear();

	//check forward matrix
	if(forward==0)
	{
		fprintf(stderr,"Forward matrix not yet computed.");
		exit(-1);
	}

	//we do sampling from back to front, the starting position is at length_x+1 and length_y+1
	int pos_x = length_x + 1;
	int pos_y = length_y + 1;
	int current_state = 0; //have to be a match state

	//given current state,sample its left state
	while (pos_x > 1 && pos_y > 1) {
		int left_pos_x = pos_x - states[current_state].emit_x;
		int left_pos_y = pos_y - states[current_state].emit_y;

		vector<LogScore> state_prob(states.size(), 0); //log probability of each possible state in the left position
		for (int s = 0; s < (int)states.size(); s++) {

			//state_prob[s] is equal to forward(s,..) * the probability of alignment of pos_x to pos_y

			if (pos_x == length_x + 1 && pos_y == length_y + 1)
				state_prob[s] = (*forward)(s, left_pos_x, left_pos_y);
			else
				state_prob[s] = (*forward)(s, left_pos_x, left_pos_y)
						+ CalculateLogProb(current_state, pos_x, pos_y, s);
		}

		vector<double> prefix_sum_of_prob(states.size() + 1, LogScore_ZERO); // the log of the prefix sum probability of the states in the left position

		//LogScore_ADD(x,y) = log (e^x + e^y)
		for (int s = 0; s < (int)states.size(); s++)
			prefix_sum_of_prob[s + 1] = LogScore_ADD(prefix_sum_of_prob[s],
					state_prob[s]);

		//prefix_sum_of_prob[states.size()] should be equal to the sum probability of all the states

		//normalize all the probabilities 
		for (int s = 0; s < (int)states.size(); s++)
			prefix_sum_of_prob[s + 1] = prefix_sum_of_prob[s + 1]
					- prefix_sum_of_prob[states.size()];

		//prefix_sum_of_prob[states.size()] should be equal to 1

//		double p = drand48(); //a random double between 0 and 1
		double p = rand_gen.genrand_real3() ;

		//randomly set a sampled state in case there are some numerical errors in probability calculation
//		int state_sampled = random() % states.size();
		int state_sampled = rand_gen.genrand_int32() % states.size();

		//LogScore_EXP (x) = e^x
		for (int s = 0; s < (int)states.size(); s++)
			if (p < LogScore_EXP(prefix_sum_of_prob[s + 1])) {
				//state s is sampled
				state_sampled = s;
				break;
			}

		//save the sampled state

		switch (states[current_state].type) 
		{
			case StateType_INSERT_X:
				alignment.push_back(pair<int, int> (left_pos_x, -left_pos_y));
				break;
			case StateType_INSERT_Y:
				alignment.push_back(pair<int, int> (-left_pos_x, left_pos_y));
				break;
			case StateType_MATCH:
				alignment.push_back(pair<int, int> (left_pos_x, left_pos_y));
				break;
			default:
				fprintf(stderr,"Unexpected value found in traceback matrix.");
				exit(-1);
		}

		//update current_state, pos_x, pos_y

		current_state = state_sampled;
		pos_x = left_pos_x;
		pos_y = left_pos_y;
	}

	//deal with head gap
	while (pos_x > 1) {
		alignment.push_back(pair<int, int> (pos_x - 1, 0));
		pos_x--;
	}

	while (pos_y > 1) {
		alignment.push_back(pair<int, int> (0, pos_y - 1));
		pos_y--;
	}

	reverse(alignment.begin(), alignment.end());

	return 0;
}


//========================= alignment fix ======================//__130408__//
//----- minor proc ----//
int CNFalign_Calc::AA26_to_AA20(int amino)
{
	switch (amino)
	{
		case 0:return 0;        //A
		case 1:return 20;       //B
		case 2:return 1;        //C
		case 3:return 2;        //D
		case 4:return 3;        //E
		case 5:return 4;        //F
		case 6:return 5;        //G
		case 7:return 6;        //H
		case 8:return 7;        //I
		case 9:return 20;       //J
		case 10:return 8;       //K
		case 11:return 9;       //L
		case 12:return 10;      //M
		case 13:return 11;      //N
		case 14:return 20;      //O
		case 15:return 12;      //P
		case 16:return 13;      //Q
		case 17:return 14;      //R
		case 18:return 15;      //S
		case 19:return 16;      //T
		case 20:return 20;      //U
		case 21:return 17;      //V
		case 22:return 18;      //W
		case 23:return 20;      //X
		case 24:return 19;      //Y
		case 25:return 20;      //Z
		default:return 20;
	}
}
void CNFalign_Calc::Get_Alignment(vector<pair<int, int> > &in,int *wali1,int *wali2,int moln1,int moln2)
{
	int i;
	int ii,jj;
	int size=(int)in.size();
	for(i=0;i<moln1;i++)wali1[i]=-1;
	for(i=0;i<moln2;i++)wali2[i]=-1;
	for(i=0;i<size;i++)
	{
		ii=in[i].first;
		jj=in[i].second;	
		if(ii>0 && jj>0)
		{
			wali1[ii-1]=jj-1;
			wali2[jj-1]=ii-1;
		}
	}
}
int CNFalign_Calc::Ali_To_Cor(int *AFP_Cor,int thres,int *ali1,int *ali2,int moln1,int moln2) 
{
	int i;
	int num;
	int ii,jj;
	int count;
	int isFirst;
	int isLast;
	int type;
	int head1,head2;
	int index;

	//init
	num=0;
	AFP_Cor[0]=0;
	head1=-1;
	head2=-1;
	isLast=0;
	isFirst=1;
	type=0;
	count=9999999;
	ii=-1;
	jj=-1;
	for(i=0;i<moln2;i++)
	{
		if(ali2[i]==-1) //purely blank
		{
			if(isFirst==0)
			{
				if(count>=thres) 
				{
					AFP_Cor[0]++;
					index=AFP_Cor[0];
					AFP_Cor[4*index+0] = type;
					AFP_Cor[4*index+1] = head1;
					AFP_Cor[4*index+2] = head2;
					AFP_Cor[4*index+3] = count;
					num+=count;
				}
				else
				{
				}
				count=0;
				isFirst=1;
			}
			continue;
		}

		if(isFirst==1)
		{
ws_init:
			isFirst=0;
			ii=ali2[i];
			type=1; // >0 mode
			jj=i;
			count=1;
			head1=ii;
			head2=jj;
			continue;
		}
		if(i==jj+1&&ali2[i]==ii+1)
		{
			ii=ali2[i];
			jj=i;
			count++;
			continue;
		}

ws_end:
		if(count>=thres) 
		{
			AFP_Cor[0]++;
			index=AFP_Cor[0];
			AFP_Cor[4*index+0] = type;
			AFP_Cor[4*index+1] = head1;
			AFP_Cor[4*index+2] = head2;
			AFP_Cor[4*index+3] = count;
			num+=count;
		}
		else
		{
		}

		if(isLast==1)goto end;
		else goto ws_init;
	}

	if(count==9999999)goto end;
	isLast=1;
	goto ws_end;
end:
	return num;
}
void CNFalign_Calc::Ali_To_Vect(vector<pair<int, int> > &alignment,int *ali1,int *ali2,int moln1,int moln2)
{
	//ali1->alignment
	//init
	alignment.clear();
	//start
	int i;
	int j;
	int wlen;
	int ii,jj;
	int pre_ii=0;
	int pre_jj=0;
	for(i=1;i<=moln1;i++)
	{
		ii=i;
		jj=ali1[i-1];  //ali1 starts from 0, correspondence also from 0
		if(jj==-1)
		{
			continue;
		}
		else
		{
			jj++;
			//previous_path
			wlen=ii-pre_ii;
			for(j=1;j<wlen;j++)
			{
				pre_ii++;
				alignment.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
			}
			wlen=jj-pre_jj;
			for(j=1;j<wlen;j++)
			{
				pre_jj++;
				alignment.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
			}
			//current_path
			alignment.push_back (pair<int,int>(ii, jj)); //Match
			//update
			pre_ii=ii;
			pre_jj=jj;
		}
	}
	//termi
	pre_ii++;
	for(i=pre_ii;i<=moln1;i++)alignment.push_back (pair<int,int>(i, -pre_jj)); //Ix
	pre_jj++;
	for(i=pre_jj;i<=moln2;i++)alignment.push_back (pair<int,int>(-moln1, i));  //Iy
}
void CNFalign_Calc::Vect_To_Str(vector<pair<int, int> > &alignment, string &str1,string &str2,string &out1,string &out2)
{
	out1.clear();
	out2.clear();
	int i;
	int ii,jj;
	int len=(int)alignment.size();
	for(i=0;i<len;i++)
	{
		ii=alignment[i].first;
		jj=alignment[i].second;
		if(ii>0 && jj>0) //match
		{
			out1.push_back(str1[ii-1]);
			out2.push_back(str2[jj-1]);
		}
		else
		{
			if(ii>0)  // IX
			{
				out1.push_back(str1[ii-1]);
				out2.push_back('-');
			}
			if(jj>0)  // IY
			{
				out1.push_back('-');
				out2.push_back(str2[jj-1]);
			}
		}
	}
}
//----- major proc ----//
void CNFalign_Calc::Fix_Minor_Error(vector<pair<int, int> > &in,string &str1,string &str2)
{
	//create
	int moln1=(int)str1.length();
	int moln2=(int)str2.length();
	int *ali1=outer_ali1;
	int *ali2=outer_ali2;
	int *AFP=outer_AFP;
	//init
	Get_Alignment(in,ali1,ali2,moln1,moln2);
	Ali_To_Cor(AFP,1,ali1,ali2,moln1,moln2);
	//process
	int i,j;
	int ii,jj,len;
	int size=AFP[0];
	for(i=1;i<=size;i++)
	{
		//start
		ii=AFP[4*i+1];
		jj=AFP[4*i+2];
		len=AFP[4*i+3];
		//head
		for(j=1;j<=moln1+moln2;j++)
		{
			int pos1=ii-j;
			int pos2=jj-j;
			if(pos1>=0 && pos2>=0)
			{
				if(ali1[pos1]==-1 && ali2[pos2]==-1)
				{
					int c1=AA26_to_AA20(str1[pos1]-'A');
					int c2=AA26_to_AA20(str2[pos2]-'A');
					if(c1!=20 && c2!=20)
					{
						if(c1==c2)
						{
							ali1[pos1]=pos2;
							ali2[pos2]=pos1;
						}
						else break;
					}
					else break;
				}
				else break;
			}
			else break;
		}
		//tail
		for(j=1;j<=moln1+moln2;j++)
		{
			int pos1=ii+len+j-1;
			int pos2=jj+len+j-1;
			if(pos1<moln1 && pos2<moln2)
			{
				if(ali1[pos1]==-1 && ali2[pos2]==-1)
				{
					int c1=AA26_to_AA20(str1[pos1]-'A');
					int c2=AA26_to_AA20(str2[pos2]-'A');
					if(c1!=20 && c2!=20)
					{
						if(c1==c2)
						{
							ali1[pos1]=pos2;
							ali2[pos2]=pos1;
						}
						else break;
					}
					else break;
				}
				else break;
			}
			else break;
		}
	}
	//final 
	Ali_To_Vect(in,ali1,ali2,moln1,moln2);
}


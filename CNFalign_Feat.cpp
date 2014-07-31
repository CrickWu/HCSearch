#include "CNFalign_Feat.h"

//----- neuron network state ----//__140226__//
const int D_ = 398;
const int D1 = 398;
const int D2 = 176;
const int D3 = 176;
const int MM_ = 8;
const int G_ = 12;
const int _G = 7;
const int WD = MM_*G_ + MM_*G_*D_; //134720
const int MM_G = MM_ * G_;
const int G_D = G_ * D_;


//------- temp structure -----//
float feat_vector[480];
float weight_vector[480];
float gate_vector[120];
float gate_weight_vector[120];

//============= WS_Total_Volumn ===========//__121030__//
const int outer_feat_dim=180;
float outer_features[480];
//--- temp matrix ---//
float *outer_mutation_matrix;    //lenx*leny*8
float *outer_stru_matrix;        //lenx*leny*26
//--- feature matrix ---//
short *outer_0_feature_data2;    //lenx*leny*177
float *outer_0_feature_data3;    //lenx*leny*177
short *outer_0_non_zero_length;  //lenx*leny
short *outer_1_feature_data2;    //lenx*177
float *outer_1_feature_data3;    //lenx*177
short *outer_1_non_zero_length;  //lenx
short *outer_2_feature_data2;    //leny*177
float *outer_2_feature_data3;    //leny*177
short *outer_2_non_zero_length;  //leny
//--vice matrix--//
float *MM_LogProb;               //lenx*leny*8
float *MI_LogProb;               //lenx+leny
float *II_LogProb;               //lenx+leny
float *MD_LogProb;               //lenx+leny
float *DD_LogProb;               //lenx+leny
float *ID_LogProb;               //lenx+leny

//----------- new array --------//
void CNFalign_Feat_Total_Volumn_New(int length_x_,int length_y_)
{
	length_x_++;
	length_y_++;
	//--- temp matrix ---//
	outer_mutation_matrix=new float[length_x_*length_y_*8];
	outer_stru_matrix=new float[length_x_*length_y_*26];
	//--- feature matrix ---//
	outer_0_feature_data2=new short[length_x_*length_y_*outer_feat_dim];
	outer_0_feature_data3=new float[length_x_*length_y_*outer_feat_dim];
	outer_0_non_zero_length=new short[length_x_*length_y_];
	outer_1_feature_data2=new short[length_x_*outer_feat_dim];
	outer_1_feature_data3=new float[length_x_*outer_feat_dim];
	outer_1_non_zero_length=new short[length_x_];
	outer_2_feature_data2=new short[length_y_*outer_feat_dim];
	outer_2_feature_data3=new float[length_y_*outer_feat_dim];
	outer_2_non_zero_length=new short[length_y_];
	//--- proc matrix ---//
	MM_LogProb=new float[length_x_*length_y_*8];
	MI_LogProb=new float[length_x_+length_y_];
	II_LogProb=new float[length_x_+length_y_];
	MD_LogProb=new float[length_x_+length_y_];
	DD_LogProb=new float[length_x_+length_y_];
	ID_LogProb=new float[length_x_+length_y_];
}



//======================== constructor and destructor ========================//
//-> default value for DISO_THRES is 0.0
CNFalign_Feat::CNFalign_Feat(SEQUENCE* t, SEQUENCE* s,float DISO_THRES) //-> pure
{
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
	CNFalign_Feat_basic();
	CNFalign_Basic_Pure(DISO_THRES);
}
CNFalign_Feat::CNFalign_Feat(TEMPLATE* t, SEQUENCE* s,float DISO_THRES) //-> norm
{
	//--- member init ---//
	temp=t;
	seq=s;
	temp_p=(PROFILE*)t;
	seq_p=(PROFILE*)s;
	length_x = t->length;
	length_y = s->length;
	tnam=t->temp_name;
	snam=s->seq_name;
	//--- basic set ---//
	CNFalign_Feat_basic();
	CNFalign_Basic_Norm(DISO_THRES);
}
CNFalign_Feat::CNFalign_Feat()
{
}
CNFalign_Feat::~CNFalign_Feat()
{
	//delete data
	if(Inner_Data_Creat==1)
	{
		if(inner_0_feature_data2)  delete [] inner_0_feature_data2;    //lenx*leny*177
		if(inner_0_feature_data3)  delete [] inner_0_feature_data3;    //lenx*leny*177
		if(inner_0_non_zero_length)delete [] inner_0_non_zero_length;  //lenx*leny
		if(inner_1_feature_data2)  delete [] inner_1_feature_data2;    //lenx*177
		if(inner_1_feature_data3)  delete [] inner_1_feature_data3;    //lenx*177
		if(inner_1_non_zero_length)delete [] inner_1_non_zero_length;  //lenx
		if(inner_2_feature_data2)  delete [] inner_2_feature_data2;    //leny*177
		if(inner_2_feature_data3)  delete [] inner_2_feature_data3;    //leny*177
		if(inner_2_non_zero_length)delete [] inner_2_non_zero_length;  //leny	             
	}
}

//================= basic setup ==================//
void CNFalign_Feat::CNFalign_Feat_basic(void)
{
	//--- macro parameters ---//__130331__//
	Inner_Data_Creat=0;   //no data_create
	//for gate
	inner_0_feature_data2=0;    //lenx*leny*177
	inner_0_feature_data3=0;    //lenx*leny*177
	inner_0_non_zero_length=0;  //lenx*leny
	inner_1_feature_data2=0;    //lenx*177
	inner_1_feature_data3=0;    //lenx*177
	inner_1_non_zero_length=0;  //lenx
	inner_2_feature_data2=0;    //leny*177
	inner_2_feature_data3=0;    //leny*177
	inner_2_non_zero_length=0;  //leny
}
//-------- ws_fast set_weight -------//
void CNFalign_Feat::SetWeights(float* x)
{
	weights = x;
}
//----- set states and traitions ----//
void CNFalign_Feat::SetStatesAndTransitions(vector<State>& s, vector<Transition>& trans)
{
	states = s;
	transitions = trans;
}
//----- set inner data --------//
void CNFalign_Feat::Set_Inner_Data(int METHOD)  //-> METHOD: 0, link outer data; 1, create inner data
{
	if(METHOD==0)
	{
		Inner_Data_Creat=0;
		inner_0_feature_data2    = outer_0_feature_data2;    //lenx*leny*177
		inner_0_feature_data3    = outer_0_feature_data3;    //lenx*leny*177
		inner_0_non_zero_length  = outer_0_non_zero_length;  //lenx*leny    
		inner_1_feature_data2    = outer_1_feature_data2;    //lenx*177     
		inner_1_feature_data3    = outer_1_feature_data3;    //lenx*177     
		inner_1_non_zero_length  = outer_1_non_zero_length;  //lenx         
		inner_2_feature_data2    = outer_2_feature_data2;    //leny*177     
		inner_2_feature_data3    = outer_2_feature_data3;    //leny*177     
		inner_2_non_zero_length  = outer_2_non_zero_length;  //leny
	}
	else
	{
		int k;
		int length_x_=length_x+1;
		int length_y_=length_y+1;
		//delete data
		if(Inner_Data_Creat==1)
		{
			if(inner_0_feature_data2)  delete [] inner_0_feature_data2;    //lenx*leny*177
			if(inner_0_feature_data3)  delete [] inner_0_feature_data3;    //lenx*leny*177
			if(inner_0_non_zero_length)delete [] inner_0_non_zero_length;  //lenx*leny
			if(inner_1_feature_data2)  delete [] inner_1_feature_data2;    //lenx*177
			if(inner_1_feature_data3)  delete [] inner_1_feature_data3;    //lenx*177
			if(inner_1_non_zero_length)delete [] inner_1_non_zero_length;  //lenx
			if(inner_2_feature_data2)  delete [] inner_2_feature_data2;    //leny*177
			if(inner_2_feature_data3)  delete [] inner_2_feature_data3;    //leny*177
			if(inner_2_non_zero_length)delete [] inner_2_non_zero_length;  //leny
		}
		//create
		Inner_Data_Creat=1;
		inner_0_feature_data2=new short[length_x_*length_y_*outer_feat_dim];
		inner_0_feature_data3=new float[length_x_*length_y_*outer_feat_dim];
		inner_0_non_zero_length=new short[length_x_*length_y_];
		inner_1_feature_data2=new short[length_x_*outer_feat_dim];
		inner_1_feature_data3=new float[length_x_*outer_feat_dim];
		inner_1_non_zero_length=new short[length_x_];
		inner_2_feature_data2=new short[length_y_*outer_feat_dim];
		inner_2_feature_data3=new float[length_y_*outer_feat_dim];
		inner_2_non_zero_length=new short[length_y_];
		//copy
		for(k=0;k<length_x_*length_y_*outer_feat_dim;k++)inner_0_feature_data2[k]=outer_0_feature_data2[k];
		for(k=0;k<length_x_*length_y_*outer_feat_dim;k++)inner_0_feature_data3[k]=outer_0_feature_data3[k];
		for(k=0;k<length_x_*length_y_;k++)inner_0_non_zero_length[k]=outer_0_non_zero_length[k];
		for(k=0;k<length_x_*outer_feat_dim;k++)inner_1_feature_data2[k]=outer_1_feature_data2[k];
		for(k=0;k<length_x_*outer_feat_dim;k++)inner_1_feature_data3[k]=outer_1_feature_data3[k];
		for(k=0;k<length_x_;k++)inner_1_non_zero_length[k]=outer_1_non_zero_length[k];
		for(k=0;k<length_y_*outer_feat_dim;k++)inner_2_feature_data2[k]=outer_2_feature_data2[k];
		for(k=0;k<length_y_*outer_feat_dim;k++)inner_2_feature_data3[k]=outer_2_feature_data3[k];
		for(k=0;k<length_y_;k++)inner_2_non_zero_length[k]=outer_2_non_zero_length[k];
	}
}

//================== generate feature ==============//
void CNFalign_Feat::GenerateFeatures()
{
	if(NORM_or_PURE==1) //temp<->seq (1), norm
	{
		GenerateFeatures_norm();
	}
	else                //seq<->seq (0), pure
	{
		GenerateFeatures_pure();
	}
}

//----------------- template-sequence features -----------------//
//compute the feature values at one alignment position (pos_x,pos_y), giving the state in the left position
//the vector<float> features store all the calculated values
void CNFalign_Feat::GenerateIniSample_Match_norm(int pos_x,int pos_y,float* mutation_features,float* stru_features)
{
	pos_x--;
	pos_y--;
	//---- mutation features ----//8
	float mScore;
	int index = 0;	
	mScore = MutationOf2Pos4(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore);
	//--mut4__--//
	float C = MutationOf2Pos4__(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore*C);
	//--mut4(3x)--//
	float s = MutationOf2Pos4(pos_x-1, pos_y-1);
	mScore += s;
	s = MutationOf2Pos4(pos_x+1, pos_y+1);
	mScore += s;
	mutation_features[index++] = (int)(1000*mScore);
	//--mut5--//
	mScore = MutationOf2Pos5(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore);
	//--mut4_--//
	mScore = MutationOf2Pos4_(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore);
	//--ss2--//
	float ss_score2 = SSOf2Pos2_norm(pos_x,pos_y);
	mutation_features[index++] = (int)(1000*ss_score2);
	//--ss3--//
	float ss_score3 = SSOf2Pos3_norm(pos_x,pos_y);
	mutation_features[index++] = (int)(1000*ss_score3);
	//--acc--//
	float acc_our_10_42 = SingleOf2Pos_ACC_our_10_42_norm(pos_x,pos_y);
	mutation_features[index++] = (int)(1000*acc_our_10_42);
	//---- structure features ----//26
	index = 0;		
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			if( temp->SS2_[pos_x][i] ==0 || seq->SS2_[pos_y][j] ==0)
				stru_features[index++] = (float)0;
			else
				stru_features[index++] = (int)(1000*temp->SS2_[pos_x][i]*seq->SS2_[pos_y][j]);	
		}
	}
	for(int i=0;i<8;i++){
		if(temp->SS8_[pos_x][i]==0 || seq->SS8_[pos_y][i] ==0)
			stru_features[index++] = (float)0;
		else
			stru_features[index++] = (int)(1000*temp->SS8_[pos_x][i]*seq->SS8_[pos_y][i]);
	}
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			if(temp->acc_[pos_x][i]==0 || seq->acc_our_10_42_[pos_y][j]==0)
				stru_features[index++] = (float)0;
			else
				stru_features[index++] = (int)(1000*temp->acc_[pos_x][i]*seq->acc_our_10_42_[pos_y][j]);
		}
	}
}
short CNFalign_Feat::GenerateIniSample_IX_norm(int currState, int leftState, int pos_x, short* data2,float* data3)
{
	pos_x--;
	int window_size = 5;
	float template_amino_freq = 0;
	short count = 1,index =1;
	data2[0]=0;data3[0]=1;


	//--- first part ---//
	for (int i = 0; i < 20; i++) 
	{
		//---- judge ----//
		if (pos_x < 0 || pos_x>length_x-1 )
			template_amino_freq = 1;
		else
			template_amino_freq = (int)(100*temp->EmissionProb_original[pos_x][i]);
		//--- record ---//
		if(template_amino_freq!=0){
			data2[index] = count; data3[index++] = template_amino_freq;
		}
		count++;
	}//1-20

	//--- second part ---//
	for(int w=-window_size;w<=window_size;w++){
		float ss_score=0;
		if (pos_x+w >= 0 && pos_x+w <= length_x-1) {
			for(int i=0;i<3;i++){
				ss_score = (int)(100*temp->SS2_[pos_x+w][i]);
				if(ss_score!=0){
					data2[index]=count; data3[index++]=ss_score;
				}
				count++;
			}
			for(int i=0;i<8;i++){
				ss_score = (int)(100*temp->SS8_[pos_x+w][i]);
				if(ss_score!=0){
					data2[index]=count;data3[index++]=ss_score;
				}
				count++;
			}
			for(int i=0;i<3;i++){
				ss_score = (int)(100*temp->acc_[pos_x+w][i]);
				if(ss_score!=0){
					data2[index]=count;data3[index++]=ss_score;
				}
				count++;
			}
		}
		else{
			count+=14;
		}
	}
	return index;
}
short CNFalign_Feat::GenerateIniSample_IY_norm(int currState, int leftState, int pos_y, short* data2, float* data3)
{
	pos_y--;
	int window_size = 5;
	float seq_amino_freq = 0;
	short count = 1,index = 1;
	data2[0]=0;data3[0]=1;	

	//--- first part ---//
	for (int i = 0; i < 20; i++) {
		//--- judge ----//
		if (pos_y < 0 || pos_y > length_y-1)
			seq_amino_freq = 1;
		else
			seq_amino_freq = (int)(100*seq->EmissionProb_original[pos_y][i]);
		//--- record ---//
		if(seq_amino_freq!=0){
			data2[index] = count; data3[index++] = seq_amino_freq;
		}
		count++;
	}//2-21

	//--- second part ---//
	for(int w=-window_size;w<=window_size;w++){
		float ss_score=0;
		if (pos_y+w >= 0 && pos_y+w <= length_y-1) {
			for (int i = 0; i < 3; i++){
				ss_score = (int)(100*seq->SS2_[pos_y+w][i]);
				if(ss_score!=0){	
					data2[index] = count; data3[index++] = ss_score;
				}
				count++;
			}
			for(int i=0;i<8;i++){
				ss_score = (int)(100*seq->SS8_[pos_y+w][i]);
				if(ss_score!=0){
					data2[index] = count; data3[index++] = ss_score;
				}
				count++;
			}
			for(int i=0;i<3;i++){
				ss_score = (int)(100*seq->acc_our_10_42_[pos_y+w][i]);
				if(ss_score!=0){
					data2[index] = count; data3[index++] = ss_score;
				}
				count++;
			}
		}
		else{
			count+=14;
		}
	}
	return index;
}
void CNFalign_Feat::GenerateOneSample_norm(int currState, int leftState,int pos_x, int pos_y, 
	float* mutation_matrix, float* stru_matrix, float* features) 
{
	int window_size = 5;
	pos_x--;
	pos_y--;

	//--- init ---//
	float seq_amino_freq = 0, template_amino_freq = 0;
	float sub = 0;
	int index = 0;
	for (int i = 0; i < 20; i++)
	{
		template_amino_freq = temp->EmissionProb_original[pos_x][i];
		seq_amino_freq = seq->EmissionProb_original[pos_y][i];
		sub = (int)(10000 * template_amino_freq * seq_amino_freq);
		features[index++] = sub;
	}//1-20

	//--- window size ---//
	for(int w = -window_size;w<=window_size;w++)
	{
		if( pos_x+w<0 || pos_y+w<0 || pos_x+w > length_x-1 || pos_y+w > length_y -1 )
		{
			float mScore=0;
			features[index++] = 0;
			features[index++] = 0;
			float s = MutationOf2Pos4(pos_x+w-1, pos_y+w-1);
			mScore += s;
			s = MutationOf2Pos4(pos_x+w+1, pos_y+w+1);
			mScore += s;
			features[index++] =(int)(1000*mScore);
			for(int i=0;i<5;i++)features[index++]=0;
		}
		else
		{
			for(int i=0;i<8;i++)features[index++]= *(mutation_matrix+((pos_x+w)*(length_y+1)+ (pos_y+w))*8+i);
		}
	}
	int wspair=AApair(pos_x, pos_y);
	if(wspair<0)
	{
		features[index++]=0;
		features[index++]=0;
	}
	else
	{
		features[index++] = (int)(1000*CC50[wspair]);
		features[index++] = (int)(1000*HDSM[wspair]);
	}
				
	//2. add secondary structure score. We can use a vector to represent this score
	for(int w = -window_size;w<=window_size;w++)
	{	
		if( pos_x+w<0 || pos_y+w<0 || pos_x+w > length_x-1 || pos_y+w > length_y -1 )
		{
			for(int i=0;i<26;i++)features[index++]=0;
		}
		else
		{
			for(int i=0;i<26;i++)features[index++]= *(stru_matrix+((pos_x+w)*(length_y+1)+ (pos_y+w))*26+i);
		}
	}
}
void CNFalign_Feat::GenerateFeatures_norm() 
{
	//match
	for(int i=0; i<=length_x;i++)
	{
		int cur_index=(i-1)*(length_y+1);
		for(int j=0; j<= length_y; j++)
		{
			int t=0;
			if (i - states[t].emit_x >= 0 && j - states[t].emit_y >= 0) 
			{
				GenerateIniSample_Match_norm(i,j,outer_mutation_matrix+(cur_index + (j-1))*8,outer_stru_matrix+(cur_index + (j-1))*26 );
			}
		}
	}
	//IX
	for (int i = 0; i <= length_x; i++)
	{
		if (i - states[1].emit_x >= 0)
		{
			short l = GenerateIniSample_IX_norm(1,0,i,outer_1_feature_data2+i*outer_feat_dim,outer_1_feature_data3+i*outer_feat_dim);
			outer_1_non_zero_length[i]=l;
		}
	}
	//IY
	for (int j = 0; j <= length_y; j++)
	{
		if(j - states[2].emit_y >= 0)
		{
			short l = GenerateIniSample_IY_norm(2,0,j,outer_2_feature_data2+j*outer_feat_dim,outer_2_feature_data3+j*outer_feat_dim);
			outer_2_non_zero_length[j]=l;
		}
	}
	//generate_one_sample
	for (int i = 0; i <= length_x; i++)
	{
		for (int j = 0; j <= length_y; j++)
		{
			int t = 0;
			int s = 0;
			if (i - states[t].emit_x >= 0 && j - states[t].emit_y >= 0)
			{
				if(t==0)
				{	
					for(int d=0;d<396;d++)outer_features[d]=0;
					GenerateOneSample_norm(t,s,i,j,outer_mutation_matrix,outer_stru_matrix,outer_features);
					int cur_index=(i*(length_y+1)+j)*outer_feat_dim;
					outer_0_feature_data2[cur_index]=0;
					outer_0_feature_data3[cur_index]=1;
					short count = 1;
					for (short d = 0; d< 396; d++) 
					{
						if(outer_features[d]!= 0 && count<outer_feat_dim)
						{
							outer_0_feature_data2[cur_index+count]=d+1;
							outer_0_feature_data3[cur_index+count]=outer_features[d];
							count++;
						}
					}
					outer_0_non_zero_length[i*(length_y+1)+j]=count;
				}
			}
		}
	}
}

//---------------- sequence-sequence features ------------------//
void CNFalign_Feat::GenerateIniSample_Match_pure(int pos_x,int pos_y,float* mutation_features,float* stru_features)
{
	pos_x--;
	pos_y--;
	//---- mutation features ----//8
	float mScore;
	int index = 0;
	//--mut4--//
	mScore = MutationOf2Pos4(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore);
	//--mut4_--//
	mScore = MutationOf2Pos4_(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore);
	//--mut5--//
	mScore = MutationOf2Pos5(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore);
	//--mut6--//
	mScore = MutationOf2Pos6(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore);
	//--mut7--//
	mScore = MutationOf2Pos7(pos_x, pos_y);
	mutation_features[index++] = (int)(1000*mScore);
	//--ss2--//
//	float ss_score2 = SSOf2Pos2_pure(pos_x,pos_y);
	float ss_score2 = 0;
	mutation_features[index++] = (int)(1000*ss_score2);
	//--blosum--//
	float blosum62 = 0.1*(BLOSUM62_Calc(pos_x,pos_y));
	mutation_features[index++] = (short)(1000*blosum62);
	//--ss4--//
//	float ss_score4 = SSOf2Pos4_pure(pos_x,pos_y);
	float ss_score4 = 0;
	mutation_features[index++] = (short)(1000*ss_score4);

	//-- other structure --//
	index = 0;
	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			stru_features[index++] = (int)0;
//			if( tseq->SS2[pos_x][i] ==0 || seq->SS2[pos_y][j] ==0)
//				stru_features[index++] = (int)0;
//			else
//				stru_features[index++] = (int)(1000*tseq->SS2[pos_x][i]*seq->SS2[pos_y][j]);
		}
	}
	for(int i=0;i<17;i++)
	{
		stru_features[index++] = (float)0;
	}
}
short CNFalign_Feat::GenerateIniSample_IX_pure(int currState, int leftState, int pos_x, short* data2,float* data3)
{
	pos_x--;
	int window_size = 5;
	float template_amino_freq = 0;
	short count = 1,index =1;
	data2[0]=0;data3[0]=1;

	//--- first part ---//
	for (int i = 0; i < 20; i++)
	{
		//---- judge ----//
		if (pos_x < 0 || pos_x>length_x-1 )
			template_amino_freq = 1;
		else
			template_amino_freq = (int)(100*temp_p->EmissionProb_original[pos_x][i]);
		//--- record ---//
		if(template_amino_freq!=0){
			data2[index] = count; data3[index++] = template_amino_freq;
		}
		count++;
	}//1-20
	return index;

	//--- second part ---//
	for(int w=-window_size;w<=window_size;w++){
		float ss_score=0;
		if (pos_x+w >= 0 && pos_x+w <= length_x-1) {
			for(int i=0;i<3;i++){
//				ss_score = (int)(100*tseq->SS2[pos_x+w][i]);
				ss_score = 0;
				if(ss_score!=0){
					data2[index]=count; data3[index++]=ss_score;
				}
				count++;
			}
			for(int i=0;i<11;i++){
				count++;
			}
		}
		else{
			count+=14;
		}
	}
	return index;
}
short CNFalign_Feat::GenerateIniSample_IY_pure(int currState, int leftState, int pos_y, short* data2, float* data3)
{
	pos_y--;
	int window_size = 5;
	float seq_amino_freq = 0;
	short count = 1,index = 1;
	data2[0]=0;data3[0]=1;	

	//--- first part ---//
	for (int i = 0; i < 20; i++) {
		//--- judge ----//
		if (pos_y < 0 || pos_y > length_y-1)
			seq_amino_freq = 1;
		else
			seq_amino_freq = (int)(100*seq->EmissionProb_original[pos_y][i]);
		//--- record ---//
		if(seq_amino_freq!=0){
			data2[index] = count; data3[index++] = seq_amino_freq;
		}
		count++;
	}//2-21
	return index;

	//--- second part ---//
	for(int w=-window_size;w<=window_size;w++){
		float ss_score=0;
		if (pos_y+w >= 0 && pos_y+w <= length_y-1) {
			for (int i = 0; i < 3; i++){
//				ss_score = (int)(100*seq->SS2[pos_y+w][i]);
				ss_score = 0;
				if(ss_score!=0){	
					data2[index] = count; data3[index++] = ss_score;
				}
				count++;
			}
			for(int i=0;i<11;i++){
				count++;
			}
		}
		else{
			count+=14;
		}
	}
	return index;
}
void CNFalign_Feat::GenerateOneSample_pure(int currState, int leftState,int pos_x, int pos_y, 
	float* mutation_matrix, float* stru_matrix, float* features) 
{
	int window_size = 5;
	pos_x--;
	pos_y--;

	//--- init ---//
	int index = 0;
	for (int i = 0; i < 20; i++) 
	{
		features[index++] = 0;
	}//1-20

	//--- window size ---//
	for(int w = -window_size;w<=window_size;w++)
	{
		if( pos_x+w<0 || pos_y+w<0 || pos_x+w > length_x-1 || pos_y+w > length_y -1 )
		{
			for(int i=0;i<8;i++)features[index++]= 0;
		}
		else
		{
			for(int i=0;i<8;i++)features[index++]= *(mutation_matrix+((pos_x+w)*(length_y+1)+ (pos_y+w))*8+i);
		}
	}
	int wspair=AApair(pos_x, pos_y);
	if(wspair<0)
	{
		features[index++]=0;
		features[index++]=0;
	}
	else
	{
		features[index++] = (int)(1000*CC50[wspair]);
		features[index++] = (int)(1000*HDSM[wspair]);
	}
		
	//2. add secondary structure score. We can use a vector to represent this score
	for(int w = -window_size;w<=window_size;w++)
	{	
		if( pos_x+w<0 || pos_y+w<0 || pos_x+w > length_x-1 || pos_y+w > length_y -1 )
		{
			for(int i=0;i<26;i++)features[index++]= 0;
		}
		else
		{
			for(int i=0;i<26;i++)features[index++]= *(stru_matrix+((pos_x+w)*(length_y+1)+ (pos_y+w))*26+i);
		}
	}
}
void CNFalign_Feat::GenerateFeatures_pure() 
{
	//match
	for(int i=0; i<=length_x;i++)
	{
		int cur_index=(i-1)*(length_y+1);
		for(int j=0; j<= length_y; j++)
		{
			int t=0;
			if (i - states[t].emit_x >= 0 && j - states[t].emit_y >= 0) 
			{
				GenerateIniSample_Match_pure(i,j,outer_mutation_matrix+(cur_index + (j-1))*8,outer_stru_matrix+(cur_index + (j-1))*26 );
			}
		}
	}
	//IX
	for (int i = 0; i <= length_x; i++)
	{
		if (i - states[1].emit_x >= 0)
		{
			short l = GenerateIniSample_IX_pure(1,0,i,outer_1_feature_data2+i*outer_feat_dim,outer_1_feature_data3+i*outer_feat_dim);
			outer_1_non_zero_length[i]=l;
		}
	}
	//IY
	for (int j = 0; j <= length_y; j++)
	{
		if(j - states[2].emit_y >= 0)
		{
			short l = GenerateIniSample_IY_pure(2,0,j,outer_2_feature_data2+j*outer_feat_dim,outer_2_feature_data3+j*outer_feat_dim);
			outer_2_non_zero_length[j]=l;
		}
	}
	//generate_one_sample
	for (int i = 0; i <= length_x; i++)
	{
		for (int j = 0; j <= length_y; j++)
		{
			int t = 0;
			int s = 0;
			if (i - states[t].emit_x >= 0 && j - states[t].emit_y >= 0)
			{
				if(t==0)
				{
					for(int d=0;d<396;d++)outer_features[d]=0;
					GenerateOneSample_pure(t,s,i,j,outer_mutation_matrix,outer_stru_matrix,outer_features);
					int cur_index=(i*(length_y+1)+j)*outer_feat_dim;
					outer_0_feature_data2[cur_index]=0;
					outer_0_feature_data3[cur_index]=1;
					short count = 1;
					for (short d = 0; d< 396; d++) 
					{
						if(outer_features[d]!= 0 && count<outer_feat_dim)
						{
							outer_0_feature_data2[cur_index+count]=d+1;
							outer_0_feature_data3[cur_index+count]=outer_features[d];
							count++;
						}
					}
					outer_0_non_zero_length[i*(length_y+1)+j]=count;
				}
			}
		}
	}
}

//----- gap and match pre-calculate -----//-> dependent on weight !!!
void CNFalign_Feat::ComputeGap_Gate_Output()
{
	int count;
	count=0;
	for(int i=0;i<=length_x;i++)
	{
		if(i >= states[1].emit_x)
		{
			MI_LogProb[count]=CalculateLogProb_1(1,i,1,0);
			II_LogProb[count]=CalculateLogProb_1(1,i,1,1);
			count++;
		}
	}
	count=0;
	for(int j=0;j<=length_y;j++)
	{
		if(j >= states[2].emit_y)
		{
			MD_LogProb[count]=CalculateLogProb_2(2,1,j,0);
			DD_LogProb[count]=CalculateLogProb_2(2,1,j,2);
			ID_LogProb[count]=CalculateLogProb_2(2,1,j,1);
			count++;
		}
	}
}
void CNFalign_Feat::Compute_Match_Protential()
{
	//----- calculate M_LogProb -----//
	for (int i = 0; i <= length_x; i++) 
	{
		for (int j = 0; j <= length_y; j++) 
		{
			for (int m = 0; m < (int) transitions.size(); m++) 
			{
				int s = transitions[m].from;
				int t = transitions[m].to;
				if (i >= states[t].emit_x && j >= states[t].emit_y) 
				{
					if(t==1)
					{
						if(s==0)
							MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j] = MI_LogProb[i-1];
						else if(s==1)
							MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j] = II_LogProb[i-1];
					}
					else if(t==2)
					{
						if(s==0)
							MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j] = MD_LogProb[j-1];
						else if(s==2)
							MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j] = DD_LogProb[j-1];
						else
							MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j] = ID_LogProb[j-1];
					}
					else if(t==0)
					{
						double aaa = CalculateLogProb_0(t,i,j,s);
						MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j]= aaa;
					}
				}
			}
		}
	}
}


/////===============   gate calculate part =========////
Score CNFalign_Feat::gate(Score x)
{
	return (Score) (1. / (1. + exp(-(double) x)));
}

//======================= calculate log probability related =====================//
//--------- calculate index --------//
int CNFalign_Feat::index_trans(int s, int t) 
{
	int ws=s*3+t;
	switch(ws)
	{
		case 0: return 0;
		case 1: return 1;
		case 2: return 2;
		case 3: return 3;
		case 4: return 4;
		case 5: return 5;
		case 6: return 6;
		case 7: return 8;
		case 8: return 7;
		default:
		{
			fprintf(stderr,"bad transition state here !!\n");
			exit(-1);
		}
	}
}

//------- calculate gate/prob for general purpose ------//
Score CNFalign_Feat::getGateOutput(int x, int y, int m, int g)
{
	int ws_idx = MM_G + m * G_D + g * D_;
	int currState = transitions[m].to;
	Score out = 0;
	switch(currState)
	{
		case 0: 
			out = getGateOutput_0(x, y, m, g,ws_idx);
			break;
		case 1: 
			out = getGateOutput_1(x, y, m, g,ws_idx);
			break;
		case 2: 
			out = getGateOutput_2(x, y, m, g,ws_idx);
			break;
		default:
		{
			fprintf(stderr,"bad transition state at getGateOutput !!\n");
			exit(-1);
		}
	}
	return out;
}
LogScore CNFalign_Feat::CalculateLogProb(int currState, int pos_x,int pos_y, int leftState) 
{
	double new_score=0;
	int s=leftState;
	int t=currState;
	int i=pos_x;
	int j=pos_y;
	int m = index_trans(leftState, currState);
	if (t != 0)  //--- gap state ---//
	{
		if(t==1)
		{
			if(s==0)
			{
				if(i - states[t].emit_x >= 0) new_score = MI_LogProb[i-1];
			}
			else if(s==1)
			{
				if(i - states[t].emit_x >= 0) new_score = II_LogProb[i-1];
			}
		}
		else if(t==2)
		{
			if(s==0)
			{
				if(j - states[t].emit_y >= 0) new_score = MD_LogProb[j-1];
			}
			else if(s==2)
			{
				if(j - states[t].emit_y >= 0) new_score = DD_LogProb[j-1];
			}
			else
			{
				if(j - states[t].emit_y >= 0) new_score = ID_LogProb[j-1];
			}
		}
	}
	else  //--- match state ---//
	{
		new_score = MM_LogProb[m*(length_x+1)*(length_y+1) + i*(length_y+1)+j];
	}
	LogScore prob = new_score;
	return prob;
}


//------- calculate gate output ------// [match:0]
Score CNFalign_Feat::getGateOutput_0(int x, int y, int m, int g,int idx) 
{
	int t = transitions[m].to;
	short* data2_f = inner_0_feature_data2+(x*(length_y+1)+y)*outer_feat_dim;
	float* data3_f = inner_0_feature_data3+(x*(length_y+1)+y)*outer_feat_dim;
	int third_dimension = inner_0_non_zero_length[x*(length_y+1)+y];
	int new_dimension = third_dimension;
	while(new_dimension%4!=0)new_dimension++;
	for(int d =0;d<third_dimension;d++)
	{
		short real_d = data2_f[d];
		feat_vector[d] = data3_f[d];
		weight_vector[d] = weights[idx + real_d];
	}
	for(int d=third_dimension;d<new_dimension;d++)
	{
		feat_vector[d] = (float)0.0;
		weight_vector[d] = (float)0.0;
	}
	Score out2 = fast_dot_product_single(feat_vector, weight_vector, new_dimension);
	if(t==0)out2 *= 0.001;  //match state
	else out2 *= 0.01;      //non-match state
	Score out = gate(out2);
	return out;
}
//------- calculate log prob ------// [match:0]
LogScore CNFalign_Feat::CalculateLogProb_0(int currState, int pos_x,int pos_y, int leftState) 
{
	int ws_m = index_trans(leftState, currState);
	int gnum = _G;
	if (ws_m == 0) gnum = G_;
	int ws_idx = MM_G + ws_m * G_D;
	int ws_idx2 = ws_m * G_;
	int inner_size = G_;
	if(ws_m!=0) inner_size = _G+1;
	for (int g = 0; g < gnum; g++)
	{	
		gate_vector[g] = getGateOutput_0(pos_x, pos_y, ws_m, g,ws_idx);
		gate_weight_vector[g] = weights[ws_idx2+g];
		ws_idx += D_;
	}
	if(ws_m!=0)
	{
		gate_vector[7] = (float)0.0;
		gate_weight_vector[7] = (float)0.0;
	}
	return fast_dot_product_single(gate_vector, gate_weight_vector,inner_size);
}

//------- calculate gate output ------// [ix:1]
Score CNFalign_Feat::getGateOutput_1(int x, int y, int m, int g,int idx) 
{
	int t = transitions[m].to;
	short* data2_f = inner_1_feature_data2+x*outer_feat_dim;
	float* data3_f = inner_1_feature_data3+x*outer_feat_dim;
	int third_dimension = inner_1_non_zero_length[x];
	int new_dimension = third_dimension;
	while(new_dimension%4!=0) new_dimension++;
	for(int d =0;d<third_dimension;d++)
	{
		short real_d = data2_f[d];
		feat_vector[d] = data3_f[d];
		weight_vector[d] = weights[idx + real_d];
	}
	for(int d=third_dimension;d<new_dimension;d++)
	{
		feat_vector[d] = (float)0.0;
		weight_vector[d] = (float)0.0;
	}
	Score out2 = fast_dot_product_single(feat_vector, weight_vector, new_dimension);
	if(t==0)out2 *= 0.001;  //match state
	else out2 *= 0.01;      //non-match state
	Score out = gate(out2);
	return out;
}
//------- calculate log prob ------// [ix:1]
LogScore CNFalign_Feat::CalculateLogProb_1(int currState, int pos_x,int pos_y, int leftState) 
{
	int ws_m = index_trans(leftState, currState);
	int gnum = _G;
	if (ws_m == 0) gnum = G_;
	int ws_idx = MM_G + ws_m * G_D;
	int ws_idx2 = ws_m * G_;
	int inner_size = G_;
	if(ws_m!=0) inner_size = _G+1;
	for (int g = 0; g < gnum; g++)
	{	
		gate_vector[g] = getGateOutput_1(pos_x, pos_y, ws_m, g,ws_idx);
		gate_weight_vector[g] = weights[ws_idx2+g];
		ws_idx += D_;
	}
	if(ws_m!=0)
	{
		gate_vector[7] = (float)0.0;
		gate_weight_vector[7] = (float)0.0;
	}
	return fast_dot_product_single(gate_vector, gate_weight_vector,inner_size);
}

//------- calculate gate output ------// [iy:2]
Score CNFalign_Feat::getGateOutput_2(int x, int y, int m, int g,int idx) 
{
	int t = transitions[m].to;
	short* data2_f = inner_2_feature_data2+y*outer_feat_dim;
	float* data3_f = inner_2_feature_data3+y*outer_feat_dim;
	int third_dimension = inner_2_non_zero_length[y];
	int new_dimension = third_dimension;
	while(new_dimension%4!=0) new_dimension++;
	for(int d =0;d<third_dimension;d++)
	{
		short real_d = data2_f[d];
		feat_vector[d] = data3_f[d];
		weight_vector[d] = weights[idx + real_d];
	}
	for(int d=third_dimension;d<new_dimension;d++)
	{
		feat_vector[d] = (float)0.0;
		weight_vector[d] = (float)0.0;
	}
	Score out2 = fast_dot_product_single(feat_vector, weight_vector, new_dimension);
	if(t==0)out2 *= 0.001;  //match state
	else out2 *= 0.01;      //non-match state
	Score out = gate(out2);
	return out;
}
//------- calculate log prob ------// [iy:2]
LogScore CNFalign_Feat::CalculateLogProb_2(int currState, int pos_x,int pos_y, int leftState) 
{
	int ws_m = index_trans(leftState, currState);
	int gnum = _G;
	if (ws_m == 0) gnum = G_;
	int ws_idx = MM_G + ws_m * G_D;
	int ws_idx2 = ws_m * G_;
	int inner_size = G_;
	if(ws_m!=0) inner_size = _G+1;
	for (int g = 0; g < gnum; g++)
	{	
		gate_vector[g] = getGateOutput_2(pos_x, pos_y, ws_m, g,ws_idx);
		gate_weight_vector[g] = weights[ws_idx2+g];
		ws_idx += D_;
	}
	if(ws_m!=0)
	{
		gate_vector[7] = (float)0.0;
		gate_weight_vector[7] = (float)0.0;
	}
	return fast_dot_product_single(gate_vector, gate_weight_vector,inner_size);
}

/*
//============= states and transitions ========//
vector<State> states;           //the states describing an alignment
vector<Transition> transitions; //the state transition matrix
void InitializeStatesAndTransitions() 
{
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
	for (int i = 0; i < states.size(); i++)
		for (int j = 0; j < states.size(); j++) {
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
*/

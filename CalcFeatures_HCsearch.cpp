#include "CalcFeatures_HCsearch.h"

namespace CalcFeature {
	//---- define ----//
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define GREAT(T, score) ( (score)>(T)? (score):0)
#define SQ(a) ((a)*(a))
#define LONG_GAP	40      // An insertion segment is long if it has more than 30 residues
#define PROB_BASE 100.0   // probability


	//======================= I/O related ==========================//
	//-------- utility ------//
	void getBaseName(string &in,string &out,char slash,char dot)
	{
		int i,j;
		int len=(int)in.length();
		for(i=len-1;i>=0;i--)
		{
			if(in[i]==slash)break;
		}
		i++;
		for(j=len-1;j>=0;j--)
		{
			if(in[j]==dot)break;
		}
		if(j==-1)j=len;
		out=in.substr(i,j-i);
	}
	void getRootName(string &in,string &out,char slash)
	{
		int i;
		int len=(int)in.length();
		for(i=len-1;i>=0;i--)
		{
			if(in[i]==slash)break;
		}
		if(i<=0)out=".";
		else out=in.substr(0,i);
	}

	//--------- FASTA I/O ------------//
	//FASTA
	int ReadToFile_FASTA(string &fn,vector<pair<int, int> > &alignment,
			string &nam1_content,string &nam2_content,
			string &nam1_full,string &nam2_full,
			string &nam1,string &nam2)
	{
		int i;
		int cur1=0;
		int cur2=0;
		int len;
		int len1,len2;
		alignment.clear();
		//init
		string seq="";  //sequence
		string tmp="";  //template
		//load
		ifstream fin;
		string buf,temp;
		fin.open(fn.c_str(), ios::in);
		if(fin.fail()!=0)
		{
			fprintf(stderr,"alignment file not found [%s] !!!\n",fn.c_str());
			return -1;
		}
		//read tmp
		for(;;)
		{
			if(!getline(fin,buf,'\n'))goto badend;
			len=(int)buf.length();
			if(len>1)
			{
				if(buf[0]=='>')
				{
					istringstream www(buf);
					www>>temp;
					len=(int)temp.length();
					nam1=temp.substr(1,len-1);
					break;
				}
			}
		}
		for(;;)
		{
			if(!getline(fin,buf,'\n'))goto badend;
			len=(int)buf.length();
			if(len==0)continue;
			if(len>1)
			{
				if(buf[0]=='>')
				{
					istringstream www(buf);
					www>>temp;
					len=(int)temp.length();
					nam2=temp.substr(1,len-1);
					break;
				}
			}
			tmp+=buf;
		}
		//read seq
		for(;;)
		{
			if(!getline(fin,buf,'\n'))break;
			len=(int)buf.length();
			if(len==0)continue;
			seq+=buf;
		}
		//process
		len1=(int)seq.length();
		len2=(int)tmp.length();
		if(len1!=len2)
		{
			fprintf(stderr,"alignment len not equal [%s] !!!\n",fn.c_str());
			return -1;
		}
		len=len1;
		nam1_content.clear();
		nam2_content.clear();
		for(i=0;i<len;i++)
		{
			if(tmp[i]!='-' && seq[i]!='-') //match
			{
				nam1_content.push_back(tmp[i]);
				nam2_content.push_back(seq[i]);
				cur1++;
				cur2++;
				alignment.push_back(pair<int,int>(cur1,cur2));
			}
			else
			{
				if(tmp[i]!='-') //Ix
				{
					nam1_content.push_back(tmp[i]);
					cur1++;
					alignment.push_back(pair<int,int>(cur1,-cur2));
				}
				if(seq[i]!='-') //Iy
				{
					nam2_content.push_back(seq[i]);
					cur2++;
					alignment.push_back(pair<int,int>(-cur1,cur2));
				}
			}
		}
		//return
		nam1_full=tmp;
		nam2_full=seq;
		return 1; //success

badend:
		fprintf(stderr,"alignment file format bad [%s] !!!\n",fn.c_str());
		return -1;
	}


	//====================================== calculate basic features =======================//
	//=================//
	//--Ori_BLOSUM----//
	//===============//
#define B62_OK		-1.01
#define B62_GOOD	0
#define B62_GREAT	1.01

	int Ori_BLOSUM_62_MD[21][21]={
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

#define B45_OK		-1.01
#define B45_GOOD	0
#define B45_GREAT	1.01

	int Ori_BLOSUM_45_MD[21][21]={
		{  5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -5 },  //A
		{ -2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -5 },  //R
		{ -1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3, -5 },  //N
		{ -2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3, -5 },  //D
		{ -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -5 },  //C
		{ -1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3, -5 },  //Q
		{ -1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3, -5 },  //E
		{  0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -5 },  //G
		{ -2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3, -5 },  //H
		{ -1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -5 },  //I
		{ -1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -5 },  //L
		{ -1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2, -5 },  //K
		{ -1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -5 },  //M
		{ -2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -5 },  //F
		{ -1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -5 },  //P
		{  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1, -5 },  //S
		{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0, -5 },  //T
		{ -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -5 },  //W
		{ -2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -5 },  //Y
		{  0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -5 },  //V
		{ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 }}; //Z
	// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z  

#define B80_OK 		-1.01
#define B80_GOOD	0
#define B80_GREAT	1.01

	int Ori_BLOSUM_80_MD[21][21]={
		{  5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -5 },  //A
		{ -2,  6, -1, -2, -4,  1, -1, -3,  0, -3, -3,  2, -2, -4, -2, -1, -1, -4, -3, -3, -5 },  //R
		{ -2, -1,  6,  1, -3,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -4, -3, -4, -5 },  //N
		{ -2, -2,  1,  6, -4, -1,  1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4, -5 },  //D
		{ -1, -4, -3, -4,  9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -5 },  //C
		{ -1,  1,  0, -1, -4,  6,  2, -2,  1, -3, -3,  1,  0, -4, -2,  0, -1, -3, -2, -3, -5 },  //Q
		{ -1, -1, -1,  1, -5,  2,  6, -3,  0, -4, -4,  1, -2, -4, -2,  0, -1, -4, -3, -3, -5 },  //E
		{  0, -3, -1, -2, -4, -2, -3,  6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -5 },  //G
		{ -2,  0,  0, -2, -4,  1,  0, -3,  8, -4, -3, -1, -2, -2, -3, -1, -2, -3,  2, -4, -5 },  //H
		{ -2, -3, -4, -4, -2, -3, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -5 },  //I
		{ -2, -3, -4, -5, -2, -3, -4, -4, -3,  1,  4, -3,  2,  0, -3, -3, -2, -2, -2,  1, -5 },  //L
		{ -1,  2,  0, -1, -4,  1,  1, -2, -1, -3, -3,  5, -2, -4, -1, -1, -1, -4, -3, -3, -5 },  //K
		{ -1, -2, -3, -4, -2,  0, -2, -4, -2,  1,  2, -2,  6,  0, -3, -2, -1, -2, -2,  1, -5 },  //M
		{ -3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4,  0,  6, -4, -3, -2,  0,  3, -1, -5 },  //F
		{ -1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4,  8, -1, -2, -5, -4, -3, -5 },  //P
		{  1, -1,  0, -1, -2,  0,  0, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2, -5 },  //S
		{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2,  1,  5, -4, -2,  0, -5 },  //T
		{ -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2,  0, -5, -4, -4, 11,  2, -3, -5 },  //W
		{ -2, -3, -3, -4, -3, -2, -3, -4,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -5 },  //Y
		{  0, -3, -4, -4, -1, -3, -3, -4, -4,  3,  1, -3,  1, -1, -3, -2,  0, -3, -2,  4, -5 },  //V
		{ -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5 }}; //Z
	// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   Z  


	//BLOSUM_Mapping//--------------ARNDCQEGHILKMFPSTWYVZ
	int Blo_AA_Map_MD[21]=
	{ 0,19, 4, 3, 6, 13,7, 8, 9, 17,11,10,12,2, 18,14,5, 1, 15,16,20};
	//A  V  C  D  E  F  G  H  I  W  K  L  M  N  Y  P  Q  R   S  T  Z
	//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17  18 19 20
	//Ori_Mapping//-----------------AVCDEFGHIWKLMNYPQRSTZ
	int Ori_AA_Map_MD[26]=
	{ 0,20,2,3,4,5,6,7,8,20,10,11,12,13,20,15,16,17,18,19,20, 1, 9,20,14,20};
	// A B C D E F G H I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
	// 0 1 2 3 4 5 6 7 8  9 10 11 12 14 14 15 16 17 18 19 20 21 22 23 24 25

	//------ calculate -------//
	int BLOSUM62_Calc(char a,char b)
	{
		int ii,jj;
		//if(a<'A' || a>'Z')a='Z';
		if(a<'A' || a>'Z') {
			cerr<<"a non standard amino acid: "<<a<<endl;
			return 0; //a non-standard amino acid, neither reward nor penalize
		}
		ii=Blo_AA_Map_MD[Ori_AA_Map_MD[a-'A']];
		//if(b<'A' || b>'Z')b='Z';
		if(b<'A' || b>'Z') {
			cerr<<"a non standard amino acid: "<<b<<endl;
			return 0; //a non-standard amino acid
		}
		jj=Blo_AA_Map_MD[Ori_AA_Map_MD[b-'A']];
		return Ori_BLOSUM_62_MD[ii][jj];
	}
	int BLOSUM45_Calc(char a,char b)
	{
		int ii,jj;
		//if(a<'A' || a>'Z')a='Z';
		if(a<'A' || a>'Z') {
			cerr<<"a non standard amino acid: "<<a<<endl;
			return 0; //a non-standard amino acid, neither reward nor penalize
		}
		ii=Blo_AA_Map_MD[Ori_AA_Map_MD[a-'A']];
		//if(b<'A' || b>'Z')b='Z';
		if(b<'A' || b>'Z') {
			cerr<<"a non standard amino acid: "<<b<<endl;
			return 0; //a non-standard amino acid
		}
		jj=Blo_AA_Map_MD[Ori_AA_Map_MD[b-'A']];
		return Ori_BLOSUM_45_MD[ii][jj];
	}
	int BLOSUM80_Calc(char a,char b)
	{
		int ii,jj;
		//if(a<'A' || a>'Z')a='Z';
		if(a<'A' || a>'Z') {
			cerr<<"a non standard amino acid: "<<a<<endl;
			return 0; //a non-standard amino acid, neither reward nor penalize
		}
		ii=Blo_AA_Map_MD[Ori_AA_Map_MD[a-'A']];
		//if(b<'A' || b>'Z')b='Z';
		if(b<'A' || b>'Z') {
			cerr<<"a non standard amino acid: "<<b<<endl;
			return 0; //a non-standard amino acid
		}
		jj=Blo_AA_Map_MD[Ori_AA_Map_MD[b-'A']];
		return Ori_BLOSUM_80_MD[ii][jj];
	}

	//------------- two important matrix ----------------//
#define HDSM_OK		-0.8
#define HDSM_GOOD	0
#define HDSM_GREAT	0.8

	double Ori_HDSM_MD[20][20]={
		{  2.09,  -0.50,  -0.57,  -0.73,   0.33,  -0.75,  -0.12,   0.27,  -1.42,  -0.97,  -0.39,  -0.38,  -0.04,  -0.76,  -0.53,   0.34,   0.13,  -0.66,  -1.25,   0.02}, //A
		{ -0.50,   2.87,   0.60,   0.13,  -1.30,   0.13,   0.99,  -0.96,   0.54,  -1.40,  -1.19,   1.42,  -0.63,  -1.40,   0.21,  -0.06,  -0.15,  -0.04,  -0.75,  -1.52}, //R
		{ -0.57,   0.60,   3.60,   1.78,  -2.08,   0.33,  -0.16,   0.79,   0.76,  -2.43,  -2.10,   0.83,  -2.01,  -2.25,  -1.10,   0.40,   0.30,  -2.89,  -0.36,  -2.17}, //N
		{ -0.73,   0.13,   1.78,   4.02,  -2.51,   0.34,   1.20,  -1.20,  -0.01,  -2.77,  -2.65,   0.66,  -2.58,  -2.19,   0.72,   0.71,  -0.75,  -1.91,  -1.21,  -2.02}, //D
		{  0.33,  -1.30,  -2.08,  -2.51,   6.99,  -0.83,  -1.97,  -2.11,  -1.50,   0.13,  -0.31,  -2.19,   1.04,   1.13,  -2.19,   0.31,  -0.59,  -0.76,   0.13,   0.34}, //C
		{ -0.75,   0.13,   0.33,   0.34,  -0.83,   2.60,   1.23,  -0.12,  -0.46,  -1.47,  -1.49,   0.92,  -0.13,  -2.31,   0.24,   1.04,   0.60,  -0.81,  -0.61,  -1.38}, //Q
		{ -0.12,   0.99,  -0.16,   1.20,  -1.97,   1.23,   2.97,  -0.41,  -0.62,  -1.81,  -2.11,   1.11,  -1.86,  -1.61,  -0.26,   0.31,  -0.21,  -2.70,  -1.64,  -1.84}, //E
		{  0.27,  -0.96,   0.79,  -1.20,  -2.11,  -0.12,  -0.41,   4.36,  -0.40,  -2.93,  -1.98,  -0.71,  -1.86,  -2.67,  -0.04,   0.29,  -0.81,  -1.21,  -1.62,  -1.96}, //G
		{ -1.42,   0.54,   0.76,  -0.01,  -1.50,  -0.46,  -0.62,  -0.40,   5.89,  -1.76,  -0.93,   0.31,  -1.04,  -0.22,  -1.44,  -0.74,  -0.52,  -1.48,  -0.12,  -0.35}, //H
		{ -0.97,  -1.40,  -2.43,  -2.77,   0.13,  -1.47,  -1.81,  -2.93,  -1.76,   2.76,   1.56,  -1.81,   0.99,   0.76,  -2.00,  -1.75,  -0.96,   0.25,   0.08,   1.94}, //I
		{ -0.39,  -1.19,  -2.10,  -2.65,  -0.31,  -1.49,  -2.11,  -1.98,  -0.93,   1.56,   2.43,  -1.96,   1.61,   1.23,  -1.56,  -2.30,  -0.86,  -0.14,   0.70,   0.81}, //L
		{ -0.38,   1.42,   0.83,   0.66,  -2.19,   0.92,   1.11,  -0.71,   0.31,  -1.81,  -1.96,   2.91,  -1.62,  -2.41,  -0.19,  -0.06,  -0.10,  -1.94,  -1.72,  -1.27}, //K
		{ -0.04,  -0.63,  -2.01,  -2.58,   1.04,  -0.13,  -1.86,  -1.86,  -1.04,   0.99,   1.61,  -1.62,   3.75,   0.80,  -1.09,  -1.34,  -1.58,   0.87,  -0.41,   0.61}, //M
		{ -0.76,  -1.40,  -2.25,  -2.19,   1.13,  -2.31,  -1.61,  -2.67,  -0.22,   0.76,   1.23,  -2.41,   0.80,   3.28,  -0.91,  -1.11,  -0.69,   2.29,   1.96,   0.51}, //F
		{ -0.53,   0.21,  -1.10,   0.72,  -2.19,   0.24,  -0.26,  -0.04,  -1.44,  -2.00,  -1.56,  -0.19,  -1.09,  -0.91,   5.45,  -0.29,   0.93,  -5.34,  -1.98,  -1.11}, //P
		{  0.34,  -0.06,   0.40,   0.71,   0.31,   1.04,   0.31,   0.29,  -0.74,  -1.75,  -2.30,  -0.06,  -1.34,  -1.11,  -0.29,   2.36,   1.20,  -1.18,  -1.56,  -1.11}, //S
		{  0.13,  -0.15,   0.30,  -0.75,  -0.59,   0.60,  -0.21,  -0.81,  -0.52,  -0.96,  -0.86,  -0.10,  -1.58,  -0.69,   0.93,   1.20,   2.04,  -0.57,  -0.41,   0.05}, //T
		{ -0.66,  -0.04,  -2.89,  -1.91,  -0.76,  -0.81,  -2.70,  -1.21,  -1.48,   0.25,  -0.14,  -1.94,   0.87,   2.29,  -5.34,  -1.18,  -0.57,   6.96,   2.15,  -1.09}, //W
		{ -1.25,  -0.75,  -0.36,  -1.21,   0.13,  -0.61,  -1.64,  -1.62,  -0.12,   0.08,   0.70,  -1.72,  -0.41,   1.96,  -1.98,  -1.56,  -0.41,   2.15,   3.95,   0.21}, //Y
		{  0.02,  -1.52,  -2.17,  -2.02,   0.34,  -1.38,  -1.84,  -1.96,  -0.35,   1.94,   0.81,  -1.27,   0.61,   0.51,  -1.11,  -1.11,   0.05,  -1.09,   0.21,   2.05}};//V
	//  A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W       Y       V       

#define CC_OK 	-0.25
#define CC_GOOD	0
#define CC_GREAT	0.25

	double Ori_CC50_MD[20][20]={
		{1.000,0.620,0.257,0.133,0.411,0.684,0.681,0.368,0.361,0.481,0.728,0.586,0.735,0.488,0.431,0.355,0.440,0.475,0.479,0.484}, //A
		{0.620,1.000,0.407,0.283,0.303,0.862,0.630,0.410,0.679,0.203,0.525,0.929,0.739,0.683,0.399,0.689,0.655,0.659,0.594,0.208}, //R
		{0.257,0.407,1.000,0.872,0.325,0.465,0.404,0.463,0.692,0.146,0.134,0.410,0.284,0.197,0.411,0.790,0.403,0.320,0.315,0.155}, //N
		{0.133,0.283,0.872,1.000,0.215,0.383,0.466,0.362,0.674,0.034,0.026,0.535,0.155,0.076,0.302,0.707,0.300,0.202,0.201,0.046}, //D
		{0.411,0.303,0.325,0.215,1.000,0.389,0.255,0.414,0.861,0.365,0.356,0.275,0.428,0.651,0.426,0.825,0.742,0.696,0.673,0.368}, //C
		{0.684,0.862,0.465,0.383,0.389,1.000,0.740,0.477,0.711,0.270,0.710,0.859,0.626,0.640,0.465,0.475,0.711,0.657,0.655,0.278}, //Q
		{0.681,0.630,0.404,0.466,0.255,0.740,1.000,0.370,0.362,0.157,0.401,0.743,0.245,0.184,0.363,0.396,0.359,0.288,0.297,0.168}, //E
		{0.368,0.410,0.463,0.362,0.414,0.477,0.370,1.000,0.468,0.282,0.273,0.394,0.396,0.333,0.465,0.470,0.455,0.423,0.412,0.288}, //G
		{0.361,0.679,0.692,0.674,0.861,0.711,0.362,0.468,1.000,0.276,0.265,0.416,0.387,0.567,0.454,0.930,0.704,0.661,0.738,0.285}, //H
		{0.481,0.203,0.146,0.034,0.365,0.270,0.157,0.282,0.276,1.000,0.499,0.168,0.465,0.491,0.364,0.263,0.380,0.443,0.449,0.948}, //I
		{0.728,0.525,0.134,0.026,0.356,0.710,0.401,0.273,0.265,0.499,1.000,0.606,0.791,0.741,0.358,0.252,0.371,0.688,0.447,0.498}, //L
		{0.586,0.929,0.410,0.535,0.275,0.859,0.743,0.394,0.416,0.168,0.606,1.000,0.512,0.444,0.363,0.687,0.643,0.292,0.307,0.174}, //K
		{0.735,0.739,0.284,0.155,0.428,0.626,0.245,0.396,0.387,0.465,0.791,0.512,1.000,0.815,0.445,0.366,0.438,0.821,0.487,0.466}, //M
		{0.488,0.683,0.197,0.076,0.651,0.640,0.184,0.333,0.567,0.491,0.741,0.444,0.815,1.000,0.400,0.550,0.650,0.918,0.919,0.492}, //F
		{0.431,0.399,0.411,0.302,0.426,0.465,0.363,0.465,0.454,0.364,0.358,0.363,0.445,0.400,1.000,0.449,0.463,0.470,0.470,0.370}, //P
		{0.355,0.689,0.790,0.707,0.825,0.475,0.396,0.470,0.930,0.263,0.252,0.687,0.366,0.550,0.449,1.000,0.791,0.646,0.724,0.271}, //S
		{0.440,0.655,0.403,0.300,0.742,0.711,0.359,0.455,0.704,0.380,0.371,0.643,0.438,0.650,0.463,0.791,1.000,0.699,0.904,0.387}, //T
		{0.475,0.659,0.320,0.202,0.696,0.657,0.288,0.423,0.661,0.443,0.688,0.292,0.821,0.918,0.470,0.646,0.699,1.000,0.742,0.444}, //W
		{0.479,0.594,0.315,0.201,0.673,0.655,0.297,0.412,0.738,0.449,0.447,0.307,0.487,0.919,0.470,0.724,0.904,0.742,1.000,0.452}, //Y
		{0.484,0.208,0.155,0.046,0.368,0.278,0.168,0.288,0.285,0.948,0.498,0.174,0.466,0.492,0.370,0.271,0.387,0.444,0.452,1.000}};//V
	//A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V      

	//============= old feature related ===========//
	int Ori_GONNET_MD[20][20]={
		{   24,   -6,   -3,   -3,    5,   -2,    0,    5,   -8,   -8,  -12,   -4,   -7,  -23,    3,   11,    6,  -36,  -22,    1},  //A
		{   -6,   47,    3,   -3,  -22,   15,    4,  -10,    6,  -24,  -22,   27,  -17,  -32,   -9,   -2,   -2,  -16,  -18,  -20},  //R
		{   -3,    3,   38,   22,  -18,    7,    9,    4,   12,  -28,  -30,    8,  -22,  -31,   -9,    9,    5,  -36,  -14,  -22},  //N
		{   -3,   -3,   22,   47,  -32,    9,   27,    1,    4,  -38,  -40,    5,  -30,  -45,   -7,    5,    0,  -52,  -28,  -29},  //D
		{    5,  -22,  -18,  -32,  115,  -24,  -30,  -20,  -13,  -11,  -15,  -28,   -9,   -8,  -31,    1,   -5,  -10,   -5,    0},  //C
		{   -2,   15,    7,    9,  -24,   27,   17,  -10,   12,  -19,  -16,   15,  -10,  -26,   -2,    2,    0,  -27,  -17,  -15},  //Q
		{    0,    4,    9,   27,  -30,   17,   36,   -8,    4,  -27,  -28,   12,  -20,  -39,   -5,    2,   -1,  -43,  -27,  -19},  //E
		{    5,  -10,    4,    1,  -20,  -10,   -8,   66,  -14,  -45,  -44,  -11,  -35,  -52,  -16,    4,  -11,  -40,  -40,  -33},  //G
		{   -8,    6,   12,    4,  -13,   12,    4,  -14,   60,  -22,  -19,    6,  -13,   -1,  -11,   -2,   -3,   -8,   22,  -20},  //H
		{   -8,  -24,  -28,  -38,  -11,  -19,  -27,  -45,  -22,   40,   28,  -21,   25,   10,  -26,  -18,   -6,  -18,   -7,   31},  //I
		{  -12,  -22,  -30,  -40,  -15,  -16,  -28,  -44,  -19,   28,   40,  -21,   28,   20,  -23,  -21,  -13,   -7,    0,   18},  //L
		{   -4,   27,    8,    5,  -28,   15,   12,  -11,    6,  -21,  -21,   32,  -14,  -33,   -6,    1,    1,  -35,  -21,  -17},  //K
		{   -7,  -17,  -22,  -30,   -9,  -10,  -20,  -35,  -13,   25,   28,  -14,   43,   16,  -24,  -14,   -6,  -10,   -2,   16},  //M
		{  -23,  -32,  -31,  -45,   -8,  -26,  -39,  -52,   -1,   10,   20,  -33,   16,   70,  -38,  -28,  -22,   36,   51,    1},  //F
		{    3,   -9,   -9,   -7,  -31,   -2,   -5,  -16,  -11,  -26,  -23,   -6,  -24,  -38,   76,    4,    1,  -50,  -31,  -18},  //P
		{   11,   -2,    9,    5,    1,    2,    2,    4,   -2,  -18,  -21,    1,  -14,  -28,    4,   22,   15,  -33,  -19,  -10},  //S
		{    6,   -2,    5,    0,   -5,    0,   -1,  -11,   -3,   -6,  -13,    1,   -6,  -22,    1,   15,   25,  -35,  -19,    0},  //T
		{  -36,  -16,  -36,  -52,  -10,  -27,  -43,  -40,   -8,  -18,   -7,  -35,  -10,   36,  -50,  -33,  -35,  142,   41,  -26},  //W
		{  -22,  -18,  -14,  -28,   -5,  -17,  -27,  -40,   22,   -7,    0,  -21,   -2,   51,  -31,  -19,  -19,   41,   78,  -11},  //Y
		{    1,  -20,  -22,  -29,    0,  -15,  -19,  -33,  -20,   31,   18,  -17,   16,    1,  -18,  -10,    0,  -26,  -11,   34}}; //V
	//   A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V      

	int Ori_Singleton_MD[20][9]={
		{  -60,  -13,  -18,    0,   58,   94,    1,   21,   33},   //A
		{   99,  -54,  -53,  129,  -34,  -58,   94,   -3,   -5},   //R
		{   87,   13,    2,   88,   25,    5,    6,  -31,  -46},   //N
		{  111,   21,  -38,  117,   33,    7,   34,  -20,  -55},   //D
		{  -33,   41,  191,  -66,    6,  121,  -70,  -23,  124},   //C
		{  107,  -32,  -97,  146,   12,  -57,  126,   43,  -29},   //Q
		{   67,  -34,  -77,  101,    4,  -51,   83,   22,  -16},   //E
		{   47,  104,  116,   19,   58,  101,  -50,  -25,  -54},   //G
		{   51,  -18,   18,   33,  -35,   -3,    3,  -30,    4},   //H
		{  -58,    5,  119,  -89,  -20,   50,  -18,   38,  135},   //I
		{  -75,  -23,   91,  -41,   19,   91,  -26,   15,  120},   //L
		{  184,   -7,  -91,  210,   -4,  -90,  196,   45,  -51},   //K
		{  -65,  -20,   74,  -27,   20,   65,  -25,   11,   74},   //M
		{  -47,    9,  139,  -64,  -19,   79,  -44,   -3,  124},   //F
		{  126,   91,   46,  131,   76,   20,  -35,  -46,  -63},   //P
		{   37,   31,    8,   32,    6,   -8,  -14,  -19,  -25},   //S
		{   29,   19,   29,   17,  -38,  -60,    1,   -8,   -9},   //T
		{  -33,  -34,  123,  -26,  -48,   65,  -32,  -15,  123},   //W
		{   -9,  -26,   97,  -26,  -70,   28,   -7,  -22,   93},   //Y
		{  -40,   19,  109,  -93,  -35,   10,   -4,   32,  101}};  //V
	// bury  medi  expo  bury  medi  expo   bury medi  expo
	//         HELIX            SHEET            LOOP

	//==== TMscore > 0.35 ======//
	////singleton matrix (substition matrix from 17W CLEF alignment) -> (more positive, better)
	double WS_Singleton_MD[20][9]={
		{ 0.4521,  0.1993,  0.1440, -0.1481, -0.4385, -0.5268, -0.0161, -0.1930, -0.2784},  //A
		{-0.3654,  0.3570,  0.4287, -0.6435,  0.2095,  0.2655, -0.3786,  0.0195,  0.1137},  //R
		{-0.5198, -0.0223,  0.1781, -0.6236, -0.1045,  0.0789,  0.0198,  0.2663,  0.3942},  //N
		{-0.6642,  0.0034,  0.3186, -0.6616, -0.1438,  0.0117, -0.0439,  0.2922,  0.3959},  //D
		{ 0.1426, -0.4187, -0.9984,  0.4265,  0.0323, -0.3947,  0.4247,  0.0587, -0.5173},  //C
		{-0.2788,  0.3363,  0.5703, -0.6899,  0.0179,  0.2306, -0.4087, -0.0564,  0.1450},  //Q
		{-0.4877,  0.3586,  0.7434, -0.8903, -0.0308,  0.3240, -0.5427, -0.1179,  0.2259},  //E
		{-0.4438, -0.5533, -0.5185, -0.1935, -0.2888, -0.2710,  0.3172,  0.3411,  0.4983},  //G
		{-0.2226,  0.0608, -0.1088, -0.1423,  0.1681,  0.0601,  0.1487,  0.1427, -0.0238},  //H
		{ 0.3407, -0.2074, -0.6775,  0.6630,  0.0352, -0.3408,  0.0132, -0.3905, -0.7773},  //I
		{ 0.5600,  0.0477, -0.3708,  0.2819, -0.2322, -0.5238,  0.0296, -0.2958, -0.6507},  //L
		{-0.5763,  0.2980,  0.5645, -0.9008,  0.1691,  0.4039, -0.5606,  0.0123,  0.3085},  //K
		{ 0.4432,  0.0740, -0.2794,  0.1453, -0.2255, -0.4119,  0.0810, -0.1618, -0.4200},  //M
		{ 0.3000, -0.1884, -0.7094,  0.4458,  0.0652, -0.3039,  0.1840, -0.1604, -0.5609},  //F
		{-0.7645, -0.4526, -0.2861, -0.6109, -0.1524,  0.2258,  0.3069,  0.3887,  0.5543},  //P
		{-0.3276, -0.1341,  0.0339, -0.3139,  0.0630,  0.2579,  0.0899,  0.1754,  0.2426},  //S
		{-0.2333, -0.1518, -0.2018, -0.0647,  0.3434,  0.5032,  0.0151,  0.1004,  0.0398},  //T
		{ 0.2079, -0.0059, -0.4754,  0.2423,  0.1374, -0.1509,  0.1852, -0.0916, -0.4484},  //W
		{ 0.0997, -0.0204, -0.4840,  0.2422,  0.2807, -0.1224,  0.0971, -0.0134, -0.3813},  //Y
		{ 0.1356, -0.3279, -0.7181,  0.7086,  0.2294, -0.0558,  0.0044, -0.3293, -0.6776}}; //V 
	// bury     medi     expo     bury     medi     expo     bury     medi     expo
	//         HELIX                       SHEET                      LOOP

	double Ori_Contact_MD[20][27]={
		{ 195.9, 309.8, -59.8,   4.0, 174.6,  -9.7,  -3.2, 109.2,  22.7,  37.2,  57.8,  10.3,  15.9,  10.6,   2.7, -27.8, -23.8,   7.9, -51.4, -68.2,  13.7, -35.1,-102.0,  15.5,  -3.1,-129.0,  13.8},  //A
		{ 168.5, 282.7, -65.9, -19.2,  77.2, -43.1, -58.4,   4.5, -23.1, -45.0, -33.0, -16.0, -16.3, -37.2,  13.4,  24.7, -30.0,  52.4,  83.4, -15.9,  94.2, 153.7,   4.3, 133.4, 237.9,  52.9, 156.8},  //R
		{ 152.9, 195.5, -82.8, -35.9,  63.7, -35.8, -55.3,  14.2,  -9.1, -12.3, -17.5,  -2.8,  -0.8, -21.5,  19.2,  17.4, -17.0,  48.1,  43.7, -17.8,  71.7,  60.0, -32.9,  84.1, 112.3, -48.9,  83.1},  //N
		{ 111.1, 185.9,-100.2, -72.2,  33.0, -49.5, -69.2, -15.8, -15.4,  -8.8, -27.5,   5.8,  31.4, -16.1,  34.9,  58.9, -13.7,  72.0,  95.3,  -0.4, 102.8, 128.8,  -1.1, 128.1, 150.0,   7.5, 133.1},  //D
		{ 314.1, 568.7,  70.8, 148.4, 302.2,  87.6,  99.3, 191.6,  55.8,  48.4,  94.3,   9.0,  -0.7,  27.5, -24.4, -56.2, -13.2, -27.6, -78.9, -71.6, -34.8, -78.8,-115.2, -42.5, -53.7,-168.0, -64.1},  //C
		{ 144.8, 245.8, -84.1, -42.4,  55.5, -50.1, -65.1,  -4.5, -16.1, -29.2, -28.7,  -5.8,   7.2, -25.4,  30.4,  37.5, -19.3,  63.7,  79.6, -16.5,  98.0, 106.4,  -0.8, 105.6, 167.9,   9.5, 126.5},  //Q
		{  97.9, 205.4,-113.3, -74.0,  35.7, -58.6, -70.9, -25.7, -17.6, -17.8, -44.8,   9.4,  31.7, -23.8,  52.9,  68.6,  -4.0, 103.9, 123.0,  16.4, 139.2, 173.1,  28.0, 161.8, 227.3,  50.7, 168.1},  //E
		{ 136.6, 204.3, -86.5,   5.1,  96.4, -31.0,  -1.2,  96.1,  12.8,  35.3,  65.3,  19.8,  26.1,  15.5,  28.9,  -9.4, -38.4,  32.7, -46.8, -67.9,  34.1, -47.6, -88.8,  32.7, -52.0,-105.8,  17.9},  //G
		{ 156.3, 305.5, -57.9,   6.1,  89.0, -12.4, -39.4,  32.2, -11.0, -30.8,  -5.4, -12.9, -28.3, -29.0,  -0.2,  -8.0, -36.2,  19.9,  34.7, -35.2,  55.6,  81.4, -40.3,  67.0, 158.3, -29.2,  82.2},  //H
		{ 257.9, 484.1,  15.5, 101.1, 232.0,  54.3,  39.6, 124.7,  40.4,  16.6,  53.2,  -9.1, -18.9,  11.7, -19.7, -47.8, -31.1, -16.0, -53.0, -69.8, -20.8, -51.3,-100.8, -20.4, -20.6,-131.7, -23.8},  //I
		{ 287.8, 449.7,  14.6,  87.0, 242.9,  43.4,  26.3, 122.3,  27.0,   1.6,  44.7, -16.8, -32.5,  -2.5, -23.1, -48.9, -41.0, -15.3, -40.7, -67.1,  -9.9, -22.7, -91.8,  -5.4,  21.1,-115.7,  10.6},  //L
		{ 110.0, 237.7,-100.8, -54.4,  42.8, -54.4, -70.0, -20.5, -25.2, -35.9, -44.9,  -7.5,   4.0, -39.2,  39.5,  64.6, -17.3,  95.2, 144.6,  19.3, 151.2, 227.1,  66.8, 200.6, 306.9, 108.8, 256.5},  //K
		{ 281.4, 397.2, -32.9,  71.2, 182.6,  15.5,   5.0,  85.8,  12.0,  -7.5,  35.5,  -7.6, -30.0,  -5.4, -11.2, -36.7, -38.2,   2.8, -36.0, -67.4,   5.4, -11.8, -83.3,  21.1,  38.0, -96.8,  12.1},  //M
		{ 239.4, 398.5,  18.2,  83.9, 196.1,  41.2,  24.8, 105.6,  13.5, -16.7,  37.4, -24.4, -47.6,  -8.0, -32.9, -51.1, -48.2, -19.3, -25.1, -75.3,   2.7,  21.2, -83.2,  35.2,  97.9, -72.4,  65.4},  //F
		{  68.4, 258.6, -81.0, -79.6,  34.4, -25.1,   5.3,  -9.1,   5.4,  25.5, -15.2,   3.7,   1.4,  -6.7,  14.7,  13.9, -12.7,  33.1,  19.3, -21.3,  43.3,  39.8,  -8.0,  47.3,  47.4, -45.3,  54.1},  //P
		{ 137.8, 219.6, -85.0, -36.9,  76.7, -38.6, -34.6,  21.0,  -1.4,   9.3,  -3.9,   6.8,  14.2, -15.7,  27.3,   2.3, -18.2,  45.6,  -0.9, -29.9,  56.1,  11.6, -46.1,  60.4,  36.5, -73.4,  55.4},  //S
		{ 164.7, 297.9, -69.0,  -4.2, 113.5, -28.7, -25.3,  15.7,   2.1,  -8.4, -23.4,   0.7,   4.7, -27.8,  10.6,  -7.5, -14.7,  30.6,  -7.8, -23.4,  48.5,  -7.0, -36.4,  51.6,  14.8, -60.7,  35.4},  //T
		{ 224.0, 555.6,   9.3,  68.1, 194.7,  34.6,  15.8,  81.4,  12.4, -38.9,  14.5, -24.2, -53.6, -14.3, -34.6, -37.9, -45.2,  -9.5,  -7.1, -63.0,   2.5,  68.0, -82.8,  34.3, 168.0, -54.1,  78.9},  //W
		{ 246.7, 418.8,   7.2,  70.6, 179.0,  33.8,  13.5, 108.7,  11.2, -14.9,  34.8, -26.1, -50.1, -15.2, -34.9, -48.2, -44.9, -13.9, -17.3, -70.6,  15.4,  35.7, -83.2,  39.2, 102.5, -74.6,  79.6},  //Y
		{ 238.3, 466.6,  -6.9,  85.9, 218.4,  37.9,  32.8, 114.0,  36.3,  30.7,  54.0,  -3.8,  -8.8,  14.4, -11.2, -44.3, -25.7, -11.7, -55.0, -63.4, -14.7, -58.9,-105.8, -19.3, -38.3,-137.9, -23.5}}; //V
	// HELIX  SHEET  LOOP  HELIX  SHEET  LOOP    HELIX  SHEET  LOOP    HELIX  SHEET  LOOP   HELIX  SHEET  LOOP  HELIX  SHEET  LOOP    HELIX  SHEET  LOOP  HELIX  SHEET  LOOP   HELIX  SHEET  LOOP
	//          0                    1                    2                    3                    4                    5                    6                    7                    8      -> contact number

	//mutation matrix for 3-class secondary structure types
	//this is a matrix from template to the predicted secondary structure of the sequence
	//it is build on the randomly sampled structure alignments with TMscore 0.34 and 0.83
	/*
	   double SSMutation_MD[3][3]={
	   {0.986896,	-2.67166,	-1.25069}, //HELIX
	   {-2.29363,	1.55292,	-0.530856}, //SHEET
	   {-0.717838,	-0.610854,	0.22509}}; //LOOP
	// HELIX      SHEET      LOOP$
	*/

	double SSMutation_MD[3][3]={
		{ 0.941183, -2.32536,  -0.87487},   //HELIX
		{-2.11462,   1.41307,  -0.401386},  //SHEET
		{-0.760861, -0.540041,  0.269711}}; //LOOP
	// HELIX     SHEET       LOOP
	//
#define SS3_OK		-0.5
#define SS3_GOOD	0
#define SS3_GREAT	0.5

	//mutation matrix for 3-state solvent accessibility
	//these 3 states are divided by equal-frequency method
	//It is build on the randomly sampled structure alignments with TMscore between 0.38 and 0.81
	//this matrix row represents templates and column are for sequence
	/*
	   double ACCMutation_MD[3][3]={
	   {0.812954, 	-0.396657, 	-1.42395}, //buried
	   {-0.0824025, 	0.328262, 	-0.112098}, //medium
	   {-1.33101,	-0.310132,	0.424815}}; //exposed
	// buried        medium           exposed
	*/

	double ACCMutation_MD[3][3]={
		{  0.760885,  -0.0701501,  -0.903965},  //buried
		{ -0.0798508,  0.218512,   -0.0623839}, //medium
		{ -1.14008,   -0.099655,    0.233613}}; //exposed
	//  buried     medium       exposed 

	//8-state secondary structure mutation score
	////from template to the predicted SS8 of the target
	double SS8Mutation_MD[6][6]={
		{ 0.9823,	-0.1944, -1.905,  -0.5508, -1.051,   -1.163},   //H(I)
		{-0.07923, 1.139,	 -0.7431,  0.2849,  0.07965, -0.1479},  //G
		{-1.868,  -0.7317,  1.274,  -0.7433, -0.2456,  -0.07621}, //E(B)
		{-0.4469,  0.2968, -0.8554,  0.9231,  0.2446,  -0.1803},  //T
		{-1.064,   0.0251, -0.3282,  0.3049,  0.6813,   0.2468 }, //S
		{-1.327,  -0.3154, -0.2324, -0.2839,  0.1512,   0.3150}}; //L
	// H(I)     G        E(B)     T        S         L

	//------------- minor calculation ------------//
	//calculate profile score by inner product of PSP vectors
#define PP_OK	-0.20
#define PP_GOOD	0
#define PP_GREAT	0.25


	//profile-profile similarity score
	double MutationOf2Pos4(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"MutationOf2Pos4: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"MutationOf2Pos4: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		//calculate profile score type_i
		double m=0;
		for(int i=0;i<20;i++)
		{
			m+=seq->EmissionProb[sPos][i]*temp->EmissionProb[tPos][i]*pow(2.0,HMMNull[i]/1000.0);
		}
		if(m==0)
		{
			cerr<<"WARNING: possible contaminated data in the HMM emission prob!"<<endl;
			cerr<<"In MutationOfPos4: tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		return log(m);
	}

	//calculate profile score by inner product of PSP vectors, not used any more
	double MutationOf2Pos4_(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"MutationOf2Pos4_: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"MutationOf2Pos4_: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		//calculate profile score type_i
		double m=0;
		for(int i=0;i<20;i++)
		{
			m+=seq->EmissionProb[sPos][i]*temp->EmissionProb[tPos][i]*pow(4.0,HMMNull[i]/1000.0);
		}
		if(m==0)
		{
			cerr<<"WARNING: possible contaminated data in the HMM emission prob!"<<endl;
			cerr<<"In MutationOfPos4_: tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		return log(m);
	}
#define PM_OK	8.0
#define PM_GOOD	9.5
#define PM_GREAT	11

	//calculate profile score by PSP*PSM, frequencey * position specific matrix. seq-template and template-seq
	double MutationOf2Pos5(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"MutationOf2Pos5: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"MutationOf2Pos5: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		//calculate profile score type_ii
		double m=0;
		for(int i=0;i<20;i++)
		{
			m+=seq->EmissionProb[sPos][i]*(temp->EmissionScore[tPos][i]+HMMNull[i])/1000.;
			m+=temp->EmissionProb[tPos][i]*(seq->EmissionScore[sPos][i]+HMMNull[i])/1000.;
		}
		return m;
	}

#define SP_OK	8.0
#define SP_GOOD	9.5
#define SP_GREAT	11

	// similarity score between the primary sequence and the profile: two-way
	double MutationOf2Pos6(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{

		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"MutationOf2Pos6: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"MutationOf2Pos6: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		int tAA=temp->residue[tPos]; //the template residue at tPos
		int sAA=seq->residue[sPos];  //the sequence residue at sPos

		//tAA and sAA are the three letter code orders, so we need to convert them to the alphabetical order
		//to access the HMM information
		int x=AA2SUB[ThreeLetterOrder[tAA]-'A'];
		int y=AA2SUB[ThreeLetterOrder[sAA]-'A'];

		//need to double check to make sure the below method is correct
		//EmissionScore is not substracted by log(background probability), so here we add HMMNull to make it up
		double m=0; 

		//here we silently skip the non-standard amino acids in the sequence or template
		//Maybe we shall report a WARNING???
		if (y>=0 && y<20)
			m+=(temp->EmissionScore[tPos][y]+HMMNull[y])/1000.; // score between the template profile and sequence residue

		if (x<20 && x>=0)
			m+=(seq->EmissionScore[sPos][x]+HMMNull[x])/1000.; //score between the sequence profile and template residue
		return m;
	}
#define SP_ST_OK	4.0
#define SP_ST_GOOD	4.75
#define SP_ST_GREAT	5.5

	// similarity score between the primary sequence and the profile, one-way: primary used for target, profile for template
	double MutationOf2Pos6_ST(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{

		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"MutationOf2Pos6: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"MutationOf2Pos6: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		int tAA=temp->residue[tPos]; //the template residue at tPos
		int sAA=seq->residue[sPos];  //the sequence residue at sPos

		//tAA and sAA are the three letter code orders, so we need to convert them to the alphabetical order
		//to access the HMM information
		int x=AA2SUB[ThreeLetterOrder[tAA]-'A'];
		int y=AA2SUB[ThreeLetterOrder[sAA]-'A'];

		//need to double check to make sure the below method is correct
		//EmissionScore is not substracted by log(background probability), so here we add HMMNull to make it up
		double m=0; 

		//here we silently skip the non-standard amino acids in the sequence or template
		//Maybe we shall report a WARNING???
		if (y>=0 && y<20) m+=(temp->EmissionScore[tPos][y]+HMMNull[y])/1000.; // score between the template profile and sequence residue
		return m;
	}

	// similarity score between the primary sequence and the profile, primary seq for template and profile for sequence
	double MutationOf2Pos6_TS(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{

		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"MutationOf2Pos6: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"MutationOf2Pos6: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		int tAA=temp->residue[tPos]; //the template residue at tPos
		int sAA=seq->residue[sPos];  //the sequence residue at sPos

		//tAA and sAA are the three letter code orders, so we need to convert them to the alphabetical order
		//to access the HMM information
		int x=AA2SUB[ThreeLetterOrder[tAA]-'A'];
		int y=AA2SUB[ThreeLetterOrder[sAA]-'A'];

		//need to double check to make sure the below method is correct
		//EmissionScore is not substracted by log(background probability), so here we add HMMNull to make it up
		double m=0; 

		//here we silently skip the non-standard amino acids in the sequence or template
		if (x<20 && x>=0)
			m+=(seq->EmissionScore[sPos][x]+HMMNull[x])/1000.; //score between the sequence profile and template residue
		return m;
	}

	//calculate secondary structure mutation score
	double SSOf2Pos1(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{

		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"SSOf2Pos1: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"SSOf2Pos1: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		int type = temp->SS[tPos];

		if (type <0 || type > 2){
			cerr<<"SSOf2Pos1: unknown secondary structure type at position "<<tPos<<" in template "<<temp->temp_name<<endl;
			return 0; //shall we terminate the program or just skip it ???
		}

		double score = 0;
		for(int i=0;i<3;i++)
		{
			//sum over the secondary structure mutation score by predicted probability
			score += SSMutation_MD[type][i]*seq->SS2[sPos][i];
		}
		return score;
	}

	//calculate secondary structure match score using a simple edit distance
	double SSOf2Pos2(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{

		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"SSOf2Pos2: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"SSOf2Pos2: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		int type = temp->SS[tPos];
		if(type >2 || type < 0) {
			cerr<<"SSOf2Pos2: unknown secondary structure type at position "<<tPos<<" in template "<<temp->temp_name<<endl;
			return 0; //shall we exit or skip it quietly???
		}

		if(type==seq->SS[sPos])
		{
			if(type==2)return 0;  //do not count loop consistency here?
			else return 1;
		}
		else
		{
			if(type==2 || seq->SS[sPos]==2)return -1;  // if one of them is loop, then penalize by -1
			return -2; //otherwise penalize by -2
		}
	}

	//calculate a 3*3 secondary structure mutation frequency, not used by Jinbo Xu
	void SSOf2Pos3(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq,double out[9])
	{

		if(tPos<0 || sPos<0) 
		{
			cerr<<"In SSOf2Pos3: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		if(tPos>=temp->length || sPos >= seq->length )
		{
			cerr<<"In SSOf2Pos3: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		//calculate 9d map
		for(int i=0;i<9;i++)out[i]=0;

		int temp_type = temp->SS[tPos];
		if(temp_type<0||temp_type>2)
		{
			cerr<<"Unkown template secondary structure type in SSOf2Pos3"<<endl;
			return; //skip it ??? or terminate the program
		}
		//output 9 score
		for(int i=0;i<3;i++)out[temp_type*3+i] = seq->SS2[sPos][i];	
	}

	//mapping the 8-state SS to 6-state indices in the mutation matrix
	int SS82SS6(int type){
		int SS6;
		if (type==0 || type==2) SS6=0;  //map H and I to H(I)
		else if (type==1) SS6=1; //G is still mapped to G
		else if (type==3 || type==4) SS6=2; //E and B to E(B)
		else if (type==5) SS6=3; // T to T
		else if (type==6) SS6=4; // S to S
		else if (type==7) SS6=5; //L to L
		else {
			cerr<<"Unknown 8-state secondary structure type in SS82SS6(int type)!"<<endl;
			exit(-1);
		}

		return SS6;
	}

	double SS8Of2Pos(int tPos, int sPos, TEMPLATE* temp, SEQUENCE* seq)
	{
		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"SSOf2Pos2: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"SSOf2Pos2: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}
		//output score
		double score = 0;
		for(int i=0;i<8;i++)
		{
			if (temp->SS8[tPos][i]<=0.001) continue;  //only one element in SS8[tPos] is equal to 1, others shall be 0
			int tSS6=SS82SS6(i); //mapping

			for(int j=0; j< 8; j++){
				//sum over the secondary structure mutation score by predicted probability
				int sSS6=SS82SS6(j);
				score += SS8Mutation_MD[tSS6][sSS6]*seq->SS8[sPos][j]*temp->SS8[tPos][i];
			}
		}
		return score;
	}


	//calculate the 8-state secondary structure similarity, not used by Jinbo Xu
	double SS8_New_Score1(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
		if(tPos<0 || sPos<0) 
		{
			cerr<<"In SS8_New_Score1: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		if(tPos>=temp->length || sPos >= seq->length )
		{
			cerr<<"In SS8_New_Score1: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		//output probability
		double ss8_tmp=0;
		for(int k=0;k<8;k++)
			ss8_tmp+=temp->SS8[tPos][k]*seq->SS8[sPos][k];
		return ss8_tmp;
	}

	//calculate the 8-state seconadry structure similarity, not used by Jinbo Xu
	void SS8_New_Score3(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq,double out[8])
	{
		if(tPos<0 || sPos<0) 
		{
			cerr<<"In SS8_New_Score3: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		if(tPos>=temp->length || sPos >= seq->length )
		{
			cerr<<"In SS8_New_Score3: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		//calculate diagonal
		for(int i=0;i<8;i++)out[i]=0;
		//output 8 score
		for(int i=0;i<8;i++)out[i]=temp->SS8[tPos][i]*seq->SS8[sPos][i];
	}

	//caclulate solvent accessbility mutation score
	double ACC_New_Score2(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"ACC_New_Score2: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"ACC_New_Score2: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		int tACC=temp->ACC[tPos];
		if (tACC>2 || tACC <0){
			cerr<<"ACC_New_Score2: unknown solvent accessibility status at position "<<tPos<< "of template "<<temp->temp_name<<endl;
			exit(-1);
		}

		int sACC=seq->ACC[sPos];

		if (sACC>2 || sACC<0){
			cerr<<"ACC_New_Score2: Unknown solvent accessibility status at position "<<tPos<< "of sequence "<<seq->seq_name<<endl;
			exit(-1);
		}
		return ACCMutation_MD[tACC][sACC];
	}

	//this function is not used by Jinbo Xu
	void ACC_New_Score3(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq,double out[9])
	{
		if(tPos<0 || sPos<0) 
		{
			cerr<<"In ACC_New_Score3: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		if(tPos>=temp->length || sPos >= seq->length )
		{
			cerr<<"In ACC_New_Score3: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
			exit(-1);
		}
		//calculate 9d map
		for(int i=0;i<9;i++)out[i]=0;

		int temp_type = temp->ACC[tPos];
		if(temp_type<0||temp_type>2)
		{
			cerr<<"In ACC_New_Score3: Unknown template solvent accessibility status at position "<<tPos<<endl;
			exit(-1);
		}
		//output 9 score
		for(int i=0;i<3;i++)out[temp_type*3+i] = seq->acc_our_10_42[sPos][i];	
	}


	//============== old feature ==========//
	/*
	//this score is not used by Jinbo Xu
	double MutationOf2Pos_old(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
	if(tPos<0 || sPos<0) 
	{
	cerr<<"MutationOf2Pos_old: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
	exit(-1);
	}
	if(tPos>=temp->length || sPos >= seq->length )
	{
	cerr<<"MutationOf2Pos_old: out of range by tPos="<<tPos<<" sPos="<<sPos<<endl;
	exit(-1);
	}
	//calculate old mutation
	float* sPSP=seq->PSP[sPos];
	float* tPSM=temp->PSM[tPos];
	double m=0;
	for(int i=0;i<20;i++)
	{
	if (sPSP[i]==0) continue;
	m+=sPSP[i]/PROB_BASE*(tPSM[i]-PSM_aver[i]);
	}
	return m;
	}
	//[contact capacity]
	double ContactCapacityOf2Pos_old(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
	if (tPos<0 || tPos>=temp->length)
	{	
	cerr<<"ContactCapacityOf2Pos_old: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
	exit(-1);
	}
	if (sPos<0 || sPos>=seq->length)
	{
	cerr<<"ContactCapacityOf2Pos_old: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
	exit(-1);
	}

	//judge
	float* sPSP=seq->PSP[sPos];
	int ss=temp->SS[tPos];

	//0, 1 and 2 are for HELIX, SHEET and LOOP, respectively?
	if(ss<0||ss>2){
	cerr<<"ContactCapacityOf2Pos_old: unknown secondary structure type at position "<<tPos<<" of template"<<temp->temp_name<<endl;
	exit(-1);
	}
	int cNum=temp->contactNum_B[tPos];
	if (cNum<0)cNum=0;
	if (cNum>8)cNum=8;
	//calc
	double c=0;
	for(int j=0;j<20;j++)
	{
	if (sPSP[j]==0) continue;
	double contact=Ori_Contact_MD[j][cNum*3+ss];
	c+=sPSP[j]/PROB_BASE*contact;
	}
	return c;
	}
	*/

	//[ws_singleton]
	double SingleOf2Pos_WS(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"ContactCapacityOf2Pos_WS: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"ContactCapacityOf2Pos_WS: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		int ss=temp->SS[tPos];
		int ac=temp->ACC[tPos];
		if(ss<0||ss>2){
			cerr<<"Unknown secondary structure type at position "<<tPos<<" of template "<<temp->temp_name<<endl;
			exit(-1);
		}
		if(ac<0||ac>2){
			cerr<<"Unknown solvent accessibility type at position "<<tPos<<" of template "<<temp->temp_name<<endl;
			exit(-1);
		}
		int j=seq->residue[sPos];
		if (j>=20 || j < 0) return 0;

		return WS_Singleton_MD[j][ss*3+ac];
	}

	//[singleton]
	double SingleOf2Pos_old(int tPos,int sPos,TEMPLATE *temp,SEQUENCE *seq)
	{
		if (tPos<0 || tPos>=temp->length)
		{
			cerr<<"SingleOf2Pos_old: out of range by tPos="<<tPos<<" in template="<<temp->temp_name<<endl;
			exit(-1);
		}
		if (sPos<0 || sPos>=seq->length)
		{
			cerr<<"SingleOf2Pos_old: out of range by sPos="<<sPos<<" in sequence="<<seq->seq_name<<endl;
			exit(-1);
		}

		//judge
		//	float* sPSP=seq->PSP[sPos];
		int ss=temp->SS[tPos];
		int ac=temp->ACC[tPos];
		if(ss<0||ss>2){
			cerr<<"Unknown secondary structure type at position "<<tPos<<" of template "<<temp->temp_name<<endl;
			exit(-1);
		}
		if(ac<0||ac>2){
			cerr<<"Unknown solvent accessibility type at position "<<tPos<<" of template "<<temp->temp_name<<endl;
			exit(-1);
		}
		//calc
		double s=0;
		for(int j=0;j<20;j++)
		{
			s += (seq->EmissionProb[sPos][j]*Ori_Singleton_MD[j][ss*3+ac]);
		}
		return s;
	}


	//========================= calculate singleton features ==============================//
	void CalcFeatures_Singleton(vector<pair<int, int> > &alignment,vector<double> &features, TEMPLATE *temp,SEQUENCE *seq)
	{	
		//features for similarity at the sequence level, use the following features to describe the relationship from high similarity to low
		double seqId=0;  //sequence identity, the number of identical residues in the alignment, maybe normalized by the alignment length or sequence length?
		double blosum80=0; //blosum80 score
		double blosum62=0; // blosum62 score
		double blosum45=0;	//blosum45
		double spScore=0; // mutation score calculated by sequence and profile (i.e., one protein uses sequence and the other profile)
		double spScore_ST=0; // mutation score calculated by primary sequence fo the target and profile of the template
		double pmScore=0; // profile frequency times profile score
		double ppScore=0; //inner product of seq profile normalized by background probability
		double cc=0; //mutation score calculated from structure information
		double hdsm=0; //amino acid mutation score calculated from structure alignment of distantly-related proteins

		//secondary structure similarity,
		double SS3=0;  //3-state
		double SS8=0; //8-state

		//solvent accessibility similarity
		double ACC = 0;

		//contact capacity score, for the relatinship between residue and the number of contacts
		double contactScore=0;	

		//environmental fitness score, for the relationship between residue and structure environement (secondary structure and solvent accessbility)
		double envScore=0;
		double wsEnvScore=0;

		//the following entries only count the positive scores
		double goodB80=0; //only count the positive blosum scores
		double goodB62=0;
		double goodB45=0;
		double goodSP=0;
		double goodSP_ST=0;
		double goodPP=0;
		double goodPM=0;
		double goodCC=0;
		double goodHDSM=0;
		double goodSS3=0;
		double goodSS8=0;
		double goodACC=0;
		double goodContact=0;
		double goodEnv=0;
		double goodWSEnv=0;

		//only count the scores above the threshold XXX_OK
		double okB80=0;
		double okB62=0;
		double okB45=0;
		double okSP=0;
		double okSP_ST=0;
		double okPP=0;
		double okPM=0;
		double okCC=0;
		double okHDSM=0;
		double okSS3=0;

		//only count the scores above the threshold XXX_GREAT
		double gtB80=0;
		double gtB62=0;
		double gtB45=0;
		double gtSP=0;
		double gtSP_ST=0;
		double gtPP=0;
		double gtPM=0;
		double gtCC=0;
		double gtHDSM=0;
		double gtSS3=0;

		int numX = 0; //the number of HIS-tags and missing residues
		int wsend=-1; //the position of the last match state in then alignment
		int wsstart=alignment.size(); //the position of the first match state in the alignment

		//tail gap is not considered
		for(int i=alignment.size()-1;i>=0;i--)
		{
			int x = alignment[i].first; //first is for template 
			int y = alignment[i].second; //second for sequence 
			wsend=i;
			if(x>0 && y>0)	break;
		}

		//also skip the head gap
		for(int i=0; i< alignment.size(); i++)
		{
			int x = alignment[i].first; //first is for template 
			int y = alignment[i].second; //second for sequence 
			wsstart=i;
			if(x>0 && y>0)	break;
		}

		int seqStart=alignment[wsstart].second;
		int seqEnd=alignment[wsend].second;
		int tempStart=alignment[wsstart].first;
		int tempEnd=alignment[wsend].first;

		//lengths
		int tempLength = 0;
		int seqLength = 0;
		int alignLength = 0;

		for(int i=wsstart; i <= wsend; i++)
		{
			int x = alignment[i].first; 	//in the alignment, residue index starts from 1
			int y = alignment[i].second;

			x--,y--; //x and y are residue indices, in the sequence and template class, residue index starts from 0 

			//template length, alignlength and seqlength
			if(alignment[i].first>0 && alignment[i].second>0) alignLength++;
			if(alignment[i].first>0) tempLength++;
			if(alignment[i].second>0) seqLength++;

			if (alignment[i].first<=0 || alignment[i].second<=0) continue;

			// non-standard amino acids
			if(temp->residue[x] >=20 || seq->residue[y] >=20){ 
				numX++; continue;
			}

			//his-tag
			if(temp->isMultiHIS[x] || seq->isMultiHIS[y]){
				numX++; continue;
			}

			//missing residues
			if(temp->isMissing[x]){
				numX++; continue;
			}

			//tAA and sAA are the three letter code orders of amino acids
			int tAA=temp->residue[x]; //template residue type
			int sAA=seq->residue[y]; //sequence residue type
			if (tAA <0){
				cerr<<"Unknown template residue type at position "<< x <<" of template "<<temp->temp_name<<endl;
				numX++; continue;
				//			exit(-1);
			}
			if (sAA<0){
				cerr<<"Unknown sequence residue type at position "<< y <<" of target "<<seq->seq_name<<endl;
				numX++; continue;
				//			exit(-1);
			}

			//seq identity
			seqId+=(tAA==sAA);

			//need to double check to make sure the residue[x] and residue[y] represents the amino acid id
			char tOneLetterCode=ThreeLetterOrder[tAA];
			char sOneLetterCode=ThreeLetterOrder[sAA];
			blosum80+=BLOSUM80_Calc(tOneLetterCode, sOneLetterCode);
			blosum62+=BLOSUM62_Calc(tOneLetterCode, sOneLetterCode);
			blosum45+=BLOSUM45_Calc(tOneLetterCode, sOneLetterCode);
			spScore+=MutationOf2Pos6(x, y, temp, seq);
			spScore_ST+=MutationOf2Pos6_ST(x, y, temp, seq);
			pmScore+=MutationOf2Pos5(x, y, temp, seq);
			ppScore+=MutationOf2Pos4(x, y, temp, seq);
			cc+=(Ori_CC50_MD[tAA][sAA]-0.5); 	//0.5 is the expected value?
			hdsm+=Ori_HDSM_MD[tAA][sAA];

			//secondary structure score
			SS3+=SSOf2Pos1(x,y,temp,seq);
			SS8+=SS8Of2Pos(x,y,temp,seq);
			//acc
			ACC+=ACC_New_Score2(x,y,temp,seq);
			//singleton score
			envScore+=SingleOf2Pos_old(x,y,temp,seq);
			wsEnvScore+=SingleOf2Pos_WS(x,y,temp,seq);


			//-> goodScore
			goodB80+=MAX(0, BLOSUM80_Calc(tOneLetterCode, sOneLetterCode));
			goodB62+=MAX(0, BLOSUM62_Calc(tOneLetterCode, sOneLetterCode));
			goodB45+=MAX(0, BLOSUM45_Calc(tOneLetterCode, sOneLetterCode));
			goodSP+=GREAT(SP_GOOD,MutationOf2Pos6(x, y, temp, seq));
			goodSP_ST+=GREAT(SP_ST_GOOD,MutationOf2Pos6_ST(x, y, temp, seq));
			goodPP+=GREAT(PP_GOOD,MutationOf2Pos4(x, y, temp, seq));
			goodPM+=GREAT(PM_GOOD,MutationOf2Pos5(x, y, temp, seq));
			goodCC+=MAX(0,(Ori_CC50_MD[tAA][sAA]-0.5)); 
			goodHDSM+=MAX(0,Ori_HDSM_MD[tAA][sAA]);
			goodSS3+=MAX(0, SSOf2Pos1(x,y,temp,seq));
			goodSS8+=MAX(0, SS8Of2Pos(x,y,temp,seq));
			goodACC+=MAX(0, ACC_New_Score2(x,y,temp,seq));
			goodEnv+=MIN(0,SingleOf2Pos_old(x,y,temp,seq)); //the smaller, the better
			goodWSEnv+=MAX(0,SingleOf2Pos_WS(x,y,temp,seq)); //the higher, the better


			gtB80+=GREAT(B80_GREAT, BLOSUM80_Calc(tOneLetterCode, sOneLetterCode));
			gtB62+=GREAT(B62_GREAT, BLOSUM62_Calc(tOneLetterCode, sOneLetterCode));
			gtB45+=GREAT(B45_GREAT, BLOSUM45_Calc(tOneLetterCode, sOneLetterCode));
			gtSP+=GREAT(SP_GREAT,MutationOf2Pos6(x, y, temp, seq));
			gtSP_ST+=GREAT(SP_ST_GREAT,MutationOf2Pos6_ST(x, y, temp, seq));
			gtPP+=GREAT(PP_GREAT,MutationOf2Pos4(x, y, temp, seq));
			gtPM+=GREAT(PM_GREAT,MutationOf2Pos5(x, y, temp, seq));
			gtCC+=GREAT(CC_GREAT,(Ori_CC50_MD[tAA][sAA]-0.5)); 
			gtHDSM+=GREAT(HDSM_GREAT,Ori_HDSM_MD[tAA][sAA]);
			gtSS3+=GREAT(SS3_GREAT, SSOf2Pos1(x,y,temp,seq));

			okB80+=GREAT(B80_OK, BLOSUM80_Calc(tOneLetterCode, sOneLetterCode));
			okB62+=GREAT(B62_OK, BLOSUM62_Calc(tOneLetterCode, sOneLetterCode));
			okB45+=GREAT(B45_OK, BLOSUM45_Calc(tOneLetterCode, sOneLetterCode));
			okSP+=GREAT(SP_OK,MutationOf2Pos6(x, y, temp, seq));
			okSP_ST+=GREAT(SP_ST_OK,MutationOf2Pos6_ST(x, y, temp, seq));
			okPP+=GREAT(PP_OK,MutationOf2Pos4(x, y, temp, seq));
			okPM+=GREAT(PM_OK,MutationOf2Pos5(x, y, temp, seq));
			okCC+=GREAT(CC_OK,(Ori_CC50_MD[tAA][sAA]-0.5)); 
			okHDSM+=GREAT(HDSM_OK,Ori_HDSM_MD[tAA][sAA]);
			okSS3+=GREAT(SS3_OK, SSOf2Pos1(x,y,temp,seq));
		}

		alignLength-=numX;
		seqLength-=numX;
		tempLength-=numX;
		if(seqLength<=0)seqLength=1;
		if(tempLength<=0)tempLength=1;
		if(alignLength<=0)alignLength=1;

		features.clear();

		features.push_back(-8); // separator, 8 is the index of the next field

		//features for lengths, [1-4]
		features.push_back(seq->length); // the length of the whole sequence
		features.push_back(seqLength);  // the number of residues in sequence, excluding head and tail insertion, his-tags and non-standard amino acids
		features.push_back(tempLength); // the number of residues in template, excluding head and tail insertion, the missing residues, his-tags and non-standard AAs
		features.push_back(temp->length-numX); 	// the number of residues in the template, excluding missing residues, his-tag and non-standard

		//start and end positions of the alignment, not real features
		features.push_back(seqStart);
		features.push_back(seqEnd);
		features.push_back(tempStart);
		features.push_back(tempEnd);

		features.push_back(-17); // separator, 17 is the index of the next field in the feature file

		//features for sequence similarity, [5-14]
		features.push_back(alignLength); 	// the number of the match states
		features.push_back(seqId); 		//
		features.push_back(blosum80);		//
		features.push_back(blosum62);		//
		features.push_back(blosum45);		//
		features.push_back(spScore);		//
		features.push_back(ppScore);		//
		features.push_back(pmScore);		//
		features.push_back(cc);			//
		features.push_back(hdsm);		//

		features.push_back(-28); //separator

		//features for structure similarity, [14-18]
		features.push_back(SS3);   		//
		features.push_back(SS8);   		//
		features.push_back(ACC);  		//
		features.push_back(envScore);		//
		features.push_back(contactScore);	//



		//-------------------------- the following procedures are for gap features --------------------//
		vector<int> templateGaps; //record the insertion gap lengths for the template, head and tail gaps are not counted
		vector<int> sequenceGaps; //record the insertion gap lengths for the sequence, head and tail gaps are not counted

		//features for gaps, use the below five basic functions to fit the relationship between the alignment quality and the gaps
		double s_gapOpens=0; //the number of gap openings, can be interpreted as \sum_g g^0 where g is the length of one gap
		double s_gapExt=0; //gap extesion, can be interpreted as \sum_g g 
		double s_gapRoot=0; // \sum_g g^{1/2} where g is the gap length
		double s_gapLog=0; // \sum_g log2(g)
		double s_gapRootLog=0; // \sum_g glog2(g)

		double t_gapOpens=0; //the number of gap openings, can be interpreted as \sum_g g^0 where g is the length of one gap
		double t_gapExt=0; //gap extesion, can be interpreted as \sum_g g 
		double t_gapRoot=0; // \sum_g g^{1/2} where g is the gap length
		double t_gapLog=0; // \sum_g log2(g)
		double t_gapRootLog=0; // \sum_g glog2(g)

		double s_realGapExt=0;
		double t_realGapExt=0;

		//detect sequence insertion gaps, a sequence position is an insertion if its corresponding template positipon >0
		int gapLength=0;
		for(int i=wsstart; i<=wsend; i++){ //head and tail gaps are not counted
			int x=alignment[i].first;
			int y=alignment[i].second;
			if (x<=0 && y<=0){
				cerr<<"Unexpected alignment status at alignment position "<<i<<endl;
				cerr<<"Both sequence and template aligned positions are <=0"<<endl;
				exit(-1);
			}
			if (x<0){
				//continue with the gap
				gapLength++;
			}else if (x>0){
				if (gapLength>0) {
					//end current gap
					sequenceGaps.push_back(gapLength);
					gapLength=0; //reset gap length for the next gap
				}
			}else{
				cerr<<"Unexpected alignment status at alignment position "<<i<<endl;
				cerr<<"One aligned template position is 0"<<endl;
				exit(-1);
			}		
		}

		//detect template insertion gaps
		gapLength=0;
		for(int i=wsstart; i<wsend; i++){ //head and tail gaps are not counted
			int x=alignment[i].first;
			int y=alignment[i].second; //y is the index of sequence residue
			if (x<=0 && y<=0){
				cerr<<"Unexpected alignment status at alignment position "<<i<<endl;
				cerr<<"Both sequence and template aligned positions are <=0"<<endl;
				exit(-1);
			}
			if (y<0) {
				gapLength++;
			}else if (y>0){
				if (gapLength>0){
					//end current gap
					templateGaps.push_back(gapLength);
					gapLength=0; //reset gapLength for the next gap
				}
			}else{
				cerr<<"Unexpected status at alignment position "<<i<<endl;
				cerr<<"One aligned sequence position is 0"<<endl;
				exit(-1);
			}
		}

		//calculate features for gaps
		t_gapOpens=templateGaps.size();
		s_gapOpens=sequenceGaps.size();

		//-> template gap
		for(int i=0; i<templateGaps.size(); i++){
			if (templateGaps[i]<1){
				cerr<<"One template gap length is 0 in the alignment between template "<<temp->temp_name<<" and sequence "<<seq->seq_name<<endl;
				exit(-1);
			}
			t_realGapExt+=templateGaps[i];

			//reduce the size of a long gap to base + (g-base)^{2/3}
#define SHORTEN(g)	( (g)>LONG_GAP?(LONG_GAP+SQ(cbrt((g)-LONG_GAP))):(g))

			double effGapLen=SHORTEN(templateGaps[i]);

			t_gapExt+=effGapLen;
			t_gapRoot+=sqrt(effGapLen);
			t_gapRootLog+=sqrt(effGapLen)*log2(effGapLen);
			t_gapLog+=log2(effGapLen);
		}

		//-> sequence gap
		for(int i=0; i<sequenceGaps.size(); i++){
			if (sequenceGaps[i]<1){
				cerr<<"One sequence gap length is 0 in the alignment between template "<<temp->temp_name<<" and sequence "<<seq->seq_name<<endl;
				exit(-1);
			}
			s_realGapExt+=sequenceGaps[i];

			double effGapLen=SHORTEN(sequenceGaps[i]);

			s_gapExt+=effGapLen;
			s_gapRoot+=sqrt(effGapLen);
			s_gapRootLog+=sqrt(effGapLen)*log2(effGapLen);
			s_gapLog+=log2(effGapLen);
		}

		//--------- push_back gap features ---------//

		//this feature is only used as separator	
		features.push_back(-34); //separator
		//features for gap
		features.push_back(s_gapOpens);		//
		features.push_back(s_gapExt);		
		features.push_back(s_gapRoot);		//
		features.push_back(s_gapRootLog);	//
		features.push_back(s_gapLog);		//
		features.push_back(t_gapOpens);		//
		features.push_back(t_gapExt);		//
		features.push_back(t_gapRoot);		//
		features.push_back(t_gapRootLog);	//
		features.push_back(t_gapLog);		//

		features.push_back(-45); //45 indicates the index of the next field
		//features for good similairty scores
		features.push_back(goodB80);
		features.push_back(goodB62);
		features.push_back(goodB45);
		features.push_back(goodSP);
		features.push_back(goodPP);
		features.push_back(goodPM);
		features.push_back(goodCC);
		features.push_back(goodHDSM);
		features.push_back(goodSS3);
		features.push_back(goodSS8);
		features.push_back(goodACC);
		features.push_back(goodEnv);
		features.push_back(goodContact);

		features.push_back(-59); //separator
		//features for great similarity scores
		features.push_back(gtB80);
		features.push_back(gtB62);
		features.push_back(gtB45);
		features.push_back(gtSP);
		features.push_back(gtPP);
		features.push_back(gtPM);
		features.push_back(gtCC);
		features.push_back(gtHDSM);
		features.push_back(gtSS3);

		features.push_back(-69); //separator
		//features for ok similairty score
		features.push_back(okB80);
		features.push_back(okB62);
		features.push_back(okB45);
		features.push_back(okSP);
		features.push_back(okPP);
		features.push_back(okPM);
		features.push_back(okCC);
		features.push_back(okHDSM);
		features.push_back(okSS3);

		features.push_back(-79); //separator, 79 indicates the index of the next field
		//features for wang sheng's environmental fitness score
		features.push_back(wsEnvScore);
		features.push_back(goodWSEnv);

		features.push_back(-82); //separator, 82 is the index of the next field
		//the real gap extensions, may not be used as features
		features.push_back(s_realGapExt);	//
		features.push_back(t_realGapExt);	//

		features.push_back(-85); //separator, 85 is the index of the next field
		//SP_ST scores, calculate from the primary sequence of the target and profile of the template
		features.push_back(spScore_ST);
		features.push_back(okSP_ST);
		features.push_back(goodSP_ST);
		features.push_back(gtSP_ST);

	}



	//========================= calculate pairwise features ==============================//
	const double d1[12][12]={
		{0.0622186,0.0705997,0.03042,0.0165638,0.0117084,0.00614276,0.00438287,0.00361513,0.00289583,0.002541,0.00228398,0.00163402},
		{0.0353955,0.294967,0.117701,0.0467576,0.0285148,0.0144187,0.011343,0.0082282,0.00588411,0.00481453,0.0038883,0.00119637},
		{0.0192365,0.141213,0.239797,0.0926083,0.0403234,0.0202416,0.0147372,0.011901,0.00831808,0.00640026,0.00495212,0.00153353},
		{0.0128414,0.0599759,0.10212,0.161021,0.0743214,0.0317909,0.0203743,0.0159278,0.0113996,0.00881565,0.00649663,0.00195537},
		{0.0167513,0.0615853,0.0747943,0.127638,0.185411,0.0813447,0.0430069,0.0321503,0.022596,0.017589,0.0132976,0.00356997},
		{0.0173854,0.048069,0.0594827,0.0882908,0.133093,0.21426,0.12213,0.0645857,0.0425284,0.031167,0.0225033,0.00573965},
		{0.0196927,0.0508668,0.0581221,0.0757459,0.0941273,0.16669,0.236824,0.127083,0.0707162,0.0495832,0.0346026,0.00810546},
		{0.0185753,0.0346239,0.0450074,0.0576078,0.0680623,0.0853465,0.124869,0.162596,0.111991,0.0678178,0.0454404,0.00988924},
		{0.019857,0.0270043,0.0353276,0.0463251,0.0532623,0.0630318,0.0782179,0.12757,0.177967,0.111814,0.0678046,0.0127051},
		{0.0215785,0.0248274,0.0295837,0.0388548,0.0451474,0.0500268,0.0600359,0.0847039,0.123292,0.157213,0.113726,0.0169458},
		{0.0241316,0.0225234,0.0255679,0.0323806,0.0384542,0.0407529,0.0470875,0.0636283,0.0840021,0.128753,0.180507,0.0231893},
		{0.732336,0.163744,0.182077,0.216206,0.227575,0.225954,0.236992,0.298012,0.33841,0.413492,0.504498,0.913536}
	};

	const double d_ref[12] = {0.0014019,0.00364422,0.00435592,0.00494332,0.00861316,0.0152909,0.0213476,0.0197301,0.0229548,0.0245999,0.0282358,0.844883};
	const double pot_ref[12] = {0.437798044,0.088590826,0.142051864,0.178072077,0.198164495,0.119791129,0.110496167,0.218379552,0.219465592,0.269886421,0.266215469,-0.115108831};

	//----- map_distrance2index ---//
	int map_distrance2index(double dis)
	{
		int abs_d = (int)dis;
		if(abs_d<4)return 0;
		if(abs_d>=15)return 11;
		return abs_d-4;
	}

	//------ Compute_Pairwise_Poteintial2 ------//
	//-> note that, TEMPLATE *temp should conduct distance calculation -> "temp->Compute_All_Distance();"
	//-> while SEQUENCE *seq should load epad_prob file -> "seq->Read_distance_file(filename);"
	void CalcFeatures_Pairwise(vector<pair<int,int> > & alignment,vector<double> &features, TEMPLATE *temp,SEQUENCE *seq)
	{
		double pairwise_score = 0;
		vector <double> pairwise_part_score (12,0);

		features.clear();
		for(int i=0;i<alignment.size();i++)
		{
			pair<int,int> i_pos = alignment[i];
			if(i_pos.first>0 && i_pos.second>0 && temp->isMissing[i_pos.first-1]==0)
			{
				double sum_pair_pot = 0;
				vector <double> sum_pair_part_pot (12,0);
				int sum_pair_num = 0;
				for(int j=0;j<alignment.size();j++)
				{
					pair<int,int> j_pos = alignment[j];
					if(j_pos.first>0 && j_pos.second>0 && abs(i_pos.second-j_pos.second)>=6 && temp->isMissing[j_pos.first-1]==0 && i!=j )
					{
						double t_d = temp->dis_matrix[i_pos.first-1][j_pos.first-1];
						int t_index = map_distrance2index(t_d);
						double prob_d = 0;
						double prob_part_d;
						for(int n=0;n<12;n++)
						{
							prob_d += d1[t_index][n] * seq->pair_dis[i_pos.second-1][j_pos.second-1][n];
							prob_part_d = d1[t_index][n] * seq->pair_dis[i_pos.second-1][j_pos.second-1][n];
							sum_pair_part_pot[n] += log( prob_part_d / d_ref[t_index] )- pot_ref[t_index];
						}
						double pot = log( prob_d / d_ref[t_index]) - pot_ref[t_index];
						sum_pair_pot += pot;
						sum_pair_num ++;
					}
				}
				if(sum_pair_num!=0)
				{
					pairwise_score += sum_pair_pot/(1.0*sum_pair_num);
					for(int n=0;n<12;n++)pairwise_part_score[n] += sum_pair_part_pot[n]/(1.0*sum_pair_num);
				}
			}
		}

		//---- push_back features -----//
		features.push_back(pairwise_score);
		for(int n=0;n<12;n++)features.push_back(pairwise_part_score[n]);

		//	return pairwise_score;
	}


	//===================== CalcFeatures ==================//
	void CalcFeatures(vector<pair<int,int> > & alignment,vector<double> &features, TEMPLATE *temp,SEQUENCE *seq)
	{
		features.clear();
		vector <double> features_singleton;
		vector <double> features_pairwise;
		CalcFeatures_Singleton(alignment,features_singleton,temp,seq);
		CalcFeatures_Pairwise(alignment,features_pairwise,temp,seq);
		//final
		features.insert( features.end(), features_singleton.begin(), features_singleton.end() );
		features.insert( features.end(), features_pairwise.begin(), features_pairwise.end() );
	}
}

/*
//---------- main ---------//
int main(int argc,char **argv)
{
	//---- PDB Merge ----//
	{
		if(argc<5)
		{
			printf("CalcFeatures <tpl_file> <tgt_file> <epad_file> <align_file> \n");
			exit(-1);
		}
		string tpl_file=argv[1];
		string tgt_file=argv[2];
		string epad_file=argv[3];
		string fasta_file=argv[4];

		//load
		//-> tpl
		string tpl_name,tpl_root;
		CalcFeature::getBaseName(tpl_file,tpl_name,'/','.');
		CalcFeature::getRootName(tpl_file,tpl_root,'/');
		TEMPLATE* t;
		t = new TEMPLATE(tpl_name,tpl_root,1);
		//-> tgt
		string tgt_name,tgt_root;
		CalcFeature::getBaseName(tgt_file,tgt_name,'/','.');
		CalcFeature::getRootName(tgt_file,tgt_root,'/');
		SEQUENCE* s;
		s = new SEQUENCE(tgt_name,tgt_root,1);
		//only need to calculate it once
		//-> distance related
		t->Compute_All_Distance();
		s->Read_distance_file(epad_file);
		//-> alignment
		vector<pair<int, int> > alignment;
		string nam1_content,nam2_content;
		string nam1_full,nam2_full;
		string nam1,nam2;
		int retv=CalcFeature::ReadToFile_FASTA(fasta_file,alignment,nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2);

		//calc
		vector<double> features;
		CalcFeature::CalcFeatures(alignment,features,t,s);
		for(int i=0;i<(int)features.size();i++)printf("%lf ",features[i]);
		printf("\n");
		//exit
		exit(0);
	}
*/

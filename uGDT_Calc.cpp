#include "uGDT_Calc.h"
namespace uGDT_Calc {
	//======================= I/O related ==========================//
	//-------- utility ------//
	/*
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
	   */

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

	//---------- dynamic programming ----------//
	int WWW_Advance_Align_Dyna_Prog_Double(int n1,int n2,const vector<double> &score,
			double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
			double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
			vector<pair<int,int> > & alignment,double &ali_sco)
	{
		int i,j;
		//input
		int m = n1 + 1;  // +1 to account for the extra row,col in
		int n = n2 + 1;  // the DP matrices corresponding to gaps
		int DP_maximal=n;
		int IN_maximal=n2;
		//const value
		const int _H_  = 0;
		const int _S_  = 1;
		const int _V_  = 2;

		//create D and M
		vector <int> D[3];      // the path (directions) matrix
		vector <double> M[3];   // the current scores (values) matrix
		//resize(m,n)
		for (i = 0; i < 3; ++i) 
		{
			D[i].resize(m*n);
			M[i].resize(m*n);
		}
		//init()
		double WS_MIN=-1000000;
		D[_S_][0*DP_maximal+ 0] = -1;
		D[_H_][0*DP_maximal+ 0] = -1;
		D[_V_][0*DP_maximal+ 0] = -1;
		M[_S_][0*DP_maximal+ 0] = 0;
		M[_H_][0*DP_maximal+ 0] = WS_MIN;
		M[_V_][0*DP_maximal+ 0] = WS_MIN;
		for (i = 1; i < m; i++) 
		{
			D[_S_][i*DP_maximal+ 0] = _V_;
			D[_H_][i*DP_maximal+ 0] = _V_;
			D[_V_][i*DP_maximal+ 0] = _V_;
			M[_S_][i*DP_maximal+ 0] = WS_MIN;
			M[_H_][i*DP_maximal+ 0] = WS_MIN;
			M[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
		}
		for (j = 1; j < n; j++) 
		{
			D[_S_][0*DP_maximal+ j] = _H_;
			D[_H_][0*DP_maximal+ j] = _H_;
			D[_V_][0*DP_maximal+ j] = _H_;
			M[_S_][0*DP_maximal+ j] = WS_MIN;
			M[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
			M[_V_][0*DP_maximal+ j] = WS_MIN;
		}
		//fill(firstSeq, secondSeq, distFunc);
		double gap_open;
		double gap_ext;
		double v1,v2,v3;
		double dist;
		for (i = 1; i < m; i++) 
		{
			for (j = 1; j < n; j++) 
			{
				//condition upper
				if(j==n-1)
				{
					gap_open=GAP_TAIL1;
					gap_ext=GAP_TAIL1;
				}
				else
				{
					gap_open=GAP_OPEN1;
					gap_ext=GAP_EXT1;
				}
				v1 = M[_V_][(i-1)*DP_maximal+ j] + gap_ext;
				v2 = M[_S_][(i-1)*DP_maximal+ j] + gap_open;
				v3 = M[_H_][(i-1)*DP_maximal+ j] + gap_open;
				M[_V_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
				if (M[_V_][i*DP_maximal+ j] == v1) D[_V_][i*DP_maximal+ j] = _V_;
				else if(M[_V_][i*DP_maximal+ j] == v2) D[_V_][i*DP_maximal+ j] = _S_;
				else D[_V_][i*DP_maximal+ j] = _H_;
				//condition left
				if(i==m-1)
				{
					gap_open=GAP_TAIL2;
					gap_ext=GAP_TAIL2;
				}
				else
				{
					gap_open=GAP_OPEN2;
					gap_ext=GAP_EXT2;
				}
				v1 = M[_H_][i*DP_maximal+ j-1] + gap_ext;
				v2 = M[_S_][i*DP_maximal+ j-1] + gap_open;
				v3 = M[_V_][i*DP_maximal+ j-1] + gap_open;
				M[_H_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
				if (M[_H_][i*DP_maximal+ j] == v1) D[_H_][i*DP_maximal+ j] = _H_;
				else if(M[_H_][i*DP_maximal+ j] == v2) D[_H_][i*DP_maximal+ j] = _S_;
				else D[_H_][i*DP_maximal+ j] = _V_;
				//condition diag
				dist = score.at((i-1)*IN_maximal+ j-1);  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
				v1 = M[_V_][(i-1)*DP_maximal+ j-1] + dist;
				v2 = M[_H_][(i-1)*DP_maximal+ j-1] + dist;
				v3 = M[_S_][(i-1)*DP_maximal+ j-1] + dist;
				M[_S_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
				if (M[_S_][i*DP_maximal+ j] == v3) D[_S_][i*DP_maximal+ j] = _S_;
				else if (M[_S_][i*DP_maximal+ j] == v1) D[_S_][i*DP_maximal+ j] = _V_;
				else D[_S_][i*DP_maximal+ j] = _H_;
			}
		}
		//build(ali, firstSeq, secondSeq, distFunc);
		i = m-1;
		j = n-1;
		v1=M[_V_][i*DP_maximal+ j];
		v2=M[_H_][i*DP_maximal+ j];
		v3=M[_S_][i*DP_maximal+ j];
		double maximal = std::max(v1, std::max(v2, v3));
		int k = -1;
		if(v3==maximal)k = _S_;
		else if(v2==maximal)k = _H_;
		else k = _V_;
		//trace_back
		alignment.clear();
		int count = 0;
		int matches = 0;
		int cur_case=k;
		int pre_case;
		for(;;)
		{
			if(i==0||j==0)break;
			pre_case=D[cur_case][i*DP_maximal+ j];
			switch (cur_case)
			{
				case _S_:
					alignment.push_back(pair<int,int>(i,j)); 
					i--;
					j--;
					++matches;
					break;
				case _V_:
					alignment.push_back(pair<int,int>(i,-j)); 
					i--;
					break;
				case _H_:
					alignment.push_back(pair<int,int>(-i,j)); 
					j--;
					break;
				default:
					cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
						<< D[k][i*DP_maximal+ j] << endl;
					exit(-1);
			}
			cur_case=pre_case;
			count++;
		}
		while (j> 0) alignment.push_back(pair<int,int>(-i,j)),j--;
		while (i> 0) alignment.push_back(pair<int,int>(i,0)), i--;
		reverse(alignment.begin(), alignment.end());
		ali_sco=maximal;
		return matches;
	}

	//============= Exttract Alignment ========//
	//-> get mapping alignment detail
	void process_oriami_record_simp(const char *seq_,const char *ami_,
			vector<pair<int,int> > &WWW_alignment)
	{
		int i,j;
		int n1,n2;
		//--[1]dynamic_programming	
		n1=(int)strlen(seq_);
		n2=(int)strlen(ami_);
		vector <double> WWW_score;
		WWW_score.resize(n1*n2);
		for(i=0;i<n1;i++)
		{
			for(j=0;j<n2;j++)
			{
				if(seq_[i]==ami_[j])
				{
					if(seq_[i]=='X'||seq_[i]=='Z'||seq_[i]=='.')WWW_score[i*n2+j]=0;
					else if(ami_[j]=='X'||ami_[j]=='Z'||ami_[j]=='.')WWW_score[i*n2+j]=0;
					else WWW_score[i*n2+j]=10;
				}
				else
				{
					if(seq_[i]=='X'||seq_[i]=='Z'||seq_[i]=='.')WWW_score[i*n2+j]=0;
					else if(ami_[j]=='X'||ami_[j]=='Z'||ami_[j]=='.')WWW_score[i*n2+j]=0;
					else WWW_score[i*n2+j]=-15;
				}
			}
		}
		double sco;
		int matchs;
		matchs=WWW_Advance_Align_Dyna_Prog_Double(n1,n2,WWW_score,-11,-1,-11,-1,0,0,0,0,
				WWW_alignment,sco);
	}

	//-> extract alignment
	//[note]: we regard namX_pdb as alignment_out
	//        the namX_seq as alignment_in
	//[example]:
	/*
	//         AAAAABBBBBCCCCDDDDD  -> 1 (nam1_pdb) \ -> ali1 \
	//              BBBBBCCCC       -> 2 (nam1_seq)  \         \
	//               ||||||                          -> ali2   -> ali4
	//               BBBBCC         -> 3 (nam2_seq)  /         /
	//          AAAABBBBBCCCCDD     -> 4 (nam2_pdb) / -> ali3 /
	*/
	//[note]: ali2 is from alignment_in, and ali4 is from alignment_out
	void Extract_Alignment(string &nam1_seq,string &nam1_pdb,string &nam2_seq,string &nam2_pdb,
			vector<pair<int,int> > &alignment_in,vector<pair<int,int> > &alignment_out,
			vector <int> &ali1_seq_pdb, vector <int> &ali2_seq_pdb)
	{
		//init
		int i;
		int l1,l2,l3,l4;
		l1=(int)nam1_pdb.length();
		l2=(int)nam1_seq.length();
		l3=(int)nam2_seq.length();
		l4=(int)nam2_pdb.length();
		int *ali1=new int[l1];
		int *ali2=new int[l2];
		int *ali3=new int[l3];
		int *ali4=new int[l1];
		vector<pair<int,int> > WWW_alignment;
		//first alignment
		ali1_seq_pdb.resize(l2);
		process_oriami_record_simp(nam1_pdb.c_str(),nam1_seq.c_str(),WWW_alignment);
		for(i=0;i<l1;i++)ali1[i]=-1;
		for(i=0;i<l2;i++)ali1_seq_pdb[i]=-1;
		for(i=0;i<(int)WWW_alignment.size();i++)
		{
			int pos1=WWW_alignment[i].first;
			int pos2=WWW_alignment[i].second;
			if(pos1>0 && pos2>0)ali1[pos1-1]=pos2-1;
			if(pos1>0 && pos2>0)ali1_seq_pdb[pos2-1]=pos1-1;
		}
		//second alignment
		for(i=0;i<l2;i++)ali2[i]=-1;
		for(i=0;i<(int)alignment_in.size();i++)
		{
			int pos1=alignment_in[i].first;
			int pos2=alignment_in[i].second;
			if(pos1>0 && pos2>0)ali2[pos1-1]=pos2-1;
		}
		//third alignment
		ali2_seq_pdb.resize(l3);
		process_oriami_record_simp(nam2_seq.c_str(),nam2_pdb.c_str(),WWW_alignment);
		for(i=0;i<l3;i++)ali3[i]=-1;
		for(i=0;i<l3;i++)ali2_seq_pdb[i]=-1;
		for(i=0;i<(int)WWW_alignment.size();i++)
		{
			int pos1=WWW_alignment[i].first;
			int pos2=WWW_alignment[i].second;
			if(pos1>0 && pos2>0)ali3[pos1-1]=pos2-1;
			if(pos1>0 && pos2>0)ali2_seq_pdb[pos1-1]=pos2-1;
		}
		//fourth alignment
		for(i=0;i<l1;i++)ali4[i]=-1;
		for(i=0;i<l1;i++)
		{
			int pos1=ali1[i];
			if(pos1==-1)continue;
			int pos2=ali2[pos1];
			if(pos2==-1)continue;
			int pos3=ali3[pos2];
			if(pos3==-1)continue;
			ali4[i]=pos3;
		}
		//final
		{
			alignment_out.clear();
			//start
			int j;
			int ii,jj;
			int wlen;
			int pre_ii=0;
			int pre_jj=0;
			for(i=1;i<=l1;i++)
			{
				ii=i;
				jj=ali4[i-1];  //ali1 starts from 0, correspondence also from 0
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
						alignment_out.push_back (pair<int,int>(pre_ii, -pre_jj)); //Ix
					}
					wlen=jj-pre_jj;
					for(j=1;j<wlen;j++)
					{
						pre_jj++;
						alignment_out.push_back (pair<int,int>(-pre_ii, pre_jj)); //Iy
					}
					//current_path
					alignment_out.push_back (pair<int,int>(ii, jj)); //Match
					//update
					pre_ii=ii;
					pre_jj=jj;
				}
			}
			//termi
			pre_ii++;
			for(i=pre_ii;i<=l1;i++)alignment_out.push_back (pair<int,int>(i, -pre_jj)); //Ix
			pre_jj++;
			for(i=pre_jj;i<=l4;i++)alignment_out.push_back (pair<int,int>(-l1, i));  //Iy
		}
		//delete
		delete [] ali1;
		delete [] ali2;
		delete [] ali3;
		delete [] ali4;
	}


	//========== process PDB ============//
	char WWW_Three2One_III(const char *input)
	{
		int i;
		int len;
		int result;
		//encoding
		len=(int)strlen(input);
		if(len!=3)return 'X';
		result=0;
		for(i=0;i<len;i++)result+=(input[i]-'A')*(int)pow(1.0*26,1.0*i);
		//switch
		switch(result)
		{
			case 286:return 'A';
			case 4498:return 'R';
			case 9256:return 'N';
			case 10608:return 'D';
			case 12794:return 'C';
			case 9080:return 'Q';
			case 13812:return 'E';
			case 16516:return 'G';
			case 12383:return 'H';
			case 2998:return 'I';
			case 13635:return 'L';
			case 12803:return 'K';
			case 12960:return 'M';
			case 2901:return 'F';
			case 9921:return 'P';
			case 11614:return 'S';
			case 11693:return 'T';
			case 10601:return 'W';
			case 12135:return 'Y';
			case 7457:return 'V';
			default:return 'X';
		}
	}

	//--------- PDB_To_XYZ ----------//
	int PDB_To_XYZ_NoChain(string &pdb,char *ami,XYZ *xyz)
	{
		//--- list for mapping ---//
		map<string, int > ws_mapping;
		map<string, int>::iterator iter;
		ws_mapping.clear();
		ifstream fin;
		string buf,temp,name;
		//read
		fin.open(pdb.c_str(), ios::in);
		if(fin.fail()!=0)
		{
			fprintf(stderr,"pdb_file %s not found!!\n",pdb.c_str());
			return -1;
		}
		int len;
		int count=0;
		for(;;)
		{
			if(!getline(fin,buf,'\n'))break;
			len=(int)buf.length();
			if(len<3)continue;
			//check TER
			temp=buf.substr(0,3);
			if(temp=="TER"||temp=="END")break;
			//check ATOM
			if(len<4)continue;
			temp=buf.substr(0,4);
			if(temp!="ATOM"&&temp!="HETA")continue;
			//check CA
			temp=buf.substr(13,2);
			if(temp!="CA")continue;
			//record name
			name=buf.substr(21,6);
			iter = ws_mapping.find(name);
			if(iter != ws_mapping.end())continue;
			count++;
			ws_mapping.insert(map < string, int >::value_type(name, count));
			//record ami
			temp=buf.substr(17,3);
			ami[count-1]=WWW_Three2One_III(temp.c_str());
			//record xyz
			temp=buf.substr(30,8);
			xyz[count-1].X=atof(temp.c_str());
			temp=buf.substr(38,8);
			xyz[count-1].Y=atof(temp.c_str());
			temp=buf.substr(46,8);
			xyz[count-1].Z=atof(temp.c_str());
		}
		//return 
		ami[count]='\0';
		return count;
	}

	//--------- TPL_To_XYZ ----------//
	int TPL_To_XYZ(string &tpl_file,char *ami,XYZ *xyz)
	{
		ifstream fin;
		string buf,temp;
		//read
		fin.open(tpl_file.c_str(), ios::in);
		if(fin.fail()!=0)
		{
			fprintf(stderr,"no such file! %s \n",tpl_file.c_str());
			return -1;
		}
		//skip
		char chain='A';
		for(;;)
		{
			if(!getline(fin,buf,'\n'))
			{
				fprintf(stderr,"file bad! %s \n",tpl_file.c_str());
				return -1;
			}
			istringstream www(buf);
			www>>temp;
			if(temp=="Chain")
			{
				www>>temp>>temp>>temp;
				istringstream sss(temp);
				string tomp;
				getline(sss,tomp,'.');
				if(tomp!="")chain=tomp[0];
			}
			if(buf=="//////////// Features ")break;
		}
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file bad! %s \n",tpl_file.c_str());
			return -1;
		}
		//process
		int miss;
		int resnum;
		string amino;
		double xca,xcb,yca,ycb,zca,zcb;
		const char *dummyaa;
		int count=0;
		for(;;)
		{
			if(!getline(fin,buf,'\n'))
			{
				fprintf(stderr,"file bad! %s \n",tpl_file.c_str());
				return -1;
			}
			if(buf=="")break;
			istringstream www(buf);
			www>>resnum>>amino>>miss;
			if(miss==1)continue;
			www>>temp>>temp>>temp>>temp>>temp>>temp;
			www>>xca>>yca>>zca>>xcb>>ycb>>zcb;
			//output
			ami[count]=amino[0];
			xyz[count].X=xca;
			xyz[count].Y=yca;
			xyz[count].Z=zca;
			//count
			count++;
		}
		//return 
		ami[count]='\0';
		return count;
	}

	//--------- TEMPLATE_To_XYZ ----------//
	int TEMPLATE_To_XYZ(TEMPLATE *tpl ,char *ami,XYZ *xyz)
	{
		int i;
		int count=0;
		for(i=0;i<tpl->tLength;i++)
		{
			if(tpl->isMissing[i]!=0)continue;
			xyz[count].X=tpl->CA[i][0];
			xyz[count].Y=tpl->CA[i][1];
			xyz[count].Z=tpl->CA[i][2];
			ami[count]=tpl->sequence[i];
			count++;
		}
		//return 
		ami[count]='\0';
		return count;
	}

	//========= vice process ===========//
	//-> [note]: simp version assumes that two sequences in the alignment are identical to ami1 and ami2
	double uGDT_Calc_Simp(XYZ *mol1,XYZ *mol2,char *ami1,char *ami2,
			vector<pair<int, int> > &alignment)
	{
		//moln
		int TM_MOLN1=(int)strlen(ami1);
		int TM_MOLN2=(int)strlen(ami2);
		//init
		XYZ *xyz1=new XYZ[TM_MOLN1];
		XYZ *xyz2=new XYZ[TM_MOLN2];
		int lali=0;
		for(int i=0;i<(int)alignment.size();i++)
		{
			int ii=alignment[i].first;
			int jj=alignment[i].second;
			if(ii>0 && jj>0 && ii<=TM_MOLN1 && jj<=TM_MOLN2)
			{
				xyz1[lali]=mol1[ii-1];
				xyz2[lali]=mol2[jj-1];
				lali++;
			}
		}
		//calc
		double Ret_Sco[8];
		TM_score tm_score(lali);
		tm_score.Calc_TM_d0(TM_MOLN1<TM_MOLN2?TM_MOLN1:TM_MOLN2);
		double tms=tm_score.Calc_TM_Score(xyz1,xyz2,lali,tm_score.d0,tm_score.d8,0,0,Ret_Sco);
		double ugdt=1.0*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4])/4.0;
		//delete
		delete [] xyz1;
		delete [] xyz2;
		//return
		return ugdt;
	}

	//-> [note]: complex version will run Extract_Alignment procedure
	double uGDT_Calc_Complex(XYZ *mol1,XYZ *mol2,char *ami1,char *ami2,
			vector<pair<int, int> > &alignment_in, string &nam1_content,string &nam2_content)
	{
		//extract alignment
		string nam1_pdb=ami1;
		string nam2_pdb=ami2;
		vector<pair<int, int> > alignment;
		vector <int> ali1_seq_pdb;
		vector <int> ali2_seq_pdb;
		Extract_Alignment(nam1_content,nam1_pdb,nam2_content,nam2_pdb,alignment_in,alignment,
				ali1_seq_pdb,ali2_seq_pdb);
		//moln
		int TM_MOLN1=(int)strlen(ami1);
		int TM_MOLN2=(int)strlen(ami2);
		//init
		XYZ *xyz1=new XYZ[TM_MOLN1];
		XYZ *xyz2=new XYZ[TM_MOLN2];
		int lali=0;
		for(int i=0;i<(int)alignment.size();i++)
		{
			int ii=alignment[i].first;
			int jj=alignment[i].second;
			if(ii>0 && jj>0 && ii<=TM_MOLN1 && jj<=TM_MOLN2)
			{
				xyz1[lali]=mol1[ii-1];
				xyz2[lali]=mol2[jj-1];
				lali++;
			}
		}
		//calc
		double Ret_Sco[8];
		TM_score tm_score(lali);
		tm_score.Calc_TM_d0(TM_MOLN1<TM_MOLN2?TM_MOLN1:TM_MOLN2);
		double tms=tm_score.Calc_TM_Score(xyz1,xyz2,lali,tm_score.d0,tm_score.d8,0,0,Ret_Sco);
		double ugdt=1.0*(Ret_Sco[1]+Ret_Sco[2]+Ret_Sco[3]+Ret_Sco[4])/4.0;
		//delete
		delete [] xyz1;
		delete [] xyz2;
		//return
		return ugdt;
	}

	//========= main process ===========//
	//-> uGDT calc from two TPLs
	double uGDT_Calc_From_Two_TPL_Simp(TEMPLATE *tpl1, TEMPLATE *tpl2,vector<pair<int, int> > &alignment)
	{
		//data
		int moln1=tpl1->tLength;
		int moln2=tpl2->tLength;
		XYZ *xyz1=new XYZ[moln1];
		XYZ *xyz2=new XYZ[moln2];
		char *ami1=new char[moln1+1];
		char *ami2=new char[moln2+1];
		//get to XYZ
		TEMPLATE_To_XYZ(tpl1,ami1,xyz1);
		TEMPLATE_To_XYZ(tpl2,ami2,xyz2);
		//calculate ugdt
		double ugdt=uGDT_Calc_Simp(xyz1,xyz2,ami1,ami2,alignment);
		//delete
		delete [] xyz1;
		delete [] xyz2;
		delete [] ami1;
		delete [] ami2;
		//return
		return ugdt;
	}

	//-> complex version
	double uGDT_Calc_From_Two_TPL_Complex(TEMPLATE *tpl1, TEMPLATE *tpl2,
			vector<pair<int, int> > &alignment_in,string &nam1_content,string &nam2_content)
	{
		//data
		int moln1=tpl1->tLength;
		int moln2=tpl2->tLength;
		XYZ *xyz1=new XYZ[moln1];
		XYZ *xyz2=new XYZ[moln2];
		char *ami1=new char[moln1+1];
		char *ami2=new char[moln2+1];
		//get to XYZ
		TEMPLATE_To_XYZ(tpl1,ami1,xyz1);
		TEMPLATE_To_XYZ(tpl2,ami2,xyz2);
		//calculate ugdt
		double ugdt=uGDT_Calc_Complex(xyz1,xyz2,ami1,ami2,alignment_in,nam1_content,nam2_content);
		//delete
		delete [] xyz1;
		delete [] xyz2;
		delete [] ami1;
		delete [] ami2;
		//return
		return ugdt;
	}

	//-> pure TPL version
	double uGDT_Calc_From_Two_TPL_Pure(TEMPLATE *tpl1, TEMPLATE *tpl2,vector<pair<int, int> > &alignment_in)
	{
		//data
		int moln1=tpl1->tLength;
		int moln2=tpl2->tLength;
		XYZ *xyz1=new XYZ[moln1];
		XYZ *xyz2=new XYZ[moln2];
		char *ami1=new char[moln1+1];
		char *ami2=new char[moln2+1];
		//get to XYZ
		TEMPLATE_To_XYZ(tpl1,ami1,xyz1);
		TEMPLATE_To_XYZ(tpl2,ami2,xyz2);
		//calculate ugdt
		double ugdt=uGDT_Calc_Complex(xyz1,xyz2,ami1,ami2,alignment_in,tpl1->sequence,tpl2->sequence);
		//delete
		delete [] xyz1;
		delete [] xyz2;
		delete [] ami1;
		delete [] ami2;
		//return
		return ugdt;
	}
}


//---------- main ----------//
/*
 *int main(int argc,char **argv)
 *{
 *
 *    //---- Calc_uGDT_Score from TPL ----//
 *    {
 *        if(argc<4)
 *        {
 *            fprintf(stderr,"Calc_TPL_uGDT <tpl1> <tpl2> <fasta_file> \n");
 *            exit(-1);
 *        }
 *        string tpl_file1=argv[1];
 *        string tpl_file2=argv[2];
 *        string fasta_file=argv[3];
 *        //load
 *        //-> tpl1
 *        string tpl_name1,tpl_root1;
 *        getBaseName(tpl_file1,tpl_name1,'/','.');
 *        getRootName(tpl_file1,tpl_root1,'/');
 *        TEMPLATE* tpl1;
 *        tpl1 = new TEMPLATE(tpl_name1,tpl_root1,1);
 *        if(tpl1->failure==0)exit(-1);
 *        //-> tpl2
 *        string tpl_name2,tpl_root2;
 *        getBaseName(tpl_file2,tpl_name2,'/','.');
 *        getRootName(tpl_file2,tpl_root2,'/');
 *        TEMPLATE* tpl2;
 *        tpl2 = new TEMPLATE(tpl_name2,tpl_root2,1);
 *        if(tpl2->failure==0)exit(-1);
 *        //-> alignment
 *        vector<pair<int, int> > alignment;
 *        string nam1_content,nam2_content;
 *        string nam1_full,nam2_full;
 *        string nam1,nam2;
 *        int retv=uGDT_Calc::ReadToFile_FASTA(fasta_file,alignment,nam1_content,nam2_content,nam1_full,nam2_full,nam1,nam2);
 *        if(retv!=1)exit(-1);
 *        //process
 *        double ugdt=uGDT_Calc::uGDT_Calc_From_Two_TPL_Pure(tpl1, tpl2, alignment);
 *        printf("%lf \n",ugdt);
 *        //exit
 *        exit(0);
 *    }
 *}
 */

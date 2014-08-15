#pragma once
#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "template.h"
#include "seq.h"
#include "profile.h"
#include "CNFalign_Basic.h"
using namespace std;

namespace CalcFeature {
	//===================== CalcFeatures ==================//
	void CalcFeatures(vector<pair<int,int> > & alignment,vector<double> &features, TEMPLATE *temp,SEQUENCE *seq);
}

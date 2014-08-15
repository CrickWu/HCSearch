#pragma once
#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include "TM_score.h"
#include "template.h"
#include "seq.h"
#include "profile.h"
using namespace std;

namespace uGDT_Calc {
	double uGDT_Calc_From_Two_TPL_Pure(TEMPLATE *tpl1, TEMPLATE *tpl2, vector<pair<int, int> > &alignment);

}

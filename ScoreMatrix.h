#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include "Score.h"
using namespace std;

class ScoreMatrix {

  int layers;
  int rows;
  int cols;
  vector<Score> data; 
  void PrintVal (ostream &outfile, const Score &value) const;

 public:
  ScoreMatrix (int layers, int rows, int cols);
  ScoreMatrix (const ScoreMatrix &m);
  ~ScoreMatrix();

  void Fill (const Score &value);
  void PrintLayer (ostream &outfile, int layer) const;
  void Print (ostream &outfile) const;
  Score ComputeSum() const;

//------ access matrix element -------//
inline Score &operator() (int layer, int row, int col){
//    Assert (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
//    Assert (0 <= row && row < rows, "Requested row out-of-bounds.");
//    Assert (0 <= col && col < cols, "Requested column out-of-bounds.");
    return data[(row * cols + col) * layers + layer];
} 
inline Score operator() (int layer, int row, int col) const {
//    Assert (0 <= layer && layer < layers, "Requested layer out-of-bounds.");
//    Assert (0 <= row && row < rows, "Requested row out-of-bounds.");
//    Assert (0 <= col && col < cols, "Requested column out-of-bounds.");
    return data[(row * cols + col) * layers + layer];
}
  int GetNumLayers() const;
  int GetNumRows() const;
  int GetNumCols() const;
};

#endif

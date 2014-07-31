#include "ScoreMatrix.h"

//----- constructor and destructor ----- //
ScoreMatrix::ScoreMatrix (int layers, int rows, int cols) : 
  layers(layers), rows(rows), cols(cols) , data (layers * rows * cols) {
//  Assert (layers >= 0, "Number of layers in matrix must be positive.");
//  Assert (rows >= 0, "Number of rows in matrix must be positive.");
//  Assert (cols >= 0, "Number of columns in matrix must be positive.");
}
ScoreMatrix::~ScoreMatrix(){	
}

//------ copy constructor -------//
ScoreMatrix::ScoreMatrix (const ScoreMatrix &m) :
  layers (m.layers), rows(m.rows), cols(m.cols), data(m.data){
//  Assert (layers >= 0, "Number of layers in matrix must be positive.");
//  Assert (rows >= 0, "Number of rows in matrix must be positive.");
//  Assert (cols >= 0, "Number of columns in matrix must be positive.");
}


//-------- Fill all matrix values --------//
void ScoreMatrix::Fill (const Score &value){
  for (int i = 0; i < layers * rows * cols; i++)
    data[i] = value;
}

//-------- Printing utility function -----//
void ScoreMatrix::PrintVal (ostream &outfile, const Score &value) const {
  outfile << setw(13) << setprecision(10) << setiosflags(ios::fixed | ios::showpoint) 
	  << value << resetiosflags(ios::fixed | ios::showpoint);
}
void ScoreMatrix::PrintLayer (ostream &outfile, int layer) const {
  for (int i = 0; i < rows; i++){
    outfile << (i == 0 ? "[" : " ") << "[";
    for (int j = 0; j < cols; j++){
      if (j > 0) outfile << ", ";
      PrintVal (outfile, operator()(layer,i,j));
    }
    outfile << "]" << (i == rows-1 ? "]" : "") << endl;
  }
  outfile << endl;
}
void ScoreMatrix::Print (ostream &outfile) const {
  for (int i = 0; i < layers; i++)
    PrintLayer (outfile, i);
}

//-------- Compute sum of all entries in matrix ----//
Score ScoreMatrix::ComputeSum() const {
  Score total = 0;
  for (int i = 0; i < layers * rows * cols; i++)
    total += data[i];
  return total;
}
/* Return number of matrix layers */
int ScoreMatrix::GetNumLayers() const {
  return layers;
}
/* Return number of matrix rows */
int ScoreMatrix::GetNumRows() const {
  return rows;
}
/* Return number of matrix columns */
int ScoreMatrix::GetNumCols() const {
  return cols;
}

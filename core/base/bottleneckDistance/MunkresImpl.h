#pragma once

#include "Munkres.h"

namespace ttk {

int Munkres::run(std::vector<matchingTuple> &matchings) {
  int step = 1;
  int iter = 0;
  int maxIter = 100000;
  bool done = false;
  Timer t;

  std::vector<std::vector<dataType>> inputMatrix(
    rowSize, std::vector<dataType>(colSize));
  copyInputMatrix(inputMatrix);

  while(!done) {
    ++iter;
    this->printMsg(
      "Step " + std::to_string(step) + ", Iteration " + std::to_string(iter),
      debug::Priority::DETAIL);

    if(iter > 20 && (iter % (int)std::round((double)maxIter / 5.0) == 0)) {
      double progress = std::round(100.0 * (double)iter / (double)maxIter);
      this->printMsg("Progress", progress / 100.0, t.getElapsedTime());
    }

    if(iter > maxIter) {
      // showCostMatrix();
      // showMaskMatrix();

      this->printMsg("Failed to converge after " + std::to_string(maxIter)
                     + " iterations. Aborting.");

      step = 7;
      // Abort. Still found something
      // though not optimal.
    }

    // Show intermediary matrices:
    // showCostMatrix();
    // showMaskMatrix();

    switch(step) {
      case 1:
        stepOne(step);
        break;
      case 2:
        stepTwo(step);
        break;
      case 3:
        stepThree(step);
        break;
      case 4:
        stepFour(step);
        break;
      case 5:
        stepFive(step);
        break;
      case 6:
        stepSix(step);
        break;
      case 7:
        stepSeven(step);
        done = true;
        break;
      default:
        break;
    }
  }

  this->computeAffectationCost(inputMatrix);
  this->affect(matchings, inputMatrix);
  this->clear();

  return 0;
}

// Preprocess cost matrix.
int Munkres::stepOne(int &step) // ~ 0% perf
{
  double minInCol;
  std::vector<std::vector<dataType>> *C
    = (std::vector<std::vector<dataType>> *)Cptr;

  // Benefit from the matrix sparsity.
  dataType maxVal = std::numeric_limits<dataType>::max();
  for(int r = 0; r < rowSize - 1; ++r) {
    rowLimitsPlus[r] = -1;
    rowLimitsMinus[r] = -1;
  }
  for(int c = 0; c < colSize - 1; ++c) {
    colLimitsPlus[c] = -1;
    colLimitsMinus[c] = -1;
  }

  int droppedMinus = 0;
  int droppedPlus = 0;

  for(int r = 0; r < rowSize - 1; ++r) {
    for(int c = 0; c < colSize - 1; ++c)
      if((*C)[r][c] != maxVal) {
        rowLimitsMinus[r] = c; // Included
        break;
      }
    if(rowLimitsMinus[r] == -1) {
      ++droppedMinus;
      rowLimitsMinus[r] = 0;
    } // Included

    for(int c = colSize - 2; c >= 0; --c)
      if((*C)[r][c] != maxVal) {
        rowLimitsPlus[r] = c + 1; // Not included
        break;
      }
    if(rowLimitsPlus[r] == -1) {
      ++droppedPlus;
      rowLimitsPlus[r] = colSize - 1;
    } // Not included
  }

  if(droppedMinus > 0) {
    this->printMsg(
      "Unexpected non-assignable row [minus], dropping optimisation for "
        + std::to_string(droppedMinus) + " row(s).",
      debug::Priority::DETAIL);
  }

  if(droppedPlus > 0) {
    this->printMsg(
      "Unexpected non-assignable row [plus], dropping optimisation for "
        + std::to_string(droppedPlus) + " row(s).",
      debug::Priority::DETAIL);
  }

  droppedMinus = 0;
  droppedPlus = 0;

  for(int c = 0; c < colSize - 1; ++c) {
    for(int r = 0; r < rowSize - 1; ++r)
      if((*C)[r][c] != maxVal) {
        colLimitsMinus[c] = r; // Inclusive
        break;
      }
    for(int r = rowSize - 1; r >= 0; --r)
      if((*C)[r][c] != maxVal) {
        colLimitsPlus[c] = r + 1; // Exclusive.
        break;
      }
    if(colLimitsPlus[c] == -1) {
      ++droppedPlus;
      colLimitsMinus[c] = 0;
    }
    if(colLimitsMinus[c] == -1) {
      ++droppedMinus;
      colLimitsMinus[c] = rowSize;
    }
  }

  if(droppedMinus > 0) {
    this->printMsg(
      "Unexpected non-assignable row [minus], dropping optimisation for "
        + std::to_string(droppedMinus) + " row(s).",
      debug::Priority::DETAIL);
  }

  if(droppedPlus > 0) {
    this->printMsg(
      "Unexpected non-assignable row [plus], dropping optimisation for "
        + std::to_string(droppedPlus) + " row(s).",
      debug::Priority::DETAIL);
  }

  rowLimitsMinus[rowSize - 1] = 0;
  rowLimitsPlus[rowSize - 1] = colSize - 1;

  // Remove last column (except the last element) from all other columns.
  // The last column will then be ignored during the solving.
  for(int r = 0; r < rowSize - 1; ++r) {
    dataType lastElement = (*C)[r][colSize - 1];
    for(int c = 0; c < colSize - 1; ++c) {
      (*C)[r][c] -= lastElement;
    }
  }

  // Substract minimum value in every column except the last.
  for(int c = 0; c < colSize - 1; ++c) {
    minInCol = (*C)[0][c];

    for(int r = 0; r < rowSize; ++r)
      if((*C)[r][c] < minInCol)
        minInCol = (*C)[r][c];

    for(int r = 0; r < rowSize; ++r)
      (*C)[r][c] -= minInCol;
  }

  step = 2;
  return 0;
}

// Find a zero in the matrix,
// star it if it is the only one in its row and col.
int Munkres::stepTwo(int &step) // ~ 0% perf
{
  std::vector<std::vector<dataType>> *C
    = (std::vector<std::vector<dataType>> *)Cptr;

  for(int r = 0; r < rowSize - 1; ++r) {
    for(int c = 0; c < colSize - 1; ++c) {
      if(!rowCover[r] && !colCover[c] && isZero((*C)[r][c])) {
        M[r][c] = 1;
        // Temporarily cover row and column to find independent zeros.
        rowCover[r] = true;
        colCover[c] = true;
      }
    }

    // Don't account for last column.
  }

  for(int c = 0; c < colSize - 1; ++c)
    if(isZero((*C)[rowSize - 1][c]) && !colCover[c]) {
      M[rowSize - 1][c] = 1;
      // Don't ban last row where elements are all independent.
      colCover[c] = true;
    }

  // Remove coverings (temporarily used to find independent zeros).
  for(int r = 0; r < rowSize; ++r)
    rowCover[r] = false;

  for(int c = 0; c < colSize - 1; ++c)
    colCover[c] = false;

  step = 3;
  return 0;
}

// Check column coverings.
// If all columns are starred (1 star only per column is possible)
// then the algorithm is terminated.
int Munkres::stepThree(int &step) // ~ 10% perf
{
  for(int r = 0; r < rowSize; ++r) {
    int start = rowLimitsMinus[r];
    int end = rowLimitsPlus[r];
    for(int c = start; c < end; ++c)
      if(M[r][c] == 1)
        colCover[c] = true;
  }

  int processedCols = 0;

  for(int c = 0; c < colSize - 1; ++c)
    if(colCover[c])
      ++processedCols;

  if(processedCols >= colSize - 1)
    step = 7; // end algorithm
  else
    step = 4; // follow prime scheme
  return 0;
}

// Find a non covered zero, prime it
// . if current row is last or has no starred zero -> step 5
// . else, cover row and uncover the col with a star
// Repeat until there are no uncovered zero left
// Save smallest uncovered value then -> step 6
int Munkres::stepFour(int &step) // ~ 45% perf
{
  int row = -1;
  int col = -1;
  bool done = false;

  while(!done) {
    findZero(row, col);

    if(row == -1) {
      done = true;
      step = 6;
    }

    else {
      M[row][col] = 2;
      int colOfStarInRow = findStarInRow(row);
      // If a star was found and it is not in the last row
      if(colOfStarInRow > -1 && row < rowSize - 1) {
        rowCover[row] = true;
        colCover[colOfStarInRow] = false;
      }

      else {
        done = true;
        step = 5;
        pathRow0 = row;
        pathCol0 = col;
      }
    }
  }

  return 0;
}

int Munkres::findStarInRow(int row) {
  int start = rowLimitsMinus[row];
  int end = rowLimitsPlus[row];
  for(int c = start; c < end; ++c)
    if(M[row][c] == 1)
      return c;
  return -1;
}

int Munkres::findZero(int &row, int &col) {
  std::vector<std::vector<dataType>> *C
    = (std::vector<std::vector<dataType>> *)Cptr;

  row = -1;
  col = -1;

  while(createdZeros.size() > 0) {
    std::pair<int, int> zero = createdZeros.back();
    int f = zero.first;
    int s = zero.second;
    createdZeros.pop_back();
    if(!rowCover[f] && !colCover[s]) {
      row = f;
      col = s;
      return 0;
    }
  }

  for(int r = 0; r < rowSize; ++r) {
    int start = rowLimitsMinus[r];
    int end = rowLimitsPlus[r];
    if(rowCover[r])
      continue;

    for(int c = start; c < end; ++c) {
      if(colCover[c])
        continue;
      if((*C)[r][c] == (dataType)0) {
        row = r;
        col = c;
        return 0;
      }
    }
  }

  this->printMsg("Zero not found.", debug::Priority::DETAIL);

  return 0;
}

// Make path of alternating primed and starred zeros
// 1. uncovered primed found at step 4
// 2. same column, starred (if any)
// 3. same row, primed (always one)
// 4. continue until a primed zero has no starred zero in its column
// Unstar each starred zero in the series, star each primed zero
// in the series,
// erase all primes, uncover every line, return to step 3.
int Munkres::stepFive(int &step) // ~ 10% perf
{
  {
    int r;
    int c;

    pathCount = 1;
    path[pathCount - 1][0] = pathRow0;
    path[pathCount - 1][1] = pathCol0;

    bool done = false;
    while(!done) {
      r = findStarInCol(path[pathCount - 1][1]);
      if(r == -1)
        done = true;

      else {
        ++pathCount;
        path[pathCount - 1][0] = r;
        path[pathCount - 1][1] = path[pathCount - 2][1];

        c = findPrimeInRow(path[pathCount - 1][0]);
        if(c == -1) {
          this->printWrn("Did not find an expected prime.");
        }
        ++pathCount;
        path[pathCount - 1][0] = path[pathCount - 2][0];
        path[pathCount - 1][1] = c;
      }
    }
  }

  // process path
  for(int p = 0; p < pathCount; ++p) {
    if(M[path[p][0]][path[p][1]] == 1)
      M[path[p][0]][path[p][1]] = 0;
    else
      M[path[p][0]][path[p][1]] = 1;
  }

  // clear covers
  for(int r = 0; r < rowSize; ++r)
    rowCover[r] = false;
  for(int c = 0; c < colSize - 1; ++c)
    colCover[c] = false;

  // erase primes
  for(int r = 0; r < rowSize; ++r) {
    int start = rowLimitsMinus[r];
    int end = rowLimitsPlus[r];
    for(int c = start; c < end; ++c)
      if(M[r][c] == 2)
        M[r][c] = 0;
  }

  step = 3;
  return 0;
}

int Munkres::findStarInCol(int col) {
  int start = colLimitsMinus[col];
  int end = colLimitsPlus[col];
  for(int r = start; r < end; ++r)
    if(M[r][col] == 1)
      return r;

  if(M[rowSize - 1][col] == 1)
    return (rowSize - 1);
  return -1;
}

int Munkres::findPrimeInRow(int row) {
  int start = rowLimitsMinus[row];
  int end = rowLimitsPlus[row];
  for(int c = start; c < end; ++c)
    if(M[row][c] == 2)
      return c;
  return -1;
}

// Add smallest value to every element of each covered row,
// subtract it from every element of each uncovered col.
// Return to step 4 without altering any stars/primes/covers.
int Munkres::stepSix(int &step) // ~ 35% perf
{
  std::vector<std::vector<dataType>> *C
    = (std::vector<std::vector<dataType>> *)Cptr;

  dataType minVal = std::numeric_limits<dataType>::max();

  // find smallest
  for(int r = 0; r < rowSize; ++r) {
    if(rowCover[r])
      continue;

    int start = rowLimitsMinus[r];
    int end = rowLimitsPlus[r];

    for(int c = start; c < end; ++c) {
      if(colCover[c])
        continue;
      if((*C)[r][c] < minVal)
        minVal = (*C)[r][c];
    }
  }

  createdZeros.clear();

  // add and subtract
  for(int r = 0; r < rowSize; ++r) {

    int start = rowLimitsMinus[r];
    int end = rowLimitsPlus[r];

    for(int c = start; c < end; ++c) {
      if(rowCover[r])
        (*C)[r][c] = (*C)[r][c] + minVal;
      if(!colCover[c]) {
        (*C)[r][c] = (*C)[r][c] - minVal;
        if(isZero((*C)[r][c])) {
          createdZeros.emplace_back(r, c);
        }
      }
    }
  }

  step = 4;
  return 0;
}

int Munkres::stepSeven(int &step) {
  this->printMsg("Step 7 over.", debug::Priority::DETAIL);
  return 0;
}

int Munkres::affect(std::vector<matchingTuple> &matchings,
                    const std::vector<std::vector<dataType>> &C) {
  int nbC = colSize;
  int nbR = rowSize;

  matchings.clear();

  for(int r = 0; r < nbR; ++r)
    for(int c = 0; c < nbC; ++c)
      if(M[r][c] == 1) {
        matchingTuple t = std::make_tuple(r, c, C[r][c]);
        matchings.push_back(t);
        // Use row cover to match to last column diagonal.
        if(r < nbR - 1)
          rowCover[r] = true;
      }

  // Clear row cover
  for(int r = 0; r < nbR - 1; ++r) {
    // Match to diagonal.
    if(!rowCover[r]) {
      matchingTuple t = std::make_tuple(r, nbC - 1, C[r][nbC - 1]);
      matchings.push_back(t);
    }
    // Ensure row covers are cleared.
    else {
      rowCover[r] = false;
    }
  }

  return 0;
}

int Munkres::computeAffectationCost(
  const std::vector<std::vector<dataType>> &C) {
  int nbC = colSize;
  int nbR = rowSize;

  dataType total = 0;

  for(int r = 0; r < nbR; ++r)
    for(int c = 0; c < nbC; ++c)
      if(M[r][c] == 1) {
        total += C[r][c];
      }

  this->printMsg("Total cost: " + std::to_string(total));

  return 0;
}

} // namespace ttk

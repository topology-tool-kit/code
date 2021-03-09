/// \ingroup base
/// \class ttk::Munkres
/// \author Maxime Soler <soler.maxime@total.com>

#pragma once

#include <Debug.h>
#include <cmath>
#include <iostream>
#include <limits>

typedef double dataType;

namespace ttk {

  class Munkres : public Debug {

  public:
    Munkres() {
      this->setDebugMsgPrefix("Munkres");
    }

    ~Munkres(){};

    int run(std::vector<matchingTuple> &matchings);

    inline void clear() {
      pathCount = 0;
      rowSize = 0;
      colSize = 0;
      createdZeros.clear();
    }

    inline void clearMatrix() {
      std::vector<std::vector<dataType>> C
        = *((std::vector<std::vector<dataType>> *)Cptr);
      for(int r = 0, rS0 = rowSize; r < rS0; ++r)
        for(int c = 0, cS0 = colSize; c < cS0; ++c)
          C[r][c] = 0.0;
    }

    inline int setInput(int rowSize_, int colSize_, void *C_) {
      rowSize = rowSize_;
      colSize = colSize_;

      createdZeros.clear();

      Cptr = C_;

      auto r = (unsigned long)rowSize_;
      auto c = (unsigned long)colSize_;

      rowCover.resize(r);
      colCover.resize(c);

      rowLimitsMinus.resize(r);
      rowLimitsPlus.resize(r);
      colLimitsMinus.resize(c);
      colLimitsPlus.resize(c);

      M.resize(r);
      for(int row = 0; row < rowSize_; ++row)
        M[row].resize(c);

      int nbPaths = 1 + colSize_ + rowSize_;
      path.resize((unsigned long)nbPaths);
      for(int p = 0; p < nbPaths; ++p)
        path[p].resize(2);

      resetMasks();

      return 0;
    }

    inline void showCostMatrix() {
      auto C = (std::vector<std::vector<dataType>> *)Cptr;
      std::stringstream msg;
      for(int r = 0; r < rowSize; ++r) {
        msg << std::endl << "  ";
        for(int c = 0; c < colSize; ++c)
          msg << std::fixed << std::setprecision(3) << (*C)[r][c] << " ";
      }

      this->printMsg(msg.str(), debug::Priority::DETAIL);
    }

    inline void showMaskMatrix() {
      std::stringstream msg;
      msg << std::endl << "   | ";
      for(int c = 0; c < colSize; ++c) {
        msg << colCover[c] << "  ";
      }
      msg << std::endl << "   | ";
      for(int c = 0; c < colSize; ++c) {
        msg << "---";
      }
      msg << std::endl;
      for(int r = 0; r < rowSize; ++r) {
        msg << " " << rowCover[r] << " | ";
        for(int c = 0; c < colSize; ++c)
          msg << M[r][c] << "  ";
        msg << std::endl;
      }
      this->printMsg(msg.str(), debug::Priority::DETAIL);
    }

  private:
    void *Cptr;

    std::vector<std::vector<int>> M;
    std::vector<bool> rowCover;
    std::vector<bool> colCover;

    std::vector<int> rowLimitsMinus;
    std::vector<int> rowLimitsPlus;
    std::vector<int> colLimitsMinus;
    std::vector<int> colLimitsPlus;

    std::vector<std::vector<int>> path;
    std::vector<std::pair<int, int>> createdZeros;

    int pathRow0;
    int pathCol0;
    int pathCount = 0;

    int rowSize = 0;
    int colSize = 0;

    // internal methods
    int stepOne(int &step);

    int stepTwo(int &step);

    int stepThree(int &step);

    int stepFour(int &step);
    int findZero(int &row, int &col);
    int findStarInRow(int row);

    int stepFive(int &step);
    int findStarInCol(int col);
    int findPrimeInRow(int row);

    int stepSix(int &step);

    int stepSeven(int &step);

    int affect(std::vector<matchingTuple> &matchings,
               const std::vector<std::vector<dataType>> &C);

    int computeAffectationCost(const std::vector<std::vector<dataType>> &C);

    inline bool isZero(dataType t) {
      // return std::abs((double) t) < 1e-15;
      return t == 0.;
    }

    inline int resetMasks() {
      for(int r = 0; r < rowSize; ++r) {
        rowCover[r] = false;
        for(int c = 0; c < colSize; ++c)
          M[r][c] = 0;
      }
      for(int c = 0; c < colSize; ++c)
        colCover[c] = false;
      return 0;
    }

    inline int copyInputMatrix(std::vector<std::vector<dataType>> &saveInput) {
      auto C = (std::vector<std::vector<dataType>> *)Cptr;
      int rS = rowSize;
      int cS = colSize;

      for(int r = 0; r < rS; ++r)
        for(int c = 0; c < cS; ++c)
          saveInput[r][c] = (*C)[r][c];

      return 0;
    }
  };

} // namespace ttk

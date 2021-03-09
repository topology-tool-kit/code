#pragma once

#include "BottleneckDistance.h"

namespace ttk {

constexpr unsigned long long str2int(const char *str, int h = 0) {
  return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}

int BottleneckDistance::execute(const bool usePersistenceMetric) {
  Timer t;

  bool fromParaView = pvAlgorithm_ >= 0;
  if(fromParaView) {
    switch(pvAlgorithm_) {
      case 0:
        this->printMsg("Solving with the TTK approach");
        this->computeBottleneck(
          *static_cast<const std::vector<diagramTuple> *>(outputCT1_),
          *static_cast<const std::vector<diagramTuple> *>(outputCT2_),
          *static_cast<std::vector<matchingTuple> *>(matchings_),
          usePersistenceMetric);
        break;
      case 1: {
        std::stringstream msg;
        this->printMsg("Solving with the legacy Dionysus exact approach.");
        this->printErr("Not supported");
      } break;
      case 2: {
        this->printMsg(
          "Solving with the approximate Dionysus geometric approach.");
        this->printErr("Not supported");
      } break;
      case 3: {
        this->printMsg("Solving with the parallel TTK approach");
        this->printErr("Not supported");
      } break;
      case 4: {
        this->printMsg("Benchmarking");
        this->printErr("Not supported");
      } break;
      default: {
        this->printErr("You must specify a valid assignment algorithm.");
      }
    }

  } else {
    switch(str2int(algorithm_.c_str())) {
      case str2int("0"):
      case str2int("ttk"):
        this->printMsg("Solving with the TTK approach");
        this->computeBottleneck(
          *static_cast<const std::vector<diagramTuple> *>(outputCT1_),
          *static_cast<const std::vector<diagramTuple> *>(outputCT2_),
          *static_cast<std::vector<matchingTuple> *>(matchings_),
          usePersistenceMetric);
        break;
      case str2int("1"):
      case str2int("legacy"): {
        this->printMsg("Solving with the legacy Dionysus exact approach.");
        this->printErr("Not supported");
      } break;
      case str2int("2"):
      case str2int("geometric"): {
        this->printMsg(
          "Solving with the approximate Dionysus geometric approach.");
        this->printErr("Not supported");
      } break;
      case str2int("3"):
      case str2int("parallel"): {
        this->printMsg("Solving with the parallel TTK approach");
        this->printErr("Not supported");
      } break;
      case str2int("bench"): {
        this->printMsg("Benchmarking");
        this->printErr("Not supported");
      } break;
      default: {
        this->printErr("You must specify a valid assignment algorithm.");
      }
    }
  }

  this->printMsg("Complete", 1, t.getElapsedTime(), threadNumber_);

  return 0;
}

double BottleneckDistance::computeGeometricalRange(
  const std::vector<diagramTuple> &CTDiagram1,
  const std::vector<diagramTuple> &CTDiagram2,
  const int d1Size,
  const int d2Size) const {
  float minX1, maxX1, minY1, maxY1, minZ1, maxZ1;
  float minX2, maxX2, minY2, maxY2, minZ2, maxZ2;
  float minX, minY, minZ, maxX, maxY, maxZ;
  minX1 = minY1 = minZ1 = minX2 = minY2 = minZ2
    = std::numeric_limits<float>::max();
  maxX1 = maxY1 = maxZ1 = maxX2 = maxY2 = maxZ2
    = std::numeric_limits<float>::min();

  for(int i = 0; i < d1Size; ++i) {
    const diagramTuple &t = CTDiagram1[i];
    float xa = t.birthPoint.x, ya = t.birthPoint.y, za = t.birthPoint.z;
    float xb = t.deathPoint.x, yb = t.deathPoint.y, zb = t.deathPoint.z;
    minX1 = std::min(std::min(minX1, xa), xb);
    minY1 = std::min(std::min(minY1, ya), yb);
    minZ1 = std::min(std::min(minZ1, za), zb);
    maxX1 = std::max(std::max(maxX1, xa), xb);
    maxY1 = std::max(std::max(maxY1, ya), yb);
    maxZ1 = std::max(std::max(maxZ1, za), zb);
  }

  for(int i = 0; i < d2Size; ++i) {
    const diagramTuple &t = CTDiagram2[i];
    float xa = t.birthPoint.x, ya = t.birthPoint.y, za = t.birthPoint.z;
    float xb = t.deathPoint.x, yb = t.deathPoint.y, zb = t.deathPoint.z;
    minX2 = std::min(std::min(minX2, xa), xb);
    minY2 = std::min(std::min(minY2, ya), yb);
    minZ2 = std::min(std::min(minZ2, za), zb);
    maxX2 = std::max(std::max(maxX2, xa), xb);
    maxY2 = std::max(std::max(maxY2, ya), yb);
    maxZ2 = std::max(std::max(maxZ2, za), zb);
  }

  minX = std::min(minX1, minX2);
  maxX = std::max(maxX1, maxX2);
  minY = std::min(minY1, minY2);
  maxY = std::max(maxY1, maxY2);
  minZ = std::min(minZ1, minZ2);
  maxZ = std::max(maxZ1, maxZ2);

  return sqrt(Geometry::pow(maxX - minX, 2) + Geometry::pow(maxY - minY, 2)
              + Geometry::pow(maxZ - minZ, 2));
}

double BottleneckDistance::computeMinimumRelevantPersistence(
  const std::vector<diagramTuple> &CTDiagram1,
  const std::vector<diagramTuple> &CTDiagram2,
  const int d1Size,
  const int d2Size) const {
  double sp = zeroThreshold_;
  double s = sp > 0.0 && sp < 100.0 ? sp / 100.0 : 0;

  std::vector<dataType> toSort;
  for(int i = 0; i < d1Size; ++i) {
    const diagramTuple &t = CTDiagram1[i];
    dataType persistence = abs<dataType>(t.persistence);
    toSort.push_back(persistence);
  }
  for(int i = 0; i < d2Size; ++i) {
    const diagramTuple &t = CTDiagram2[i];
    dataType persistence = abs<dataType>(t.persistence);
    toSort.push_back(persistence);
  }
  sort(toSort.begin(), toSort.end());

  double minVal = toSort.at(0);
  double maxVal = toSort.at(toSort.size() - 1);
  s *= (maxVal - minVal);

  // double epsilon = 0.0000001;
  // int largeSize = 2000;
  // dataType zeroThresh = (dataType) epsilon;
  // if (d1Size + d2Size > largeSize + 1) {
  //   zeroThresh = toSort.at(d1Size + d2Size - largeSize);
  //   if (toSort.at(d1Size + d2Size - (largeSize+1)) == zeroThresh)
  //     zeroThresh += (dataType) epsilon;
  // }
  // if (zeroThresh < epsilon) zeroThresh = epsilon;

  return s;
}

void BottleneckDistance::computeMinMaxSaddleNumberAndMapping(
  const std::vector<diagramTuple> &CTDiagram,
  int dSize,
  int &nbMin,
  int &nbMax,
  int &nbSaddle,
  std::vector<int> &minMap,
  std::vector<int> &maxMap,
  std::vector<int> &sadMap,
  const dataType zeroThresh) {
  for(int i = 0; i < dSize; ++i) {
    const diagramTuple &t = CTDiagram[i];
    BNodeType nt1 = t.birthType;
    BNodeType nt2 = t.deathType;
    dataType dt = t.persistence;
    if(abs<dataType>(dt) < zeroThresh)
      continue;

    if(nt1 == BLocalMin && nt2 == BLocalMax) {
      nbMax++;
      maxMap.push_back(i);
    } else {
      if(nt1 == BLocalMax || nt2 == BLocalMax) {
        nbMax++;
        maxMap.push_back(i);
      }
      if(nt1 == BLocalMin || nt2 == BLocalMin) {
        nbMin++;
        minMap.push_back(i);
      }
      if((nt1 == BSaddle1 && nt2 == BSaddle2)
         || (nt1 == BSaddle2 && nt2 == BSaddle1)) {
        nbSaddle++;
        sadMap.push_back(i);
      }
    }
  }
}

void BottleneckDistance::buildCostMatrices(
  const std::vector<diagramTuple> &CTDiagram1,
  const std::vector<diagramTuple> &CTDiagram2,
  const int d1Size,
  const int d2Size,
  std::function<dataType(const diagramTuple, const diagramTuple)>
    &distanceFunction,
  std::function<dataType(const diagramTuple)> &diagonalDistanceFunction,
  const double zeroThresh,
  std::vector<std::vector<dataType>> &minMatrix,
  std::vector<std::vector<dataType>> &maxMatrix,
  std::vector<std::vector<dataType>> &sadMatrix,
  const bool reverseMin,
  const bool reverseMax,
  const bool reverseSad,
  const int wasserstein) {
  int maxI = 0, minI = 0;
  int maxJ = 0, minJ = 0;
  int sadI = 0, sadJ = 0;

  for(int i = 0; i < d1Size; ++i) {
    const diagramTuple &t1 = CTDiagram1[i];
    if(abs<dataType>(t1.persistence) < zeroThresh)
      continue;

    BNodeType t11 = t1.birthType;
    BNodeType t13 = t1.deathType;
    bool isMin1 = (t11 == BLocalMin || t13 == BLocalMin);
    bool isMax1 = (t11 == BLocalMax || t13 == BLocalMax);
    bool isSad1 = (t11 == BSaddle1 && t13 == BSaddle2)
                  || (t11 == BSaddle2 && t13 == BSaddle1);
    if(t11 == BLocalMin && t13 == BLocalMax) {
      isMin1 = false;
      isMax1 = true;
    }

    minJ = 0;
    maxJ = 0;
    sadJ = 0;

    for(int j = 0; j < d2Size; ++j) {
      const diagramTuple &t2 = CTDiagram2[j];
      if(abs<dataType>(t2.persistence) < zeroThresh)
        continue;

      BNodeType t21 = t2.birthType;
      BNodeType t23 = t2.deathType;
      bool isMin2 = (t21 == BLocalMin || t23 == BLocalMin);
      bool isMax2 = (t21 == BLocalMax || t23 == BLocalMax);
      bool isSad2 = (t21 == BSaddle1 && t23 == BSaddle2)
                    || (t21 == BSaddle2 && t23 == BSaddle1);
      if(t21 == BLocalMin && t23 == BLocalMax) {
        isMin2 = false;
        isMax2 = true;
      }
      if((isMin1 && !isMin2) || (isMax1 && !isMax2) || (isSad1 && !isSad2))
        continue;

      dataType distance = distanceFunction(t1, t2);
      dataType diag1 = diagonalDistanceFunction(t1);
      dataType diag2 = diagonalDistanceFunction(t2);

      if(distance > diag1 + diag2)
        distance = std::numeric_limits<dataType>::max();

      if(isMin1 && isMin2) {
        if(reverseMin)
          minMatrix[minJ++][minI] = distance;
        else
          minMatrix[minI][minJ++] = distance;
      } else if(isMax1 && isMax2) {
        if(reverseMax)
          maxMatrix[maxJ++][maxI] = distance;
        else
          maxMatrix[maxI][maxJ++] = distance;
      } else if(isSad1 && isSad2) {
        if(reverseSad)
          sadMatrix[sadJ++][sadI] = distance;
        else
          sadMatrix[sadI][sadJ++] = distance;
      }
    }

    dataType distanceToDiagonal = diagonalDistanceFunction(t1);
    if(isMin1) {
      if(reverseMin)
        minMatrix[minJ++][minI] = distanceToDiagonal;
      else
        minMatrix[minI][minJ++] = distanceToDiagonal;
    }
    if(isMax1) {
      if(reverseMax)
        maxMatrix[maxJ++][maxI] = distanceToDiagonal;
      else
        maxMatrix[maxI][maxJ++] = distanceToDiagonal;
    }
    if(isSad1) {
      if(reverseSad)
        sadMatrix[sadJ++][sadI] = distanceToDiagonal;
      else
        sadMatrix[sadI][sadJ++] = distanceToDiagonal;
    }

    if(isMin1)
      ++minI;
    if(isMax1)
      ++maxI;
    if(isSad1)
      ++sadI;
  }

  minJ = 0;
  maxJ = 0;
  sadJ = 0;

  // Last row: match remaining J components with diagonal.
  for(int j = 0; j < d2Size; ++j) {
    const diagramTuple &t2 = CTDiagram2[j];
    if(abs<dataType>(t2.persistence) < zeroThresh)
      continue;

    BNodeType t21 = t2.birthType;
    BNodeType t23 = t2.deathType;
    bool isMin2 = (t21 == BLocalMin || t23 == BLocalMin);
    bool isMax2 = (t21 == BLocalMax || t23 == BLocalMax);
    bool isSad2 = (t21 == BSaddle1 && t23 == BSaddle2)
                  || (t21 == BSaddle2 && t23 == BSaddle1);
    if(t21 == BLocalMin && t23 == BLocalMax) {
      isMin2 = false;
      isMax2 = true;
    }

    dataType distanceToDiagonal = diagonalDistanceFunction(t2);
    if(isMin2) {
      if(reverseMin)
        minMatrix[minJ++][minI] = distanceToDiagonal;
      else
        minMatrix[minI][minJ++] = distanceToDiagonal;
    }
    if(isMax2) {
      if(reverseMax)
        maxMatrix[maxJ++][maxI] = distanceToDiagonal;
      else
        maxMatrix[maxI][maxJ++] = distanceToDiagonal;
    }
    if(isSad2) {
      if(reverseSad)
        sadMatrix[sadJ++][sadI] = distanceToDiagonal;
      else
        sadMatrix[sadI][sadJ++] = distanceToDiagonal;
    }
  }

  // Last cell
  {
    if(reverseMin)
      minMatrix[minJ][minI] = std::numeric_limits<dataType>::max();
    else
      minMatrix[minI][minJ] = std::numeric_limits<dataType>::max();
  }
  {
    if(reverseMax)
      maxMatrix[maxJ][maxI] = std::numeric_limits<dataType>::max();
    else
      maxMatrix[maxI][maxJ] = std::numeric_limits<dataType>::max();
  }
  {
    if(reverseSad)
      sadMatrix[sadJ][sadI] = std::numeric_limits<dataType>::max();
    else
      sadMatrix[sadI][sadJ] = std::numeric_limits<dataType>::max();
  }
}

void BottleneckDistance::solvePWasserstein(
  const int nbRow,
  const int nbCol,
  std::vector<std::vector<dataType>> &matrix,
  std::vector<matchingTuple> &matchings,
  Munkres &solver) {
  solver.setInput(nbRow, nbCol, (void *)&matrix);
  solver.run(matchings);
  solver.clearMatrix();
}

void BottleneckDistance::solveInfinityWasserstein(
  const int nbRow,
  const int nbCol,
  const int nbRowToCut,
  const int nbColToCut,
  std::vector<std::vector<dataType>> &matrix,
  std::vector<matchingTuple> &matchings,
  GabowTarjan &solver) {
  std::vector<std::vector<dataType>> bottleneckMatrix(
    nbRow, std::vector<dataType>(nbCol));

  // Copy input matrix.
  for(int i = 0; i < nbRow; ++i)
    for(int j = 0; j < nbCol; ++j)
      bottleneckMatrix[i][j] = matrix[i][j];

  // Solve.
  solver.setInput(nbRow, nbCol, (void *)&bottleneckMatrix);
  solver.run(matchings);
  solver.clear();
}

dataType BottleneckDistance::buildMappings(
  const std::vector<matchingTuple> &inputMatchings,
  const bool transposeGlobal,
  const bool transposeLocal,
  std::vector<matchingTuple> &outputMatchings,
  const std::vector<int> &m1,
  const std::vector<int> &m2,
  int wasserstein) {
  // Input map permutation (so as to ignore transposition later on)
  const std::vector<int> map1 = transposeLocal ? m2 : m1;
  const std::vector<int> map2 = transposeLocal ? m1 : m2;

  std::stringstream msg;
  dataType addedPersistence = 0;
  for(int i = 0, s = (int)inputMatchings.size(); i < s; ++i) {
    matchingTuple t = inputMatchings.at(i);
    dataType val = abs<dataType>(std::get<2>(t));

    int p1 = std::get<0>(t);
    int p2 = std::get<1>(t);

    if(p1 >= (int)map1.size() || p1 < 0) {
      addedPersistence = (wasserstein > 0 ? addedPersistence + val
                                          : std::max(val, addedPersistence));
    } else if(p2 >= (int)map2.size() || p2 < 0) {
      addedPersistence = (wasserstein > 0 ? addedPersistence + val
                                          : std::max(val, addedPersistence));
    } else {
      int point1 = map1.at((unsigned long)p1);
      int point2 = map2.at((unsigned long)p2);
      bool doTranspose = transposeGlobal ^ transposeLocal;

      matchingTuple newT = doTranspose ? std::make_tuple(point2, point1, val)
                                       : std::make_tuple(point1, point2, val);

      outputMatchings.push_back(newT);
    }
  }

  return addedPersistence;
}

} // namespace ttk

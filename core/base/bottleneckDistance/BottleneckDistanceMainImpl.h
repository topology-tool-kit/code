#pragma once

#include "BottleneckDistance.h"

namespace ttk {


//  vector <   -- diagram
//    tuple <    -- pair of critical points
//      int , NodeType
//      int, NodeType
//      dataType  -- persistance of pair
//      int  -- type (0/min, 1/saddle, 2/max)
//      dataType            -- scalar value at vertex 1
//      float, float, float -- vertex 1 coordinates
//      dataType            -- scalar value at vertex 2
//      float, float, float -- vertex 2 coordinates
int BottleneckDistance::computeBottleneck(const std::vector<diagramTuple> &d1,
                                          const std::vector<diagramTuple> &d2,
                                          std::vector<matchingTuple> &matchings,
                                          const bool usePersistenceMetric) {
  auto d1Size = (int)d1.size();
  auto d2Size = (int)d2.size();

  bool transposeOriginal = d1Size > d2Size;
  const std::vector<diagramTuple> &CTDiagram1 = transposeOriginal ? d2 : d1;
  const std::vector<diagramTuple> &CTDiagram2 = transposeOriginal ? d1 : d2;
  if(transposeOriginal) {
    int temp = d1Size;
    d1Size = d2Size;
    d2Size = temp;
  }

  if(transposeOriginal) {
    this->printMsg("The first persistence diagram is larger than the second.");
    this->printMsg("Solving the transposed problem.");
  }

  // Check user parameters.
  const int wasserstein = (wasserstein_ == "inf") ? -1 : stoi(wasserstein_);
  if(wasserstein < 0 && wasserstein != -1)
    return -4;

  // Needed to limit computation time.
  const dataType zeroThresh = this->computeMinimumRelevantPersistence(
    CTDiagram1, CTDiagram2, d1Size, d2Size);

  // Initialize solvers.
  std::vector<matchingTuple> minMatchings;
  std::vector<matchingTuple> maxMatchings;
  std::vector<matchingTuple> sadMatchings;

  // Initialize cost matrices.
  int nbRowMin = 0, nbColMin = 0;
  int maxRowColMin = 0, minRowColMin = 0;
  int nbRowMax = 0, nbColMax = 0;
  int maxRowColMax = 0, minRowColMax = 0;
  int nbRowSad = 0, nbColSad = 0;
  int maxRowColSad = 0, minRowColSad = 0;

  // Remap for matchings.
  std::vector<int> minMap1;
  std::vector<int> minMap2;
  std::vector<int> maxMap1;
  std::vector<int> maxMap2;
  std::vector<int> sadMap1;
  std::vector<int> sadMap2;

  this->computeMinMaxSaddleNumberAndMapping(CTDiagram1, d1Size, nbRowMin,
                                            nbRowMax, nbRowSad, minMap1,
                                            maxMap1, sadMap1, zeroThresh);
  this->computeMinMaxSaddleNumberAndMapping(CTDiagram2, d2Size, nbColMin,
                                            nbColMax, nbColSad, minMap2,
                                            maxMap2, sadMap2, zeroThresh);

  // Automatically transpose if nb rows > nb cols
  maxRowColMin = std::max(nbRowMin + 1, nbColMin + 1);
  maxRowColMax = std::max(nbRowMax + 1, nbColMax + 1);
  maxRowColSad = std::max(nbRowSad + 1, nbColSad + 1);

  minRowColMin = std::min(nbRowMin + 1, nbColMin + 1);
  minRowColMax = std::min(nbRowMax + 1, nbColMax + 1);
  minRowColSad = std::min(nbRowSad + 1, nbColSad + 1);

  std::vector<std::vector<dataType>> minMatrix(
    (unsigned long)minRowColMin, std::vector<dataType>(maxRowColMin));
  std::vector<std::vector<dataType>> maxMatrix(
    (unsigned long)minRowColMax, std::vector<dataType>(maxRowColMax));
  std::vector<std::vector<dataType>> sadMatrix(
    (unsigned long)minRowColSad, std::vector<dataType>(maxRowColSad));

  double px = px_;
  double py = py_;
  double pz = pz_;
  double pe = pe_;
  double ps = ps_;

  std::function<dataType(const diagramTuple, const diagramTuple)>
    distanceFunction
    = [wasserstein, px, py, pz, pe, ps](
        const diagramTuple a, const diagramTuple b) -> dataType {
    BNodeType ta1 = a.birthType;
    BNodeType ta2 = a.deathType;
    const int w = wasserstein > 1 ? wasserstein : 1; // L_inf not managed.

    // We don't match critical points of different index.
    // This must be ensured before calling the distance function.
    // BNodeType tb1 = get<1>(b);
    // BNodeType tb2 = get<3>(b);
    bool isMin1 = ta1 == BLocalMin;
    bool isMax1 = ta2 == BLocalMax;
    // bool isBoth = isMin1 && isMax1;

    dataType rX = a.birthPoint.val;
    dataType rY = a.deathPoint.val;
    dataType cX = b.birthPoint.val;
    dataType cY = b.deathPoint.val;
    dataType x = ((isMin1 && !isMax1) ? pe : ps)
                 * Geometry::pow(abs_diff<dataType>(rX, cX), w);
    dataType y
      = (isMax1 ? pe : ps) * Geometry::pow(abs_diff<dataType>(rY, cY), w);
    double geoDistance
      = isMax1
          ? (px * Geometry::pow(abs(a.deathPoint.x - b.deathPoint.x), w)
             + py * Geometry::pow(abs(a.deathPoint.y - b.deathPoint.y), w)
             + pz * Geometry::pow(abs(a.deathPoint.z - b.deathPoint.z), w))
        : isMin1
          ? (px * Geometry::pow(abs(a.birthPoint.x - b.birthPoint.x), w)
             + py * Geometry::pow(abs(a.birthPoint.y - b.birthPoint.y), w)
             + pz * Geometry::pow(abs(a.birthPoint.z - b.birthPoint.z), w))
          : (px
               * Geometry::pow(abs(a.birthPoint.x + a.deathPoint.x) / 2
                                 - abs(b.birthPoint.x + b.deathPoint.x) / 2,
                               w)
             + py
                 * Geometry::pow(abs(a.birthPoint.y + a.deathPoint.y) / 2
                                   - abs(b.birthPoint.y + b.deathPoint.y) / 2,
                                 w)
             + pz
                 * Geometry::pow(abs(a.birthPoint.z + a.deathPoint.z) / 2
                                   - abs(b.birthPoint.z + b.deathPoint.z) / 2,
                                 w));

    double persDistance = x + y;
    double val = persDistance + geoDistance;
    val = Geometry::pow(val, 1.0 / w);
    return val;
  };

  std::function<dataType(const diagramTuple)> diagonalDistanceFunction
    = [wasserstein, px, py, pz, ps, pe](const diagramTuple a) -> dataType {
    BNodeType ta1 = a.birthType;
    BNodeType ta2 = a.deathType;
    const int w = wasserstein > 1 ? wasserstein : 1;
    bool isMin1 = ta1 == BLocalMin;
    bool isMax1 = ta2 == BLocalMax;

    dataType rX = a.birthPoint.val;
    dataType rY = a.deathPoint.val;
    double x1 = a.birthPoint.x;
    double y1 = a.birthPoint.y;
    double z1 = a.birthPoint.z;
    double x2 = a.deathPoint.x;
    double y2 = a.deathPoint.y;
    double z2 = a.deathPoint.z;

    double infDistance = (isMin1 || isMax1 ? pe : ps)
                         * Geometry::pow(abs_diff<dataType>(rX, rY), w);
    double geoDistance = (px * Geometry::pow(abs(x2 - x1), w)
                          + py * Geometry::pow(abs(y2 - y1), w)
                          + pz * Geometry::pow(abs(z2 - z1), w));
    double val = infDistance + geoDistance;
    return Geometry::pow(val, 1.0 / w);
  };

  const bool transposeMin = nbRowMin > nbColMin;
  const bool transposeMax = nbRowMax > nbColMax;
  const bool transposeSad = nbRowSad > nbColSad;

  Timer t;

  this->buildCostMatrices(
    CTDiagram1, CTDiagram2, d1Size, d2Size, distanceFunction,
    diagonalDistanceFunction, zeroThresh, minMatrix, maxMatrix, sadMatrix,
    transposeMin, transposeMax, transposeSad, wasserstein);

  if(wasserstein > 0) {

    if(nbRowMin > 0 && nbColMin > 0) {
      Munkres solverMin;
      this->printMsg("Affecting minima...");
      this->solvePWasserstein(
        minRowColMin, maxRowColMin, minMatrix, minMatchings, solverMin);
    }

    if(nbRowMax > 0 && nbColMax > 0) {
      Munkres solverMax;
      this->printMsg("Affecting maxima...");
      this->solvePWasserstein(
        minRowColMax, maxRowColMax, maxMatrix, maxMatchings, solverMax);
    }

    if(nbRowSad > 0 && nbColSad > 0) {
      Munkres solverSad;
      this->printMsg("Affecting saddles...");
      this->solvePWasserstein(
        minRowColSad, maxRowColSad, sadMatrix, sadMatchings, solverSad);
    }

  } else {

    // Launch solving for minima.
    if(nbRowMin > 0 && nbColMin > 0) {
      GabowTarjan solverMin;
      this->printMsg("Affecting minima...");
      this->solveInfinityWasserstein(minRowColMin, maxRowColMin, nbRowMin,
                                     nbColMin, minMatrix, minMatchings,
                                     solverMin);
    }

    // Launch solving for maxima.
    if(nbRowMax > 0 && nbColMax > 0) {
      GabowTarjan solverMax;
      this->printMsg("Affecting maxima...");
      this->solveInfinityWasserstein(minRowColMax, maxRowColMax, nbRowMax,
                                     nbColMax, maxMatrix, maxMatchings,
                                     solverMax);
    }

    // Launch solving for saddles.
    if(nbRowSad > 0 && nbColSad > 0) {
      GabowTarjan solverSad;
      this->printMsg("Affecting saddles...");
      this->solveInfinityWasserstein(minRowColSad, maxRowColSad, nbRowSad,
                                     nbColSad, sadMatrix, sadMatchings,
                                     solverSad);
    }
  }

  this->printMsg("TTK CORE DONE", 1, t.getElapsedTime());

  // Rebuild mappings.
  // Begin cost computation for unpaired vertices.
  // std::cout << "Min" << std::endl;
  dataType addedMinPersistence = this->buildMappings(
    minMatchings, transposeOriginal, transposeMin, matchings, minMap1, minMap2,
    wasserstein);

  // std::cout << "Max" << std::endl;
  dataType addedMaxPersistence = this->buildMappings(
    maxMatchings, transposeOriginal, transposeMax, matchings, maxMap1, maxMap2,
    wasserstein);

  // std::cout << "Sad" << std::endl;
  dataType addedSadPersistence = this->buildMappings(
    sadMatchings, transposeOriginal, transposeSad, matchings, sadMap1, sadMap2,
    wasserstein);

  // TODO [HIGH] do that for embeddings
  // Recompute matching weights for user-friendly distance.
  dataType d = 0;
  std::vector<bool> paired1(d1Size);
  std::vector<bool> paired2(d2Size);
  for(int b = 0; b < d1Size; ++b)
    paired1[b] = false;
  for(int b = 0; b < d2Size; ++b)
    paired2[b] = false;

  int numberOfMismatches = 0;
  for(int m = 0, ms = (int)matchings.size(); m < ms; ++m) {
    matchingTuple mt = matchings[m];
    int i = transposeOriginal ? std::get<1>(mt) : std::get<0>(mt);
    int j = transposeOriginal ? std::get<0>(mt) : std::get<1>(mt);
    // dataType val = std::get<2>(t);

    diagramTuple t1 = CTDiagram1[i];
    diagramTuple t2 = CTDiagram2[j];
    // dataType rX = t1.birthPoint.val; dataType rY = t1.deathPoint.val;
    // dataType cX = t2.birthPoint.val; dataType cY = t2.deathPoint.val;
    // dataType x = rX - cX; dataType y = rY - cY;
    paired1[i] = true;
    paired2[j] = true;
    // dataType lInf = std::max(abs<dataType>(x), abs<dataType>(y));

    // if (((wasserstein < 0 && lInf != val) || (wasserstein > 0 && pow(lInf,
    // wasserstein) != val)))
    //++numberOfMismatches;

    dataType partialDistance = distanceFunction(t1, t2);
    // wasserstein > 0 ? pow(lInf, wasserstein) : std::max(d, lInf);

    if(wasserstein > 0)
      d += partialDistance;
    else
      d = partialDistance;
  }

  if(numberOfMismatches > 0) {
    this->printWrn("Distance mismatch when rebuilding "
                   + std::to_string(numberOfMismatches) + " matchings");
  }

  dataType affectationD = d;
  d = wasserstein > 0
        ? Geometry::pow(
          d + addedMaxPersistence + addedMinPersistence + addedSadPersistence,
          (1.0 / (double)wasserstein))
        : std::max(
          d, std::max(addedMaxPersistence,
                      std::max(addedMinPersistence, addedSadPersistence)));

  {
    std::stringstream msg;
    this->printMsg("Computed distance:");
    this->printMsg("diagMax(" + std::to_string(addedMaxPersistence)
                   + "), diagMin(" + std::to_string(addedMinPersistence)
                   + "), diagSad(" + std::to_string(addedSadPersistence) + ")");
    this->printMsg("affAll(" + std::to_string(affectationD) + "), res("
                   + std::to_string(d) + ")");
  }

  distance_ = (double)d;
  return 0;
}

} // namespace ttk

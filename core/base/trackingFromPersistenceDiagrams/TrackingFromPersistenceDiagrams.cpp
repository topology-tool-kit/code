#include <TrackingFromPersistenceDiagrams.h>

ttk::TrackingFromPersistenceDiagrams::TrackingFromPersistenceDiagrams() {
  inputData_ = nullptr;
  numberOfInputs_ = 0;
  this->setDebugMsgPrefix("TrackingFromPersistenceDiagrams");
}

ttk::TrackingFromPersistenceDiagrams::~TrackingFromPersistenceDiagrams() {
  if(inputData_)
    free(inputData_);
}

int ttk::TrackingFromPersistenceDiagrams::execute() {
  ttk::Timer t;

  // Check the consistency of the variables
#ifndef TTK_ENABLE_KAMIKAZE
  if(!numberOfInputs_)
    return -1;
  if(!inputData_)
    return -3;

  for(int i = 0; i < numberOfInputs_; i++) {
    if(!inputData_[i])
      return -4;
  }
#endif

  this->printMsg("Complete", 1, t.getElapsedTime(), threadNumber_);

  return 0;
}

int ttk::TrackingFromPersistenceDiagrams::performSingleMatching(
  int i,
  std::vector<std::vector<diagramTuple>> &inputPersistenceDiagrams,
  std::vector<std::vector<matchingTuple>> &outputMatchings,
  std::string algorithm,
  std::string wasserstein,
  double tolerance,
  bool is3D,
  double alpha,
  double px,
  double py,
  double pz,
  double ps,
  double pe) {
  ttk::BottleneckDistance bottleneckDistance_;
  // bottleneckDistance_.setWrapper(this);
  bottleneckDistance_.setPersistencePercentThreshold(tolerance);
  bottleneckDistance_.setPX(px);
  bottleneckDistance_.setPY(py);
  bottleneckDistance_.setPZ(pz);
  bottleneckDistance_.setPS(ps);
  bottleneckDistance_.setPE(pe);
  bottleneckDistance_.setAlgorithm(algorithm);
  bottleneckDistance_.setWasserstein(wasserstein);

  bottleneckDistance_.setCTDiagram1(&inputPersistenceDiagrams[i]);
  bottleneckDistance_.setCTDiagram2(&inputPersistenceDiagrams[i + 1]);
  bottleneckDistance_.setOutputMatchings(&outputMatchings[i]);
  bottleneckDistance_.execute(false);

  return 0;
}

int ttk::TrackingFromPersistenceDiagrams::performMatchings(
  int numInputs,
  std::vector<std::vector<diagramTuple>> &inputPersistenceDiagrams,
  std::vector<std::vector<matchingTuple>> &outputMatchings,
  const std::string &algorithm,
  const std::string &wasserstein,
  double tolerance,
  bool is3D,
  double alpha,
  double px,
  double py,
  double pz,
  double ps,
  double pe) {

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
  for(int i = 0; i < numInputs - 1; ++i) {
    performSingleMatching(
      i, inputPersistenceDiagrams, outputMatchings,
      algorithm, // Not from paraview, from enclosing tracking plugin
      wasserstein, tolerance, is3D,
      alpha, // Blending
      px, py, pz, ps, pe // Coefficients
    );
  }

  // Never from PV,
  // Alaways activate output matchings.
  return 0;
}

int ttk::TrackingFromPersistenceDiagrams::performTracking(
  std::vector<std::vector<diagramTuple>> &allDiagrams,
  std::vector<std::vector<matchingTuple>> &allMatchings,
  std::vector<trackingTuple> &trackings) {
  auto numPersistenceDiagramsInput = (int)allDiagrams.size();

  for(int in = 1; in < numPersistenceDiagramsInput - 1; ++in) {
    std::vector<matchingTuple> matchings1 = allMatchings[in - 1];
    std::vector<matchingTuple> matchings2 = allMatchings[in];

    auto matchingsSize1 = (int)matchings1.size();
    auto matchingsSize2 = (int)matchings2.size();
    int endIndex = numPersistenceDiagramsInput - 2;

    for(int i = 0; i < matchingsSize1; ++i) {
      auto m1ai0 = (int)std::get<0>(matchings1[i]);
      auto m1ai1 = (int)std::get<1>(matchings1[i]);

      for(int j = 0; j < matchingsSize2; ++j) {
        auto m2aj0 = (int)std::get<0>(matchings2[j]);
        auto m2aj1 = (int)std::get<1>(matchings2[j]);

        if(m1ai1 != m2aj0)
          continue;

        // Detect in trackings and push.
        bool found = false;
        for(trackingTuple &tt : trackings) {
          int chainStart = std::get<0>(tt);
          int chainEnd = std::get<1>(tt);
          std::vector<BIdVertex> &chain = std::get<2>(tt);

          if(chainEnd == -1) {
            auto chainSize = (int)chain.size();
            if(chainSize == 0) {
              // Should not happen
              this->printErr("Brain error.");
            } else if(chainStart + chainSize == in
                      && chain.at((unsigned long)chainSize - 1) == m1ai0) {
              found = true;
              chain.push_back(m1ai1);
              int numEnd = in == endIndex ? endIndex : -1;
              if(in == endIndex) {
                chain.push_back(m2aj1);
                std::get<1>(tt) = numEnd;
              }
              std::get<2>(tt) = chain;
            }
          }
          tt = std::make_tuple(chainStart, chainEnd, chain);
        }

        if(!found) {
          std::vector<BIdVertex> chain;
          chain.push_back(m1ai0);
          chain.push_back(m1ai1);
          if(in == endIndex) {
            chain.push_back(m2aj1);
          }
          int numEnd = in == endIndex ? endIndex : -1;
          trackingTuple tt = std::make_tuple(in - 1, numEnd, chain);
          trackings.push_back(tt);
        }
        // Create new.
      }
    }

    // End non-matched chains.
    for(trackingTuple &tt : trackings) {
      int chainStart = std::get<0>(tt);
      int chainEnd = std::get<1>(tt);
      if(chainEnd == -1) {
        std::vector<BIdVertex> &chain = std::get<2>(tt);
        auto chainSize = (int)chain.size();
        if(chainStart + chainSize - 1 < in)
          std::get<1>(tt) = in - 1;
      }
    }
  }

  // Post-processing
  std::sort(trackings.begin(), trackings.end(),
            [](const trackingTuple &a, const trackingTuple &b) -> bool {
              return std::get<0>(a) < std::get<0>(b);
            });

  return 0;
}

int ttk::TrackingFromPersistenceDiagrams::performPostProcess(
  std::vector<std::vector<diagramTuple>> &allDiagrams,
  std::vector<trackingTuple> &trackings,
  std::vector<std::set<int>> &trackingTupleToMerged,
  double postProcThresh) {
  auto numPersistenceDiagramsInput = (int)allDiagrams.size();

  // Merge close connected components with threshold.
  for(unsigned int k = 0; k < trackings.size(); ++k) {
    trackingTuple tk = trackings[k];
    int startK = std::get<0>(tk);
    int endK = std::get<1>(tk);
    if(endK < 0)
      endK = numPersistenceDiagramsInput - 1;
    std::vector<BIdVertex> chainK = std::get<2>(tk);
    std::vector<diagramTuple> &diagramStartK = allDiagrams[startK];
    std::vector<diagramTuple> &diagramEndK = allDiagrams[endK];

    auto n1 = (int)chainK.at(0);
    auto n2 = (int)chainK.at(chainK.size() - 1);
    diagramTuple &tuple1 = diagramStartK[n1];
    diagramTuple &tuple2 = diagramEndK[n2];

    double x1, y1, z1, x2, y2, z2;

    BNodeType point1Type1 = tuple1.birthType;
    BNodeType point1Type2 = tuple1.deathType;
    bool t11Min = point1Type1 == BLocalMin;
    bool t11Max = point1Type1 == BLocalMax;
    bool t12Min = point1Type2 == BLocalMin;
    bool t12Max = point1Type2 == BLocalMax;
    // bool bothEx1 = t11Ex && t12Ex;
    bool t1Max = t11Max || t12Max;
    bool t1Min = !t1Max && (t11Min || t12Min);

    x1 = t1Max ? tuple1.deathPoint.x : t1Min ? tuple1.birthPoint.x : 0;
    y1 = t1Max ? tuple1.deathPoint.y : t1Min ? tuple1.birthPoint.y : 0;
    z1 = t1Max ? tuple1.deathPoint.z : t1Min ? tuple1.birthPoint.z : 0;

    BNodeType point2Type1 = tuple2.birthType;
    BNodeType point2Type2 = tuple2.deathType;
    bool t21Min = point2Type1 == BLocalMin;
    bool t21Max = point2Type1 == BLocalMax;
    bool t22Min = point2Type2 == BLocalMin;
    bool t22Max = point2Type2 == BLocalMax;
    // bool bothEx2 = t21Ex && t22Ex;
    bool t2Max = t21Max || t22Max;
    bool t2Min = !t2Max && (t21Min || t22Min);

    // if (bothEx2) {
    x2 = t2Max ? tuple2.deathPoint.x : t2Min ? tuple2.birthPoint.x : 0;
    y2 = t2Max ? tuple2.deathPoint.y : t2Min ? tuple2.birthPoint.y : 0;
    z2 = t2Max ? tuple2.deathPoint.z : t2Min ? tuple2.birthPoint.z : 0;
    // }
    // if (!bothEx1 && !bothEx2)
    //  continue;

    // Saddle-saddle matching not supported.
    if(!t1Min && !t2Min && !t1Max && !t2Max)
      continue;

    // Check every other tracking trajectory.
    for(unsigned int m = k + 1; m < trackings.size(); ++m) {
      trackingTuple &tm = trackings[m];
      int startM = std::get<0>(tm);
      int endM = std::get<1>(tm);
      std::vector<BIdVertex> &chainM = std::get<2>(tm);
      if((endK > 0 && startM > endK) || (endM > 0 && startK > endM))
        continue;

      for(int c = 0; c < (int)chainM.size(); ++c) {
        bool doMatch1 = startM + c == startK;
        bool doMatch2 = startM + c == endK;

        // if (startM + c != startK && startM + c != endK) continue;
        if(!doMatch1 && !doMatch2)
          continue;

        /// Check proximity.
        auto n3 = (int)chainM[c];
        std::vector<diagramTuple> &diagramM = allDiagrams[startM + c];
        diagramTuple &tuple3 = diagramM[n3];
        double x3, y3, z3;
        BNodeType point3Type1 = tuple3.birthType;
        BNodeType point3Type2 = tuple3.deathType;
        bool t31Min = point3Type1 == BLocalMin;
        bool t31Max = point3Type1 == BLocalMax;
        bool t32Min = point3Type2 == BLocalMin;
        bool t32Max = point3Type2 == BLocalMax;
        // bool bothEx3 = t31Ex && t32Ex;
        // if (!bothEx3)
        //  continue;
        bool t3Max = t31Max || t32Max;
        bool t3Min = !t3Max && (t31Min || t32Min);

        x3 = t3Max ? tuple3.deathPoint.x : t3Min ? tuple3.birthPoint.x : 0;
        y3 = t3Max ? tuple3.deathPoint.y : t3Min ? tuple3.birthPoint.y : 0;
        z3 = t3Max ? tuple3.deathPoint.z : t3Min ? tuple3.birthPoint.z : 0;

        double dist = 0;
        bool hasMatched = false;
        if(doMatch1 && ((t3Max && t1Max) || (t3Min && t1Min))) {
          double dist13
            = sqrt(Geometry::pow(x1 - x3, 2) + Geometry::pow(y1 - y3, 2)
                   + Geometry::pow(z1 - z3, 2));
          dist = dist13;
          if(dist13 >= postProcThresh)
            continue;
          hasMatched = true;
        }

        if(doMatch2 && ((t3Max && t2Max) || (t3Min && t2Min))) {
          double dist23
            = sqrt(Geometry::pow(x2 - x3, 2) + Geometry::pow(y2 - y3, 2)
                   + Geometry::pow(z2 - z3, 2));
          dist = dist23;
          if(dist23 >= postProcThresh)
            continue;
          hasMatched = true;
        }

        if(!hasMatched)
          continue;

        /// Merge!
        std::stringstream msg;
        msg << "Merged " << m << " with " << k << ": d = " << dist << ".";
        printMsg(msg.str());

        // Get every other tracking trajectory.
        std::set<int> &mergedM = trackingTupleToMerged[m];
        // std::set<int> mergedK = trackingTupleToMerged[k];

        // Push for others to merge.
        // for (auto& i : mergedM) mergedK.insert(i);
        // for (auto& i : mergedK) mergedM.insert(i);
        // mergedK.insert(m);
        mergedM.insert(k);
        break;
      }
    }
  }

  return 0;
}

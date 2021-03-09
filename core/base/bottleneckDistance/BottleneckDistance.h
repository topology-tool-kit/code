/// \ingroup base
/// \class ttk::BottleneckDistance
/// \author Maxime Soler <soler.maxime@total.com>
#pragma once

#ifndef BNodeType
#define BNodeType ttk::CriticalType
#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2
#define BIdVertex int
#endif

// base code includes
#include <GabowTarjan.h>
#include <Munkres.h>
#include <Triangulation.h>

#include <functional>
#include <string>
#include <tuple>

#include <PersistenceDiagram.h>

typedef double dataType;
typedef std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, dataType, int, dataType, float, float, float, dataType, float, float, float> diagramTupleOld;

namespace ttk {
  struct DecoratedDiagramTuple : public ttk::PersistencePair
  {
    DecoratedDiagramTuple(diagramTupleOld t) :
      ttk::PersistencePair(std::get<0>(t), std::get<1>(t), std::get<2>(t), std::get<3>(t), std::get<4>(t), std::get<5>(t)), birthPoint(std::get<6>(t),std::get<7>(t),std::get<8>(t), std::get<9>(t)), deathPoint(std::get<10>(t),std::get<11>(t),std::get<12>(t), std::get<13>(t))
    {}
    DecoratedDiagramTuple() : DecoratedDiagramTuple(diagramTupleOld()) {} 
    struct Coord { float x, y, z; };
    struct DataPoint : public Coord {
      DataPoint(dataType _v, float _x, float _y, float _z) : Coord({_x,_y,_z}),val(_v) {}
      DataPoint() = default;
      dataType val;
    };
    DataPoint birthPoint, deathPoint;
  };
  
} // namespace ttk

using diagramTuple = ttk::DecoratedDiagramTuple;

namespace ttk {

  using trackingTuple = std::tuple<int, int, std::vector<int>>;

  class BottleneckDistance : virtual public Debug {

  public:
    BottleneckDistance()
      : distance_(-1), wasserstein_("inf"), pvAlgorithm_(-1), zeroThreshold_(0),
        px_(0), py_(0), pz_(0), pe_(0), ps_(0) {
      this->setDebugMsgPrefix("BottleneckDistance");
    }

    ~BottleneckDistance(){};

    int execute(bool usePersistenceMetric);

    inline int setPersistencePercentThreshold(double t) {
      zeroThreshold_ = t;
      return 0;
    }

    inline int setPX(double px) {
      px_ = px;
      return 0;
    }
    inline int setPY(double py) {
      py_ = py;
      return 0;
    }
    inline int setPZ(double pz) {
      pz_ = pz;
      return 0;
    }
    inline int setPE(double pe) {
      pe_ = pe;
      return 0;
    }
    inline int setPS(double ps) {
      ps_ = ps;
      return 0;
    }

    inline int setCTDiagram1(void *diagram) {
      outputCT1_ = diagram;
      return 0;
    }

    inline int setCTDiagram2(void *diagram) {
      outputCT2_ = diagram;
      return 0;
    }

    inline int setOutputMatchings(void *matchings) {
      matchings_ = matchings;
      return 0;
    }

    inline int setAlgorithm(const std::string &algorithm) {
      algorithm_ = algorithm;
      return 0;
    }

    inline int setPVAlgorithm(const int algorithm) {
      pvAlgorithm_ = algorithm;
      return 0;
    }

    inline int setWasserstein(const std::string &wasserstein) {
      wasserstein_ = wasserstein;
      return 0;
    }

    inline void message(const char *s) {
      std::stringstream msg;
      msg << s;
      this->printMsg(msg.str());
    }

    double getDistance() {
      return distance_;
    }

    template <typename type>
    static type abs(const type var) {
      return (var >= 0) ? var : -var;
    }

    template <typename type>
    static type abs_diff(const type var1, const type var2) {
      return (var1 > var2) ? var1 - var2 : var2 - var1;
    }

  protected:
    void *outputCT1_;
    void *outputCT2_;
    void *matchings_; // ids from CT1 to CT2
    double distance_;

    std::string wasserstein_;
    std::string algorithm_;
    int pvAlgorithm_;
    double zeroThreshold_;
    double px_;
    double py_;
    double pz_;
    double pe_;
    double ps_;

  private:
    int computeBottleneck(const std::vector<diagramTuple> &d1,
                          const std::vector<diagramTuple> &d2,
                          std::vector<matchingTuple> &matchings,
                          bool usePersistenceMetric);

    double computeGeometricalRange(const std::vector<diagramTuple> &CTDiagram1,
                                   const std::vector<diagramTuple> &CTDiagram2,
                                   int d1Size,
                                   int d2Size) const;

    double computeMinimumRelevantPersistence(
      const std::vector<diagramTuple> &CTDiagram1,
      const std::vector<diagramTuple> &CTDiagram2,
      int d1Size,
      int d2Size) const;

    void computeMinMaxSaddleNumberAndMapping(
      const std::vector<diagramTuple> &CTDiagram,
      int dSize,
      int &nbMin,
      int &nbMax,
      int &nbSaddle,
      std::vector<int> &minMap,
      std::vector<int> &maxMap,
      std::vector<int> &sadMap,
      dataType zeroThresh);

    void buildCostMatrices(
      const std::vector<diagramTuple> &CTDiagram1,
      const std::vector<diagramTuple> &CTDiagram2,
      int d1Size,
      int d2Size,
      std::function<dataType(const diagramTuple, const diagramTuple)>
        &distanceFunction,
      std::function<dataType(const diagramTuple)> &diagonalDistanceFunction,
      double zeroThresh,
      std::vector<std::vector<dataType>> &minMatrix,
      std::vector<std::vector<dataType>> &maxMatrix,
      std::vector<std::vector<dataType>> &sadMatrix,
      bool reverseMin,
      bool reverseMax,
      bool reverseSad,
      int wasserstein);

    void solvePWasserstein(int nbRow,
                           int nbCol,
                           std::vector<std::vector<dataType>> &matrix,
                           std::vector<matchingTuple> &matchings,
                           Munkres &solver);

    void solveInfinityWasserstein(int nbRow,
                                  int nbCol,
                                  int nbRowToCut,
                                  int nbColToCut,
                                  std::vector<std::vector<dataType>> &matrix,
                                  std::vector<matchingTuple> &matchings,
                                  GabowTarjan &solver);

    dataType buildMappings(const std::vector<matchingTuple> &inputMatchings,
                           bool transposeGlobal,
                           bool transposeLocal,
                           std::vector<matchingTuple> &outputMatchings,
                           const std::vector<int> &m1,
                           const std::vector<int> &m2,
                           int wasserstein);
  };

} // namespace ttk


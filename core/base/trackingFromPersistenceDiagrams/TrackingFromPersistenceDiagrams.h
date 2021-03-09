/// \ingroup base
/// \class ttk::TrackingFromPersistenceDiagrams
/// \author Maxime Soler <soler.maxime@total.com>
/// \date August 2018.
#pragma once

// base code includes
#include <BottleneckDistance.h>
#include <set>

namespace ttk {

  class TrackingFromPersistenceDiagrams : virtual public Debug {

  public:
    TrackingFromPersistenceDiagrams();

    ~TrackingFromPersistenceDiagrams();

    /// Execute the package.
    /// \return Returns 0 upon success, negative values otherwise.
    int execute();

    int performSingleMatching(
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
      double pe);

    int performMatchings(
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
      double pe);

    int performTracking(std::vector<std::vector<diagramTuple>> &allDiagrams,
                        std::vector<std::vector<matchingTuple>> &allMatchings,
                        std::vector<trackingTuple> &trackings);

    int performPostProcess(std::vector<std::vector<diagramTuple>> &allDiagrams,
                           std::vector<trackingTuple> &trackings,
                           std::vector<std::set<int>> &trackingTupleToMerged,
                           double postProcThresh);

    /// Pass a pointer to an input array representing a scalarfield.
    /// The array is expected to be correctly allocated. idx in
    /// [0,numberOfInputs_[ \param idx Index of the input scalar field. \param
    /// data Pointer to the data array. \return Returns 0 upon success, negative
    /// values otherwise. \sa setNumberOfInputs() and setVertexNumber().
    inline int setInputDataPointer(int idx, void *data) {
      if(idx < numberOfInputs_)
        inputData_[idx] = data;
      else
        return -1;
      return 0;
    }

    /// Set the number of input scalar fields
    /// \param numberOfInputs Number of input scalar fields.
    /// \return Returns 0 upon success, negative values otherwise
    inline int setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
      return 0;
    }

  protected:
    int numberOfInputs_;
    void **inputData_;
  };
} // namespace ttk


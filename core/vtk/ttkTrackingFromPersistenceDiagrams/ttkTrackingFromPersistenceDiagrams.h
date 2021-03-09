#pragma once

#include <tuple>

#include <TrackingFromPersistenceDiagrams.h>
#include <ttkAlgorithm.h>
// #include <ttkUtils.h>
#include <ttkMacros.h>

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIndent.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkTrackingFromPersistenceDiagramsModule.h>

class TTKTRACKINGFROMPERSISTENCEDIAGRAMS_EXPORT
  ttkTrackingFromPersistenceDiagrams
  : public ttkAlgorithm,
    protected ttk::TrackingFromPersistenceDiagrams {

public:
  static ttkTrackingFromPersistenceDiagrams *New();

  vtkTypeMacro(ttkTrackingFromPersistenceDiagrams, ttkAlgorithm);

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);

  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);

  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);

  vtkSetMacro(PE, double);
  vtkGetMacro(PE, double);

  vtkSetMacro(PS, double);
  vtkGetMacro(PS, double);

  vtkSetMacro(WassersteinMetric, std::string);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, std::string);
  vtkGetMacro(DistanceAlgorithm, std::string);

  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);

  vtkSetMacro(UseGeometricSpacing, bool);
  vtkGetMacro(UseGeometricSpacing, bool);

  vtkSetMacro(Spacing, double);
  vtkGetMacro(Spacing, double);

  vtkSetMacro(DoPostProc, bool);
  vtkGetMacro(DoPostProc, bool);

  vtkSetMacro(PostProcThresh, double);
  vtkGetMacro(PostProcThresh, double);

  using trackingTuple = ttk::trackingTuple;

  static int
    buildMesh(std::vector<trackingTuple> &trackings,
              std::vector<std::vector<matchingTuple>> &outputMatchings,
              std::vector<std::vector<diagramTuple>> &inputPersistenceDiagrams,
              bool useGeometricSpacing,
              double spacing,
              bool DoPostProc,
              std::vector<std::set<int>> &trackingTupleToMerged,
              vtkPoints *points,
              vtkUnstructuredGrid *persistenceDiagram,
              vtkDoubleArray *costScalars,
              vtkDoubleArray *persistenceScalars,
              vtkDoubleArray *valueScalars,
              vtkIntArray *matchingIdScalars,
              vtkIntArray *lengthScalars,
              vtkIntArray *timeScalars,
              vtkIntArray *componentIds,
              vtkIntArray *pointTypeScalars);

protected:
  ttkTrackingFromPersistenceDiagrams();

  ~ttkTrackingFromPersistenceDiagrams() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  // Warn: this is a duplicate from ttkBottleneckDistance.h
  int getPersistenceDiagram(
    std::vector<diagramTuple> &diagram,
    const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram_,
    double spacing,
    int diagramNumber);

  // Warn: ditto
  int augmentPersistenceDiagrams(
    const std::vector<diagramTuple> &diagram1,
    const std::vector<diagramTuple> &diagram2,
    const std::vector<matchingTuple> &matchings,
    vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram1_,
    vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram2_);

private:
  // Input bottleneck config.
  bool UseGeometricSpacing{false};
  bool Is3D{true};
  bool DoPostProc{false};
  double PostProcThresh{0.0};
  double Spacing{1.0};
  double Alpha{1.0};
  double Tolerance{1.0};
  double PX{1};
  double PY{1};
  double PZ{1};
  double PE{1};
  double PS{1};
  std::string DistanceAlgorithm{"ttk"};
  int PVAlgorithm{-1};
  std::string WassersteinMetric{"1"};

  vtkUnstructuredGrid *outputMesh_{nullptr};
};

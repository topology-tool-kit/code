/// \ingroup vtk
/// \class ttkBottleneckDistance
/// \author Maxime Soler <soler.maxime@total.com>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the bottleneckDistance processing package.
///
/// VTK wrapping code for the @BottleneckDistance package.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the corresponding ParaView state file example for a usage example
/// within a VTK pipeline.
///
/// \sa ttk::BottleneckDistance
#pragma once

// TTK includes
#include <BottleneckDistance.h>
#include <ttkAlgorithm.h>
#include <ttkBottleneckDistanceModule.h>

class vtkUnstructuredGrid;

class TTKBOTTLENECKDISTANCE_EXPORT ttkBottleneckDistance
  : public ttkAlgorithm,
    protected ttk::BottleneckDistance {

public:
  static ttkBottleneckDistance *New();

  vtkTypeMacro(ttkBottleneckDistance, ttkAlgorithm);

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  vtkSetMacro(Tolerance, double);
  vtkGetMacro(Tolerance, double);

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

  vtkSetMacro(UseOutputMatching, bool);
  vtkGetMacro(UseOutputMatching, bool);

  vtkSetMacro(BenchmarkSize, int);
  vtkGetMacro(BenchmarkSize, int);

  vtkSetMacro(UsePersistenceMetric, bool);
  vtkGetMacro(UsePersistenceMetric, bool);

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

  double Getresult() {
    return this->distance_;
  }

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation * /*request*/,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  // Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
  int getPersistenceDiagram(std::vector<diagramTuple> &diagram,
                            vtkUnstructuredGrid *const CTPersistenceDiagram_,
                            double spacing,
                            int diagramNumber);

  int generatePersistenceDiagram(std::vector<diagramTuple> &diagram, int size);

  // Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
  int augmentPersistenceDiagrams(
    const std::vector<diagramTuple> &diagram1,
    const std::vector<diagramTuple> &diagram2,
    const std::vector<matchingTuple> &matchings,
    vtkUnstructuredGrid *const CTPersistenceDiagram1_,
    vtkUnstructuredGrid *const CTPersistenceDiagram2_);

  int translateSecondDiagram(vtkUnstructuredGrid *outputCT2, double &spacing);

  int getMatchingMesh(vtkUnstructuredGrid *const outputCT3,
                      const std::vector<diagramTuple> &diagram1,
                      const std::vector<diagramTuple> &diagram2,
                      const std::vector<matchingTuple> &matchings,
                      bool useGeometricSpacing,
                      double spacing,
                      bool is2D);

  int doBenchmark();

protected:
  ttkBottleneckDistance();

private:
  int BenchmarkSize{-1};

  bool UseOutputMatching{false};
  bool Is3D{false};
  double Spacing{0.0};
  double Alpha{1.0};
  double Tolerance{1.0};
  double PX{0.0};
  double PY{0.0};
  double PZ{0.0};
  double PE{1.0};
  double PS{1.0};
  std::string DistanceAlgorithm{};
  std::string WassersteinMetric{"2"};
  bool UsePersistenceMetric{false};
  bool UseGeometricSpacing{false};
  int PVAlgorithm{-1};
};

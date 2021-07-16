/// \ingroup vtk
/// \class ttkPlanarGraphLayout
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.12.2018
///
/// \brief TTK VTK-filter that computes a planar graph layout.
///
/// VTK wrapping code for the @PlanarGraphLayout package.
///
/// This filter computes a planar graph layout of a \b vtkUnstructuredGrid. To
/// improve the quality of the layout it is possible to pass additional field
/// data to the algorithm:\n \b 1) \b Sequences: Points are positioned along the
/// x-axis based on a sequence (e.g., time indicies or scalar values). \b 1) \b
/// Sizes: Points cover space on the y-axis based on their size. \b 1) \b
/// Branches: Points with the same branch label are positioned on straight
/// lines. \b 1) \b Levels: The layout of points with the same level label are
/// computed individually and afterwards nested based on the level hierarchy.
/// This makes it possible to draw nested graphs where each level is a layer of
/// the resulting graph.
///
/// \b Related \b publication: \n
/// 'Nested Tracking Graphs'.
/// Jonas Lukasczyk, Gunther Weber, Ross Maciejewski, Christoph Garth, and Heike
/// Leitte. Computer Graphics Forum (Special Issue, Proceedings Eurographics /
/// IEEE Symposium on Visualization). Vol. 36. No. 3. 2017.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \param Input A \b vtkUnstructuredGrid that represents a graph.
/// \param Output The input \b vtkUnstructuredGrid with an additional point data
/// field that records the computed layout. Note: to project the graph based on
/// the comptued layout use either the \b ttkProjectFromField filter or the \b
/// vtkCalculator.
///
/// \sa ttk::PlanarGraphLayout

#pragma once

// Module
#include <ttkPlanarGraphLayoutModule.h>

// VTK includes
#include <ttkAlgorithm.h>

// TTK includes
#include <PlanarGraphLayout.h>
#include <vtkUnstructuredGrid.h>

class TTKPLANARGRAPHLAYOUT_EXPORT ttkPlanarGraphLayout
  : public ttkAlgorithm,
    protected ttk::PlanarGraphLayout {

private:
  // --- Graph Planar Layout
  // optional field data
  bool UseSequences{false};
  bool UseSizes{false};
  bool UseBranches{false};
  bool UseLevels{false};

  // output field name
  std::string OutputArrayName{"Layout"};

  // --- Merge Tree Planar Layout
  bool InputIsAMergeTree = false;
  bool BranchDecompositionPlanarLayout = false;
  double BranchSpacing = 1.;
  double ImportantPairs = 10.; // important pairs threshold
  double ImportantPairsSpacing = 1.;
  double NonImportantPairsSpacing = 0.1;
  double NonImportantPairsProximity = 0.05;

public:
  // --- Graph Planar Layout
  // getters and setters for optional arrays
  vtkSetMacro(UseSequences, bool);
  vtkGetMacro(UseSequences, bool);

  vtkSetMacro(UseSizes, bool);
  vtkGetMacro(UseSizes, bool);

  vtkSetMacro(UseBranches, bool);
  vtkGetMacro(UseBranches, bool);

  vtkSetMacro(UseLevels, bool);
  vtkGetMacro(UseLevels, bool);

  // getters and setters for output array name
  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  // --- Merge Tree Planar Layout
  vtkSetMacro(InputIsAMergeTree, bool);
  vtkGetMacro(InputIsAMergeTree, bool);

  vtkSetMacro(BranchDecompositionPlanarLayout, bool);
  vtkGetMacro(BranchDecompositionPlanarLayout, bool);

  vtkSetMacro(BranchSpacing, double);
  vtkGetMacro(BranchSpacing, double);

  vtkSetMacro(ImportantPairs, double);
  vtkGetMacro(ImportantPairs, double);

  vtkSetMacro(ImportantPairsSpacing, double);
  vtkGetMacro(ImportantPairsSpacing, double);

  vtkSetMacro(NonImportantPairsSpacing, double);
  vtkGetMacro(NonImportantPairsSpacing, double);

  vtkSetMacro(NonImportantPairsProximity, double);
  vtkGetMacro(NonImportantPairsProximity, double);

  // ---
  static ttkPlanarGraphLayout *New();
  vtkTypeMacro(ttkPlanarGraphLayout, ttkAlgorithm);

protected:
  ttkPlanarGraphLayout();
  ~ttkPlanarGraphLayout();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int planarGraphLayoutCall(vtkInformation *request,
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector);
  int mergeTreePlanarLayoutCall(vtkInformation *request,
                                vtkInformationVector **inputVector,
                                vtkInformationVector *outputVector);

  template <class dataType>
  int mergeTreePlanarLayoutCallTemplate(vtkUnstructuredGrid *treeNodes,
                                        vtkUnstructuredGrid *treeArcs,
                                        vtkUnstructuredGrid *output);
};

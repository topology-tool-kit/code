#include                  "ttkBottleneckDistance.h"

vtkStandardNewMacro(ttkBottleneckDistance)

int ttkBottleneckDistance::doIt(
  vector<vtkDataSet *> &inputs, 
  vector<vtkDataSet *> &outputs)
{

  // Prepare IO
  vtkDataSet *input1 = inputs[0];
  vtkDataSet *input2 = inputs[1];
	vtkDataSet *output1 = outputs[0];
  vtkDataSet *output2 = outputs[1];
  vtkDataSet *output3 = outputs[2];

  vtkUnstructuredGrid *outputCT1 = vtkUnstructuredGrid::SafeDownCast(output1);
  vtkUnstructuredGrid *outputCT2 = vtkUnstructuredGrid::SafeDownCast(output2);
  vtkUnstructuredGrid *outputCT3 = vtkUnstructuredGrid::SafeDownCast(output3);

  // Wrap
  bottleneckDistance_.setWrapper(this); 

  CTPersistenceDiagram1_ = vtkUnstructuredGrid::SafeDownCast(input1);
  CTPersistenceDiagram2_ = vtkUnstructuredGrid::SafeDownCast(input2);

  if (!CTPersistenceDiagram1_ || !CTPersistenceDiagram2_ ||
    !outputCT3) return -1;

  int dataType1 = CTPersistenceDiagram1_->GetCellData()->GetArray("Persistence")->GetDataType();
  int dataType2 = CTPersistenceDiagram2_->GetCellData()->GetArray("Persistence")->GetDataType();
  if (dataType1 != dataType2) return -1;

  // Call package
  int status = 0;
  switch (dataType1) {
    vtkTemplateMacro(({

      vector<tuple<idVertex, NodeType, idVertex, NodeType, VTK_TT, idVertex,
        VTK_TT, float, float, float, VTK_TT, float, float, float> >*
        CTDiagram1 =
        new vector<tuple<idVertex, NodeType, idVertex, NodeType, VTK_TT, idVertex,
        VTK_TT, float, float, float, VTK_TT, float, float, float> >();

      vector<tuple<idVertex, NodeType, idVertex, NodeType, VTK_TT, idVertex,
        VTK_TT, float, float, float, VTK_TT, float, float, float> >*
        CTDiagram2 =
        new vector<tuple<idVertex, NodeType, idVertex, NodeType, VTK_TT, idVertex,
        VTK_TT, float, float, float, VTK_TT, float, float, float> >();

      status = getPersistenceDiagram<VTK_TT>(
        CTDiagram1,
        CTPersistenceDiagram1_, Spacing, 0);

      status = getPersistenceDiagram<VTK_TT>(
        CTDiagram2,
        CTPersistenceDiagram2_, Spacing, 1);

      bottleneckDistance_.setCTDiagram1(CTDiagram1);
      bottleneckDistance_.setCTDiagram2(CTDiagram2);

      string wassersteinMetric = WassersteinMetric;
      bottleneckDistance_.setWasserstein(wassersteinMetric);

      // Empty matchings.
      vector<tuple<idVertex, idVertex, VTK_TT> >* matchings
        = new vector<tuple<idVertex, idVertex, VTK_TT> >();
      bottleneckDistance_.setOutputMatchings(matchings);

      // Exec.
      bool usePersistenceMetric = UsePersistenceMetric;
      double alpha = Alpha;
      status = bottleneckDistance_.execute<VTK_TT>(usePersistenceMetric, alpha);

      // Apply results to outputs 0 and 1.
      status = augmentPersistenceDiagrams<VTK_TT>(
        CTDiagram1,
        CTDiagram2,
        matchings,
        CTPersistenceDiagram1_,
        CTPersistenceDiagram2_);

      bool useOutputMatching = UseOutputMatching;
      bool useGeometricSpacing = UseGeometricSpacing;

      // Apply results to output 2.
      if (useOutputMatching) {
        status = getMatchingMesh<VTK_TT>(
          CTDiagram1, CTDiagram2, matchings,
          useGeometricSpacing, Spacing);
      }

      if (status != 0) { return status; }
  
    }));
  }

  // Set output.
  outputCT1->ShallowCopy(CTPersistenceDiagram1_);
  outputCT2->ShallowCopy(CTPersistenceDiagram2_);
  if (UseOutputMatching)
    outputCT3->ShallowCopy(CTPersistenceDiagram3_);

  return status;
}


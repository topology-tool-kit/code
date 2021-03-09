#include "BottleneckDistance.h"
#include <ttkTrackingFromPersistenceDiagrams.h>

vtkStandardNewMacro(ttkTrackingFromPersistenceDiagrams)

  using dataType = double;

ttkTrackingFromPersistenceDiagrams::ttkTrackingFromPersistenceDiagrams() {

  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(1);
}

ttkTrackingFromPersistenceDiagrams::~ttkTrackingFromPersistenceDiagrams() {
  if(outputMesh_)
    outputMesh_->Delete();
}

int ttkTrackingFromPersistenceDiagrams::FillInputPortInformation(
  int port, vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

int ttkTrackingFromPersistenceDiagrams::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}

int ttkTrackingFromPersistenceDiagrams::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  // Output pointers informations
  vtkInformation *outInfo;

  // Unified bound fields
  outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *mesh = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkUnstructuredGrid::DATA_OBJECT()));

  // Number of input files
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();
  this->printMsg("Number of inputs: " + std::to_string(numInputs));

  // Get input data
  std::vector<vtkDataSet *> input(numInputs);
  for(int i = 0; i < numInputs; i++) {
    input[i] = vtkDataSet::GetData(inputVector[0], i);
  }

  std::vector<std::vector<diagramTuple>> inputPersistenceDiagrams(
    (unsigned long)numInputs, std::vector<diagramTuple>());

  std::vector<vtkSmartPointer<vtkUnstructuredGrid>> outputPersistenceDiagrams(
    (unsigned long)2 * numInputs - 2,
    vtkSmartPointer<vtkUnstructuredGrid>::New());

  std::vector<std::vector<matchingTuple>> outputMatchings(
    (unsigned long)numInputs - 1, std::vector<matchingTuple>());

  // Input parameters.
  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
  double alpha = Alpha;
  double tolerance = Tolerance;
  bool is3D = Is3D;
  std::string wasserstein = WassersteinMetric;

  // Transform inputs into the right structure.
  for(int i = 0; i < numInputs; ++i) {
    vtkUnstructuredGrid *grid1 = vtkUnstructuredGrid::New();
    grid1->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i]));
    this->getPersistenceDiagram(inputPersistenceDiagrams[i], grid1, spacing, 0);
  }

  this->performMatchings(
    numInputs, inputPersistenceDiagrams, outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein, tolerance, is3D,
    alpha, // Blending
    PX, PY, PZ, PS, PE // Coefficients
  );

  // Get back meshes.
  //  #pragma omp parallel for num_threads(ThreadNumber)
  for(int i = 0; i < numInputs - 1; ++i) {
    vtkUnstructuredGrid *grid1 = vtkUnstructuredGrid::New();
    grid1->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i]));
    vtkUnstructuredGrid *grid2 = vtkUnstructuredGrid::New();
    grid2->ShallowCopy(vtkUnstructuredGrid::SafeDownCast(input[i + 1]));

    vtkUnstructuredGrid *CTPersistenceDiagram1_
      = vtkUnstructuredGrid::SafeDownCast(grid1);
    vtkUnstructuredGrid *CTPersistenceDiagram2_
      = vtkUnstructuredGrid::SafeDownCast(grid2);

    int status = this->augmentPersistenceDiagrams(
      inputPersistenceDiagrams[i], inputPersistenceDiagrams[i + 1],
      outputMatchings[i], CTPersistenceDiagram1_, CTPersistenceDiagram2_);
    if(status < 0)
      return -2;

    outputPersistenceDiagrams[2 * i]->ShallowCopy(CTPersistenceDiagram1_);
    outputPersistenceDiagrams[2 * i + 1]->ShallowCopy(CTPersistenceDiagram2_);
  }

  auto numPersistenceDiagramsInput
    = (int)outputPersistenceDiagrams.size(); // numInputs;

  for(int i = 0; i < numPersistenceDiagramsInput; ++i) {
    vtkSmartPointer<vtkUnstructuredGrid> grid = outputPersistenceDiagrams[i];
    if(!grid || !grid->GetCellData()
       || !grid->GetCellData()->GetArray("Persistence")) {
      this->printErr("Inputs are not persistence diagrams");
      return 0;
    }

    // Check if inputs have the same data type and the same number of points
    if(grid->GetCellData()->GetArray("Persistence")->GetDataType()
       != outputPersistenceDiagrams[0]
            ->GetCellData()
            ->GetArray("Persistence")
            ->GetDataType()) {
      this->printErr("Inputs of different data types");
      return 0;
    }
  }

  for(int i = 0; i < numPersistenceDiagramsInput - 1; ++i) {
    vtkSmartPointer<vtkUnstructuredGrid> grid1 = outputPersistenceDiagrams[i];
    vtkSmartPointer<vtkUnstructuredGrid> grid2
      = outputPersistenceDiagrams[i + 1];
    if(i % 2 == 1 && i < numPersistenceDiagramsInput - 1
       && grid1->GetCellData()->GetNumberOfTuples()
            != grid2->GetCellData()->GetNumberOfTuples()) {
      this->printErr("Inconsistent length or order of input diagrams.");
      return 0;
    }
  }

  outputMesh_ = vtkUnstructuredGrid::New();
  vtkUnstructuredGrid *outputMesh = vtkUnstructuredGrid::SafeDownCast(mesh);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkDoubleArray> costScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> valueScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> matchingIdScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> lengthScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> timeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> componentIds
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> pointTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  costScalars->SetName("Cost");
  persistenceScalars->SetName("Persistence");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("CriticalType");

  // (+ vertex id)
  std::vector<trackingTuple>
    trackingsBase; // structure containing all trajectories
  this->performTracking(
    inputPersistenceDiagrams, outputMatchings, trackingsBase);

  std::vector<std::set<int>> trackingTupleToMerged(
    trackingsBase.size(), std::set<int>());
  if(DoPostProc)
    this->performPostProcess(inputPersistenceDiagrams, trackingsBase,
                                       trackingTupleToMerged, PostProcThresh);

  // bool Is3D = true;
  bool useGeometricSpacing = UseGeometricSpacing;
  // auto spacing = (float) Spacing;

  // std::vector<trackingTuple> trackings = *trackingsBase;
  // Row = iteration number
  // Col = (id in pd 2, arrival point id)

  // Build mesh.
  buildMesh(
    trackingsBase, outputMatchings, inputPersistenceDiagrams,
    useGeometricSpacing, spacing, DoPostProc, trackingTupleToMerged, points,
    persistenceDiagram, costScalars, persistenceScalars, valueScalars, matchingIdScalars,
    lengthScalars, timeScalars, componentIds, pointTypeScalars);

  outputMesh_->ShallowCopy(persistenceDiagram);
  outputMesh->ShallowCopy(outputMesh_);

  return 1;
}

int ttkTrackingFromPersistenceDiagrams::buildMesh(
  std::vector<trackingTuple> &trackings,
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
  vtkIntArray *pointTypeScalars) {
  int currentVertex = 0;
  for(unsigned int k = 0; k < trackings.size(); ++k) {
    trackingTuple tt = trackings.at((unsigned long)k);

    int numStart = std::get<0>(tt);
    //     int numEnd = std::get<1>(tt);
    std::vector<BIdVertex> chain = std::get<2>(tt);
    int chainLength = chain.size();

    if(chain.size() <= 1) {
      // printErr("Got an unexpected 0-size chain.");
      return 0;
    }

    for(int c = 0; c < (int)chain.size() - 1; ++c) {
      std::vector<matchingTuple> &matchings1 = outputMatchings[numStart + c];
      int d1id = numStart + c;
      int d2id = d1id + 1; // c % 2 == 0 ? d1id + 1 : d1id;
      std::vector<diagramTuple> &diagram1 = inputPersistenceDiagrams[d1id];
      std::vector<diagramTuple> &diagram2 = inputPersistenceDiagrams[d2id];

      // Insert segments
      vtkIdType ids[2];
      auto n1 = (int)chain.at((unsigned long)c);
      auto n2 = (int)chain.at((unsigned long)c + 1);

      // Search for right matching.
      double cost = 0.0;
      for(matchingTuple tup : matchings1) {
        auto d1id1 = (int)std::get<0>(tup);
        if(d1id1 == n1) {
          cost = std::get<2>(tup);
          break;
        }
      }

      diagramTuple tuple1 = diagram1[n1];
      diagramTuple tuple2 = diagram2[n2];

      ttk::DecoratedDiagramTuple::DataPoint startPoint, endPoint;
      auto barycenter = [](ttk::DecoratedDiagramTuple::DataPoint a,
                           ttk::DecoratedDiagramTuple::DataPoint b) {
        return ttk::DecoratedDiagramTuple::DataPoint((a.val + b.val) / 2,
          (a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2);
      };

      auto extractRepresentingPoint = [&](diagramTuple tup) {
        ttk::DecoratedDiagramTuple::DataPoint ret;
        // I guess this should be simplified by assuming birthType < deathType
        bool t1Min = tup.birthType == BLocalMin;
        bool t1Max = tup.birthType == BLocalMax;
        bool t2Min = tup.deathType == BLocalMin;
        bool t2Max = tup.deathType == BLocalMax;
        bool bothEx1 = (t1Min && t2Max) || (t1Max && t2Min);
        if(bothEx1)
          ret = t2Max ? tup.deathPoint : tup.birthPoint;
        else
          ret = t2Max ? tup.deathPoint
              : t1Min ? tup.birthPoint
                      : barycenter(tup.birthPoint, tup.deathPoint);
        if(useGeometricSpacing)
          ret.z += spacing * (numStart + c);

        return ret;
      };
      auto typeOfTuple = [](diagramTuple tup) {
        return tup.birthType == BLocalMax || tup.deathType == BLocalMax ? BLocalMax
             : tup.birthType == BLocalMin || tup.deathType == BLocalMin ? BLocalMin
             : tup.birthType == BSaddle2  || tup.deathType == BSaddle2  ? BSaddle2
                                                                        : BSaddle1;
      };
      startPoint = extractRepresentingPoint(tuple1);

      points->InsertNextPoint(startPoint.x, startPoint.y, startPoint.z);
      ids[0] = 2 * currentVertex;
      pointTypeScalars->InsertTuple1(ids[0], (double)(int)typeOfTuple(tuple1));
      timeScalars->InsertTuple1(ids[0], (double)numStart + c);
      persistenceScalars->InsertTuple1(ids[0], tuple1.persistence);
      componentIds->InsertTuple1(ids[0], k);

      endPoint = extractRepresentingPoint(tuple2);
      if(useGeometricSpacing)
        endPoint.z += spacing;

      points->InsertNextPoint(endPoint.x, endPoint.y, endPoint.z);
      ids[1] = 2 * currentVertex + 1;
      pointTypeScalars->InsertTuple1(ids[1], (double)(int)typeOfTuple(tuple2));
      timeScalars->InsertTuple1(ids[1], (double)numStart + c);
      persistenceScalars->InsertTuple1(ids[1], tuple2.persistence);
      componentIds->InsertTuple1(ids[1], k);

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);

      costScalars->InsertTuple1(currentVertex, cost);
      valueScalars->InsertTuple1(
        currentVertex, barycenter(startPoint, endPoint).val);
      matchingIdScalars->InsertTuple1(currentVertex, currentVertex);
      lengthScalars->InsertTuple1(currentVertex, chainLength);

      currentVertex++;
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(costScalars);
  persistenceDiagram->GetPointData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(valueScalars);
  persistenceDiagram->GetCellData()->AddArray(matchingIdScalars);
  persistenceDiagram->GetCellData()->AddArray(lengthScalars);
  persistenceDiagram->GetPointData()->AddArray(timeScalars);
  persistenceDiagram->GetPointData()->AddArray(componentIds);
  persistenceDiagram->GetPointData()->AddArray(pointTypeScalars);

  return 0;
}

// Warn: this is a duplicate from ttkBottleneckDistance
int ttkTrackingFromPersistenceDiagrams::getPersistenceDiagram(
  std::vector<diagramTuple> &diagram,
  const vtkSmartPointer<vtkUnstructuredGrid> &CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber) {
  ttkSimplexIdTypeArray *vertexIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(
      CTPersistenceDiagram_->GetPointData()->GetArray(
        ttk::VertexScalarFieldName));

  vtkIntArray *nodeTypeScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("CriticalType"));

  ttkSimplexIdTypeArray *pairIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(
      CTPersistenceDiagram_->GetCellData()->GetArray("PairIdentifier"));

  vtkIntArray *extremumIndexScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("PairType"));

  vtkDoubleArray *persistenceScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("Persistence"));

  vtkDoubleArray *birthScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Birth"));

  vtkDoubleArray *deathScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Death"));

  vtkPoints *points = (CTPersistenceDiagram_->GetPoints());
  if(!pairIdentifierScalars)
    return -2;

  auto pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();

  // Continuous indexing (no gap in indices)
  for(int pairIndex = 0; pairIndex < pairingsSize; ++pairIndex) {
    const float indexOfPair = pairIndex;
    if(*pairIdentifierScalars->GetTuple(pairIndex) != -1) // except diagonal
      pairIdentifierScalars->SetTuple(pairIndex, &indexOfPair);
  }

  auto s = (float)0.0;

  if(!deathScalars != !birthScalars)
    return -2;
  bool is2D = !deathScalars && !birthScalars;
  bool is3D = !is2D;
  if(Is3D && !is3D)
    Is3D = false;
  if(!Is3D && diagramNumber == 1)
    s = (float)spacing;

  if(pairingsSize < 1 || !vertexIdentifierScalars || !nodeTypeScalars
     || !persistenceScalars || !extremumIndexScalars || !points)
    return -2;

  diagram.resize((unsigned long)pairingsSize);
  int nbNonCompact = 0;

  for(int i = 0; i < pairingsSize; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    int index1 = 2 * i;
    double *coords1 = points->GetPoint(index1);
    auto x1 = (float)coords1[0];
    auto y1 = (float)coords1[1];
    auto z1 = (float)coords1[2];

    int index2 = index1 + 1;
    double *coords2 = points->GetPoint(index2);
    auto x2 = (float)coords2[0];
    auto y2 = (float)coords2[1];
    auto z2 = (float)coords2[2];

    dataType value1 = (!birthScalars) ? (dataType)x1
                                      : (dataType)birthScalars->GetValue(2 * i);
    dataType value2 = (!deathScalars)
                        ? (dataType)y2
                        : (dataType)deathScalars->GetValue(2 * i + 1);

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize)
      diagram.at(pairIdentifier)
        = std::make_tuple(vertexId1, (BNodeType)nodeType1, vertexId2,
                          (BNodeType)nodeType2, (dataType)persistence, pairType,
                          value1, x1, y1, z1 + s, value2, x2, y2, z2 + s);

    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        this->printWrn("Diagram pair identifiers must be compact (not exceed "
                       "the diagram size).");
      }
    }
  }

  if(nbNonCompact > 0) {
    this->printWrn("Missed " + std::to_string(nbNonCompact)
                   + " pairs due to non-compactness.");
  }

  sort(diagram.begin(), diagram.end(),
       [](const diagramTuple &a, const diagramTuple &b) -> bool {
         return a.birthPoint.val < b.birthPoint.val;
       });

  return 0;
}

// Warn: this is a duplicate from ttkBottleneckDistance
int ttkTrackingFromPersistenceDiagrams::augmentPersistenceDiagrams(
  const std::vector<diagramTuple> &diagram1,
  const std::vector<diagramTuple> &diagram2,
  const std::vector<matchingTuple> &matchings,
  vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram1_,
  vtkSmartPointer<vtkUnstructuredGrid> CTPersistenceDiagram2_) {

  auto diagramSize1 = (BIdVertex)diagram1.size();
  auto diagramSize2 = (BIdVertex)diagram2.size();
  auto matchingsSize = (BIdVertex)matchings.size();

  vtkSmartPointer<vtkIntArray> matchingIdentifiers1
    = vtkSmartPointer<vtkIntArray>::New();
  matchingIdentifiers1->SetName("MatchingIdentifier");

  vtkSmartPointer<vtkIntArray> matchingIdentifiers2
    = vtkSmartPointer<vtkIntArray>::New();
  matchingIdentifiers2->SetName("MatchingIdentifier");

  if(matchingsSize > 0) {
    ttk::SimplexId ids[2];
    matchingIdentifiers1->SetNumberOfComponents(1);
    matchingIdentifiers2->SetNumberOfComponents(1);
    matchingIdentifiers1->SetNumberOfTuples(diagramSize1);
    matchingIdentifiers2->SetNumberOfTuples(diagramSize2);

    // Unaffected by default
    for(BIdVertex i = 0; i < diagramSize1; ++i)
      matchingIdentifiers1->InsertTuple1(i, -1);
    for(BIdVertex i = 0; i < diagramSize2; ++i)
      matchingIdentifiers2->InsertTuple1(i, -1);

    // Last cell = junction
    if(diagramSize1
       < CTPersistenceDiagram1_->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers1->InsertTuple1(diagramSize1, -1);
      matchingIdentifiers1->InsertTuple1(diagramSize1 + 1, -1);
    }
    if(diagramSize2
       < CTPersistenceDiagram2_->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers2->InsertTuple1(diagramSize2, -1);
      matchingIdentifiers2->InsertTuple1(diagramSize2 + 1, -1);
    }

    // Affect bottleneck matchings
    int pairingIndex = 0;
    for(BIdVertex i = 0; i < matchingsSize; ++i) {
      matchingTuple t = matchings.at((unsigned long)i);
      ids[0] = std::get<0>(t);
      ids[1] = std::get<1>(t);
      matchingIdentifiers1->InsertTuple1(ids[0], pairingIndex);
      matchingIdentifiers2->InsertTuple1(ids[1], pairingIndex);
      pairingIndex++;
    }

    CTPersistenceDiagram1_->GetCellData()->AddArray(matchingIdentifiers1);
    CTPersistenceDiagram2_->GetCellData()->AddArray(matchingIdentifiers2);
  }

  return 0;
}

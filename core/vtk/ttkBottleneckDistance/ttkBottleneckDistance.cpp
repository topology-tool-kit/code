#include "ttkBottleneckDistance.h"
#include <ttkMacros.h>
#include <ttkUtils.h>

// VTK includes
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

// Misc.
#include <cstdlib>
#include <ctime>
#include <random>

vtkStandardNewMacro(ttkBottleneckDistance);

ttkBottleneckDistance::ttkBottleneckDistance() {
  SetNumberOfInputPorts(2);
  SetNumberOfOutputPorts(3);
}

int ttkBottleneckDistance::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(ttkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkBottleneckDistance::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0 || port == 1 || port == 2) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkBottleneckDistance::generatePersistenceDiagram(
  std::vector<diagramTuple> &diagram, const int size) {
  // srand(time(NULL));
  int vertexId1 = 1;
  int vertexId2 = 2;

  // diagram.resize(size);
  // diagram.push_back(std::make_tuple(
  //    0, BLocalMin,
  //    1, BLocalMax,
  //    1.0,
  //    2,
  //    0.0, 0.0001, 0.0001, 0.0,
  //    1.0, 0.0001, 1.0, 0.0
  //));

  std::default_random_engine generator1(1998985);
  std::default_random_engine generator2(8584584);
  std::uniform_real_distribution<> dis1(0.0, 1.0);
  std::uniform_real_distribution<> dis2(0.0, 1.0);

  for(int i = 1; i < size; ++i) {
    int r0 = 2; // (rand() % 3);
    float r1
      = 0.001f + static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    // (float) dis1(generator1); // static_cast <float> (rand()) / static_cast
    // <float> (RAND_MAX);
    float r2
      = 0.001f + static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    // (float) dis2(generator2); // static_cast <float> (rand()) / static_cast
    // <float> (RAND_MAX); r2 *= 0.1f;

    // BLocalMin BSaddle1 BSaddle2 BLocalMax
    int pairType = r0; // (0/min, 1/saddle, 2/max)
    BNodeType nodeType1; //
    BNodeType nodeType2; //
    switch(pairType) {
      // case 0:
      //   nodeType1 = BLocalMin;
      //   nodeType2 = BSaddle1;
      //   break;
      // case 1:
      //   nodeType1 = BSaddle1;
      //   nodeType2 = BSaddle2;
      //   break;
      case 2:
      default:
        nodeType1 = BSaddle2;
        nodeType2 = BLocalMax;
        break;
        // default:
        //   nodeType1 = (BNodeType) -1;
        //   nodeType2 = (BNodeType) -1;
    }

    float x1 = 0.5f * r1;
    float y1 = x1;
    float z1 = 0.f;

    float x2 = x1;
    float y2 = x1 + 0.5f * r2; // x1 + rand(0.5)
    float z2 = 0.f; // 0

    auto value1 = (dataType)x1;
    auto value2 = (dataType)y2;

    dataType persistence = y2 - x1;

    diagram.push_back(std::make_tuple(vertexId1, nodeType1, vertexId2,
                                      nodeType2, persistence, pairType, value1,
                                      x1, y1, z1, value2, x2, y2, z2));

    vertexId1++;
    vertexId2++;
  }

  sort(diagram.begin(), diagram.end(),
       [](const diagramTuple &a, const diagramTuple &b) -> bool {
         return a.birthPoint.val < b.birthPoint.val;
       });

  return 1;
}

// Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
int ttkBottleneckDistance::getPersistenceDiagram(
  std::vector<diagramTuple> &diagram,
  vtkUnstructuredGrid *const CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber) {

  auto pointData = CTPersistenceDiagram_->GetPointData();
  auto cellData = CTPersistenceDiagram_->GetCellData();

  if(pointData == nullptr || cellData == nullptr) {
    return -1;
  }

  auto vertexIdentifierScalars = ttkSimplexIdTypeArray::SafeDownCast(
    pointData->GetArray(ttk::VertexScalarFieldName));

  auto nodeTypeScalars
    = vtkIntArray::SafeDownCast(pointData->GetArray("CriticalType"));

  auto pairIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray("PairIdentifier"));

  auto extremumIndexScalars
    = vtkIntArray::SafeDownCast(cellData->GetArray("PairType"));

  auto persistenceScalars
    = vtkDoubleArray::SafeDownCast(cellData->GetArray("Persistence"));

  auto birthScalars
    = vtkDoubleArray::SafeDownCast(pointData->GetArray("Birth"));

  auto deathScalars
    = vtkDoubleArray::SafeDownCast(pointData->GetArray("Death"));

  vtkPoints *points = CTPersistenceDiagram_->GetPoints();
  if(!pairIdentifierScalars)
    return -2;

  auto pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();

  // Continuous indexing (no gap in indices)
  for(int pairIndex = 0; pairIndex < pairingsSize; ++pairIndex) {
    const float indexOfPair = pairIndex;
    if(*pairIdentifierScalars->GetTuple(pairIndex) != -1) // except diagonal
      pairIdentifierScalars->SetTuple(pairIndex, &indexOfPair);
  }

  float s{0.0};

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
        std::stringstream msg;
        msg << "Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        this->printWrn(msg.str());
      }
    }
  }

  if(nbNonCompact > 0) {
    std::stringstream msg;
    msg << "Missed " << nbNonCompact << " pairs due to non-compactness."
        << std::endl;
    this->printWrn(msg.str());
  }

  sort(diagram.begin(), diagram.end(),
       [](const diagramTuple &a, const diagramTuple &b) -> bool {
         return a.birthPoint.val < b.birthPoint.val;
       });

  return 1;
}

// Warn: this is duplicated in ttkTrackingFromPersistenceDiagrams
int ttkBottleneckDistance::augmentPersistenceDiagrams(
  const std::vector<diagramTuple> &diagram1,
  const std::vector<diagramTuple> &diagram2,
  const std::vector<matchingTuple> &matchings,
  vtkUnstructuredGrid *const CTPersistenceDiagram1,
  vtkUnstructuredGrid *const CTPersistenceDiagram2) {

  auto diagramSize1 = (BIdVertex)diagram1.size();
  auto diagramSize2 = (BIdVertex)diagram2.size();
  auto matchingsSize = (BIdVertex)matchings.size();

  vtkNew<vtkIntArray> matchingIdentifiers1{};
  matchingIdentifiers1->SetName("MatchingIdentifier");

  vtkNew<vtkIntArray> matchingIdentifiers2{};
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
       < CTPersistenceDiagram1->GetCellData()->GetNumberOfTuples()) {
      matchingIdentifiers1->InsertTuple1(diagramSize1, -1);
      matchingIdentifiers1->InsertTuple1(diagramSize1 + 1, -1);
    }
    if(diagramSize2
       < CTPersistenceDiagram2->GetCellData()->GetNumberOfTuples()) {
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

    CTPersistenceDiagram1->GetCellData()->AddArray(matchingIdentifiers1);
    CTPersistenceDiagram2->GetCellData()->AddArray(matchingIdentifiers2);
  }

  return 1;
}

int ttkBottleneckDistance::translateSecondDiagram(
  vtkUnstructuredGrid *outputCT2, double &spacing) {

  vtkNew<vtkPoints> points2{};
  vtkPoints *points = (outputCT2->GetPoints());
  auto pairIdentifierScalars = ttkSimplexIdTypeArray::SafeDownCast(
    outputCT2->GetCellData()->GetArray("PairIdentifier"));

  if(pairIdentifierScalars == nullptr) {
    return -1;
  }

  auto pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();

  for(int i = 0; i < pairingsSize; ++i) {
    int index1 = 2 * i;
    double *coords1 = points->GetPoint(index1);
    auto x1 = coords1[0];
    auto y1 = coords1[1];
    auto z1 = coords1[2] + spacing;

    int index2 = index1 + 1;
    double *coords2 = points->GetPoint(index2);
    auto x2 = coords2[0];
    auto y2 = coords2[1];
    auto z2 = coords2[2] + spacing;
    points2->InsertNextPoint(x1, y1, z1);
    points2->InsertNextPoint(x2, y2, z2);
  }

  outputCT2->SetPoints(points2);

  return 1;
}

int ttkBottleneckDistance::getMatchingMesh(
  vtkUnstructuredGrid *const outputCT3,
  const std::vector<diagramTuple> &diagram1,
  const std::vector<diagramTuple> &diagram2,
  const std::vector<matchingTuple> &matchings,
  const bool useGeometricSpacing,
  const double spacing,
  const bool is2D) {

  vtkNew<vtkPoints> points{};
  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetName("Cost");

  vtkNew<vtkIntArray> matchingIdScalars{};
  matchingIdScalars->SetName("MatchingIdentifier");

  auto matchingsSize = (BIdVertex)matchings.size();

  // Build matchings.
  if(matchingsSize > 0) {

    for(BIdVertex i = 0; i < matchingsSize; ++i) {
      vtkIdType ids[2];

      matchingTuple t = matchings.at((unsigned long)i);
      auto n1 = (int)std::get<0>(t);
      auto n2 = (int)std::get<1>(t);

      diagramTuple tuple1 = diagram1.at(n1);
      diagramTuple tuple2 = diagram2.at(n2);

      double x1, y1, z1, x2, y2, z2;

      bool linkMiddles = false;
      if(linkMiddles) {
        x1 = (tuple1.birthPoint.x + tuple1.deathPoint.x) / 2;
        y1 = (tuple1.birthPoint.y + tuple1.deathPoint.y) / 2;
        z1 = (tuple1.birthPoint.z + tuple1.deathPoint.z) / 2;
      } else {
        BNodeType t11 = tuple1.birthType;
        BNodeType t12 = tuple1.deathType;
        bool t11Max = t11 == BLocalMin || t11 == BLocalMax;
        bool t12Max = t12 == BLocalMin || t12 == BLocalMax;
        if(is2D) { // Quickchage for highlighting 2D matching
          if(t11 != BLocalMax && t12 != BLocalMax) {
            t11Max = t11 != BLocalMin;
            t12Max = t12 != BLocalMin;
          }
        }
        x1 = t12Max   ? tuple1.deathPoint.x
             : t11Max ? tuple1.birthPoint.x
                      : (tuple1.birthPoint.x + tuple1.deathPoint.x) / 2;
        y1 = t12Max   ? tuple1.deathPoint.y
             : t11Max ? tuple1.birthPoint.y
                      : (tuple1.birthPoint.y + tuple1.deathPoint.y) / 2;
        z1 = t12Max   ? tuple1.deathPoint.z
             : t11Max ? tuple1.birthPoint.z
                      : (tuple1.birthPoint.z + tuple1.deathPoint.z) / 2;
      }
      points->InsertNextPoint(x1, y1, z1);

      if(linkMiddles) {
        x2 = (tuple2.birthPoint.x + tuple2.deathPoint.x) / 2;
        y2 = (tuple2.birthPoint.y + tuple2.deathPoint.y) / 2;
        z2 = (tuple2.birthPoint.z + tuple2.deathPoint.z) / 2;
        if(useGeometricSpacing)
          z2 += spacing;
      } else {
        BNodeType t21 = tuple2.birthType;
        BNodeType t22 = tuple2.deathType;
        bool t21Max = t21 == BLocalMin || t21 == BLocalMax;
        bool t22Max = t22 == BLocalMin || t22 == BLocalMax;
        if(is2D) { // Quickchage for highlighting 2D matching
          if(t21 != BLocalMax && t22 != BLocalMax) {
            t21Max = t21 != BLocalMin;
            t22Max = t22 != BLocalMin;
          }
        }
        x2 = t22Max   ? tuple2.deathPoint.x
             : t21Max ? tuple2.birthPoint.x
                      : (tuple2.birthPoint.x + tuple2.deathPoint.x) / 2;
        y2 = t22Max   ? tuple2.deathPoint.y
             : t21Max ? tuple2.birthPoint.y
                      : (tuple2.birthPoint.y + tuple2.deathPoint.y) / 2;
        z2 = t22Max   ? tuple2.deathPoint.z
             : t21Max ? tuple2.birthPoint.z
                      : (tuple2.birthPoint.z + tuple2.deathPoint.z) / 2;
      }
      points->InsertNextPoint(x2, y2, z2);

      ids[0] = 2 * i;
      ids[1] = 2 * i + 1;

      persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);

      persistenceScalars->InsertTuple1(i, std::get<2>(t));
      matchingIdScalars->InsertTuple1(i, i);
    }
  }

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(matchingIdScalars);

  outputCT3->ShallowCopy(persistenceDiagram);

  return 1;
}

int ttkBottleneckDistance::doBenchmark() {
  std::vector<diagramTuple> CTDiagram1;
  std::vector<diagramTuple> CTDiagram2;

  int benchmarkSize = BenchmarkSize;
  int status = 0;
  status = generatePersistenceDiagram(CTDiagram1, benchmarkSize);
  if(status < 0)
    return status;
  status = generatePersistenceDiagram(CTDiagram2, 4 * benchmarkSize);
  if(status < 0)
    return status;

  this->setPersistencePercentThreshold(Tolerance);
  this->setPX(PX);
  this->setPY(PY);
  this->setPZ(PZ);
  this->setPE(PE);
  this->setPS(PS);
  this->setCTDiagram1(&CTDiagram1);
  this->setCTDiagram2(&CTDiagram2);

  std::string wassersteinMetric = WassersteinMetric;
  this->setWasserstein(wassersteinMetric);
  std::string algorithm = DistanceAlgorithm;
  this->setAlgorithm(algorithm);
  int pvAlgorithm = PVAlgorithm;
  this->setPVAlgorithm(pvAlgorithm);
  // this->setThreadNumber(thread);

  // Empty matchings.
  auto matchings = new std::vector<diagramTuple>();
  this->setOutputMatchings(matchings);

  // Exec.
  bool usePersistenceMetric = UsePersistenceMetric;
  status = this->execute(usePersistenceMetric);

  return status;
}

int ttkBottleneckDistance::RequestData(vtkInformation * /*request*/,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {
  int benchmarkSize = BenchmarkSize;
  bool benchmark = benchmarkSize > 0;
  if(benchmark) {
    return doBenchmark();
  }

  auto outputCT1 = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto outputCT2 = vtkUnstructuredGrid::GetData(outputVector, 1);
  auto outputCT3 = vtkUnstructuredGrid::GetData(outputVector, 2);

  // Wrap
  // this->setWrapper(this);
  this->setPersistencePercentThreshold(Tolerance);
  this->setPX(PX);
  this->setPY(PY);
  this->setPZ(PZ);
  this->setPE(PE);
  this->setPS(PS);

  auto CTPersistenceDiagram1 = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto CTPersistenceDiagram2 = vtkUnstructuredGrid::GetData(inputVector[1]);

  if(!CTPersistenceDiagram1 || !CTPersistenceDiagram2 || !outputCT3) {
    this->printErr("Input grids should be non-NULL");
    return 0;
  }

  int dataType1 = CTPersistenceDiagram1->GetCellData()
                    ->GetArray("Persistence")
                    ->GetDataType();
  int dataType2 = CTPersistenceDiagram2->GetCellData()
                    ->GetArray("Persistence")
                    ->GetDataType();
  if(dataType1 != dataType2) {
    this->printErr("Persistence array data type should be the same");
    return 0;
  }

  auto birthScalars1 = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram1->GetPointData()->GetArray("Birth"));
  auto deathScalars1 = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram1->GetPointData()->GetArray("Death"));
  auto birthScalars2 = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram2->GetPointData()->GetArray("Birth"));
  auto deathScalars2 = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram2->GetPointData()->GetArray("Death"));
  bool is2D1 = !deathScalars1 && !birthScalars1;
  bool is2D2 = !deathScalars2 && !birthScalars2;
  if(is2D1 != is2D2) {
    this->printErr("Diagrams should not be embedded");
    return 0;
  }
  bool is2D = is2D1;

  // Call package
  int status = 0;

  //  switch (dataType1) {
  //    vtkTemplateMacro(({
  // TODO template my methods
  std::vector<diagramTuple> CTDiagram1;
  std::vector<diagramTuple> CTDiagram2;

  status = getPersistenceDiagram(
    CTDiagram1, CTPersistenceDiagram1, Spacing, 0);
  if(status < 0) {
    this->printErr("Could not extract diagram from first input data-set");
    return 0;
  }

  status = getPersistenceDiagram(
    CTDiagram2, CTPersistenceDiagram2, Spacing, 1);
  if(status < 0) {
    this->printErr("Could not extract diagram from second input data-set");
    return 0;
  }

  this->setCTDiagram1(&CTDiagram1);
  this->setCTDiagram2(&CTDiagram2);

  this->setWasserstein(WassersteinMetric);
  this->setAlgorithm(DistanceAlgorithm);
  this->setPVAlgorithm(PVAlgorithm);

  // Empty matchings.
  std::vector<matchingTuple> matchings;
  this->setOutputMatchings(&matchings);

  // Exec.
  status = this->execute(UsePersistenceMetric);
  if(status != 0) {
    this->printErr("Base layer failed with error status "
                   + std::to_string(status));
    return 0;
  }

  // Apply results to outputs 0 and 1.
  status = augmentPersistenceDiagrams(
    CTDiagram1, CTDiagram2, matchings, CTPersistenceDiagram1,
    CTPersistenceDiagram2);
  if(status != 1) {
    this->printErr("Could not augment diagrams");
    return 0;
  }

  // Apply results to output 2.
  if(UseOutputMatching) {
    status
      = getMatchingMesh(outputCT3, CTDiagram1, CTDiagram2, matchings,
                                  UseGeometricSpacing, Spacing, is2D);

    if(status != 1) {
      this->printErr("Could not compute matchings");
      return 0;
    }
  }

  // Set output.
  outputCT1->ShallowCopy(CTPersistenceDiagram1);
  outputCT2->DeepCopy(CTPersistenceDiagram2);
  if(UseGeometricSpacing)
    translateSecondDiagram(outputCT2, Spacing);

  return 1;
}

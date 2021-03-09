#include <vtkInformation.h>

#include <ttkMacros.h>
#include <ttkTrackingFromFields.h>
#include <ttkTrackingFromPersistenceDiagrams.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTrackingFromFields);

constexpr unsigned long long str2int(const char *str, int h = 0) {
  return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}

ttkTrackingFromFields::ttkTrackingFromFields() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkTrackingFromFields::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}
int ttkTrackingFromFields::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

// (*) Persistence-driven approach
template <class dataType, class triangulationType>
int ttkTrackingFromFields::trackWithPersistenceMatching(
  vtkDataSet *input,
  vtkUnstructuredGrid *output,
  unsigned long fieldNumber,
  const triangulationType *triangulation) {

  using trackingTuple = ttk::trackingTuple;

  // 1. get persistence diagrams.
  std::vector<std::vector<diagramTuple>> persistenceDiagrams(
    fieldNumber, std::vector<diagramTuple>());

  this->performDiagramComputation<dataType, triangulationType>(
    (int)fieldNumber, persistenceDiagrams, triangulation);

  // 2. call feature tracking with threshold.
  std::vector<std::vector<matchingTuple>> outputMatchings(
    fieldNumber - 1, std::vector<matchingTuple>());

  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
  double alpha = Alpha;
  double tolerance = Tolerance;
  bool is3D = true; // Is3D;
  std::string wasserstein = WassersteinMetric;

  ttk::TrackingFromPersistenceDiagrams tfp{};
  tfp.setThreadNumber(this->threadNumber_);
  tfp.performMatchings(
    (int)fieldNumber, persistenceDiagrams, outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein, tolerance, is3D,
    alpha, // Blending
    PX, PY, PZ, PS, PE // Coefficients
  );

  vtkNew<vtkPoints> points{};
  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};

  vtkNew<vtkDoubleArray> costScalars{};
  vtkNew<vtkDoubleArray> persistenceScalars{};
  vtkNew<vtkDoubleArray> valueScalars{};
  vtkNew<vtkIntArray> matchingIdScalars{};
  vtkNew<vtkIntArray> lengthScalars{};
  vtkNew<vtkIntArray> timeScalars{};
  vtkNew<vtkIntArray> componentIds{};
  vtkNew<vtkIntArray> pointTypeScalars{};

  costScalars->SetName("Cost");
  persistenceScalars->SetName("Persistence");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("CriticalType");

  // (+ vertex id)
  std::vector<trackingTuple> trackingsBase;
  tfp.performTracking(
    persistenceDiagrams, outputMatchings, trackingsBase);

  std::vector<std::set<int>> trackingTupleToMerged(
    trackingsBase.size(), std::set<int>());

  if(DoPostProc) {
    tfp.performPostProcess(persistenceDiagrams, trackingsBase,
                                     trackingTupleToMerged, PostProcThresh);
  }

  bool useGeometricSpacing = UseGeometricSpacing;

  // Build mesh.
  ttkTrackingFromPersistenceDiagrams::buildMesh(
    trackingsBase, outputMatchings, persistenceDiagrams, useGeometricSpacing,
    spacing, DoPostProc, trackingTupleToMerged, points, persistenceDiagram,
    costScalars, persistenceScalars, valueScalars, matchingIdScalars, lengthScalars,
    timeScalars, componentIds, pointTypeScalars);

  output->ShallowCopy(persistenceDiagram);

  return 1;
}

int ttkTrackingFromFields::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);

  // Test validity of datasets
  if(input == nullptr || output == nullptr) {
    return -1;
  }

  // Get number and list of inputs.
  std::vector<vtkDataArray *> inputScalarFieldsRaw;
  std::vector<vtkDataArray *> inputScalarFields;
  const auto pointData = input->GetPointData();
  int numberOfInputFields = pointData->GetNumberOfArrays();
  if(numberOfInputFields < 3) {
    this->printErr("Not enough input fields to perform tracking.");
  }

  vtkDataArray *firstScalarField = pointData->GetArray(0);

  for(int i = 0; i < numberOfInputFields; ++i) {
    vtkDataArray *currentScalarField = pointData->GetArray(i);
    if(currentScalarField == nullptr
       || currentScalarField->GetName() == nullptr) {
      continue;
    }
    std::string sfname{currentScalarField->GetName()};
    if(sfname.rfind("_Order") == (sfname.size() - 6)) {
      continue;
    }
    if(firstScalarField->GetDataType() != currentScalarField->GetDataType()) {
      this->printErr("Inconsistent field data type or size between fields `"
                     + std::string{firstScalarField->GetName()} + "' and `"
                     + sfname + "'");
      return -1;
    }
    inputScalarFieldsRaw.push_back(currentScalarField);
  }

  std::sort(inputScalarFieldsRaw.begin(), inputScalarFieldsRaw.end(),
            [](vtkDataArray *a, vtkDataArray *b) {
              std::string s1 = a->GetName();
              std::string s2 = b->GetName();
              return std::lexicographical_compare(
                s1.begin(), s1.end(), s2.begin(), s2.end());
            });

  numberOfInputFields = inputScalarFieldsRaw.size();
  int end = EndTimestep <= 0 ? numberOfInputFields
                             : std::min(numberOfInputFields, EndTimestep);
  for(int i = StartTimestep; i < end; i += Sampling) {
    vtkDataArray *currentScalarField = inputScalarFieldsRaw[i];
    // Print scalar field names:
    // std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  // Input -> persistence filter.
  std::string algorithm = DistanceAlgorithm;
  int pvalg = PVAlgorithm;
  bool useTTKMethod = false;

  if(pvalg >= 0) {
    switch(pvalg) {
      case 0:
      case 1:
      case 2:
      case 3:
        useTTKMethod = true;
        break;
      case 4:
        break;
      default:
        this->printMsg("Unrecognized tracking method.");
        break;
    }
  } else {
    switch(str2int(algorithm.c_str())) {
      case str2int("0"):
      case str2int("ttk"):
      case str2int("1"):
      case str2int("legacy"):
      case str2int("2"):
      case str2int("geometric"):
      case str2int("3"):
      case str2int("parallel"):
        useTTKMethod = true;
        break;
      case str2int("4"):
      case str2int("greedy"):
        break;
      default:
        this->printMsg("Unrecognized tracking method.");
        break;
    }
  }

  // 0. get data
  int fieldNumber = inputScalarFields.size();
  std::vector<void *> inputFields(fieldNumber);
  for(int i = 0; i < fieldNumber; ++i) {
    inputFields[i] = ttkUtils::GetVoidPointer(inputScalarFields[i]);
  }
  this->setInputScalars(inputFields);

  // 0'. get offsets
  std::vector<ttk::SimplexId *> inputOrders(fieldNumber);
  for(int i = 0; i < fieldNumber; ++i) {
    this->SetInputArrayToProcess(0, 0, 0, 0, inputScalarFields[i]->GetName());
    auto orderArray = this->GetOrderArray(input, 0, 0, false);
    inputOrders[i]
      = static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(orderArray));
  }
  this->setInputOffsets(inputOrders);

  int status = 0;
  if(useTTKMethod) {
    ttkVtkTemplateMacro(
      inputScalarFields[0]->GetDataType(), triangulation->getType(),
      (status = this->trackWithPersistenceMatching<VTK_TT, TTK_TT>(
         input, output, fieldNumber, (TTK_TT *)triangulation->getData())));
  } else {
    this->printMsg("The specified matching method is not supported.");
  }

  return status;
}

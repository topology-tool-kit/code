#include "ttkContourForests.h"

vtkStandardNewMacro(ttkContourForests);

ttkContourForests::ttkContourForests()
  :  // Base //
    FieldId{0},
    InputOffsetFieldId{-1},
    inputOffsetScalarFieldName_{"OutputOffsetScalarField"},
    isLoaded_{},
    lessPartition_{true},
    tree_{},
    // Here the given number of core only serve for preprocess,
    // a clean tree append before the true process and re-set
    // the good number of threads
    contourTree_{new ContourForests()},
    skeletonNodes_{vtkPolyData::New()},
    skeletonArcs_{vtkPolyData::New()},
    segmentation_{},

    // Void //
    voidUnstructuredGrid_{vtkUnstructuredGrid::New()},
    voidPolyData_{vtkPolyData::New()},

    // Configuration //
    varyingMesh_{},
    varyingDataValues_{},
    treeType_{TreeType::Contour},
    showMin_{true},
    showMax_{true},
    showSaddle1_{true},
    showSaddle2_{true},
    showArc_{true},
    arcResolution_{1},
    partitionNum_{-1},
    skeletonSmoothing_{},
    simplificationType_{},
    simplificationThreshold_{},
    simplificationThresholdBuffer_{},

    // Computation handles //
    toUpdateVertexSoSoffsets_{true},
    toComputeContourTree_{true},
    toUpdateTree_{true},
    toComputeSkeleton_{true},
    toComputeSegmentation_{true},

    // Convenient storage //
    deltaScalar_{},
    numberOfVertices_{},
    vertexSoSoffsets_{new vector<int>},
    criticalPoints_{new vector<int>},
    vertexScalars_{nullptr},
    inputScalars_{new vector<vector<double>>},
    inputScalarsName_{new vector<string>},
    samples_{new vector<vector<vector<vector<int>>>>},
    barycenters_{new vector<vector<vector<vector<double>>>>}
{
  contourTree_->setWrapper(this);
  contourTree_->setDebugLevel(debugLevel_);
  UseAllCores = false;

  // VTK Interface //
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);

  triangulation_ = NULL;
  vtkWarningMacro("[ContourForests]: DEPRECATED This plugin will be removed in a futur release, please use FTM instead");
}

ttkContourForests::~ttkContourForests()
{
  // Base //
  delete contourTree_;
  if(skeletonNodes_)
    skeletonNodes_->Delete();
  if(skeletonArcs_)
    skeletonArcs_->Delete();
  // segmentation_->Delete();

  // Void //
  voidUnstructuredGrid_->Delete();
  voidPolyData_->Delete();

  // Convenient storage //
  delete vertexSoSoffsets_;
  delete criticalPoints_;
  delete inputScalars_;
  delete inputScalarsName_;
  delete samples_;
  delete barycenters_;
}

void ttkContourForests::clearSkeleton()
{
  samples_->clear();
  barycenters_->clear();

  skeletonNodes_->Delete();
  skeletonNodes_ = vtkPolyData::New();
  skeletonArcs_->Delete();
  skeletonArcs_ = vtkPolyData::New();
}

void ttkContourForests::clearSegmentation()
{
  segmentation_->Delete();
}

void ttkContourForests::clearTree()
{
  tree_ = nullptr;
  delete contourTree_;
  contourTree_ = new ContourForests();
  contourTree_->setWrapper(this);
  contourTree_->setDebugLevel(debugLevel_);
}

// transmit abort signals -- to copy paste in other wrappers
bool ttkContourForests::needsToAbort()
{
  return GetAbortExecute();
}

// transmit progress status -- to copy paste in other wrappers
int ttkContourForests::updateProgress(const float& progress)
{
  {
    stringstream msg;
    msg << "[ttkContourForests] " << progress * 100 << "% processed...." << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkContourForests::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
  return 1;
}

int ttkContourForests::FillOutputPortInformation(int port, vtkInformation* info)
{
  switch (port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
      break;

    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
      break;

    case 2:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
      break;
  }

  return 1;
}

void ttkContourForests::SetThreads()
{
  if (!UseAllCores)
    threadNumber_ = ThreadNumber;
  else
    threadNumber_ = OsCall::getNumberOfCores();

  if(partitionNum_ != -1){
    toComputeContourTree_ = true;
    toComputeSkeleton_    = true;
    toUpdateTree_         = true;
  }

  Modified();
}

void ttkContourForests::SetThreadNumber(int threadNumber)
{
  ThreadNumber = threadNumber;
  SetThreads();
}

void ttkContourForests::SetDebugLevel(int d)
{
  Debug::setDebugLevel(d);
  debugLevel_ = d;
}

void ttkContourForests::SetUseAllCores(bool onOff)
{
  UseAllCores = onOff;
  SetThreads();
}

void ttkContourForests::SetUseInputOffsetScalarField(bool onOff)
{
  toUpdateVertexSoSoffsets_ = true;
  toComputeContourTree_     = true;
  toUpdateTree_             = true;
  toComputeSkeleton_        = true;
  toComputeSegmentation_    = true;

  useInputOffsetScalarField_ = onOff;
  Modified();
}

void ttkContourForests::SetScalarField(string scalarField)
{
  toComputeContourTree_  = true;
  toComputeSkeleton_     = true;
  toComputeSegmentation_ = true;

  scalarField_ = scalarField;
  Modified();
}

void ttkContourForests::SetTreeType(int treeType)
{
  if (treeType >= 0 && treeType <= 2) {
    toUpdateTree_ = true;
    if (treeType_ == TreeType::Contour ||
        static_cast<TreeType>(treeType) == TreeType::Contour) {
      toComputeContourTree_ = true;
    }

    toComputeSkeleton_     = true;
    toComputeSegmentation_ = true;

    treeType_ = static_cast<TreeType>(treeType);

    Modified();
  }
}

void ttkContourForests::ShowMin(bool state)
{
  toComputeSkeleton_ = true;

  showMin_ = state;
  Modified();
}

void ttkContourForests::ShowMax(bool state)
{
  toComputeSkeleton_ = true;

  showMax_ = state;
  Modified();
}

void ttkContourForests::ShowSaddle1(bool state)
{
  toComputeSkeleton_ = true;

  showSaddle1_ = state;
  Modified();
}

void ttkContourForests::ShowSaddle2(bool state)
{
  toComputeSkeleton_ = true;

  showSaddle2_ = state;
  Modified();
}

void ttkContourForests::ShowArc(bool state)
{
  toComputeSkeleton_ = true;

  showArc_ = state;
  Modified();
}

void ttkContourForests::SetArcResolution(int arcResolution)
{
  if (arcResolution >= 0) {
    toComputeSkeleton_ = true;

    arcResolution_ = arcResolution;
    Modified();
  }
}

void ttkContourForests::SetPartitionNumber(int partitionNum)
{
  partitionNum_ = partitionNum;

  toComputeContourTree_ = true;
  toComputeSkeleton_ = true;
  toUpdateTree_ = true;
  Modified();
}

void ttkContourForests::SetLessPartition(bool l)
{
  lessPartition_ = l;
  Modified();
}

void ttkContourForests::SetSkeletonSmoothing(double skeletonSmoothing)
{
  if (skeletonSmoothing >= 0) {
    toComputeSkeleton_ = true;

    skeletonSmoothing_ = skeletonSmoothing;
    Modified();
  }
}

void ttkContourForests::SetSimplificationType(int type)
{
  simplificationType_ = type;
  Modified();
}

void ttkContourForests::SetSimplificationThreshold(double simplificationThreshold)
{
  if (simplificationThreshold >= 0.0 && simplificationThreshold <= 1.0) {
    simplificationThresholdBuffer_ = simplificationThreshold;
    toComputeContourTree_    = true;
    toUpdateTree_            = true;
    toComputeSkeleton_       = true;
    toComputeSegmentation_   = true;

    Modified();
  }
}

int ttkContourForests::vtkDataSetToStdVector(vtkDataSet* input)
{

  triangulation_ = ttkTriangulation::getTriangulation(input);

  if(!triangulation_)
    return -1;

  varyingMesh_ = false;
  if(triangulation_->isEmpty())
    varyingMesh_ = true;
  if(ttkTriangulation::hasChangedConnectivity(triangulation_, input, this))
    varyingMesh_ = true;

  // init
  if (varyingMesh_ || !numberOfVertices_){
    numberOfVertices_ = input->GetNumberOfPoints();
  }

  if(varyingMesh_){
    segmentation_ = input->NewInstance();
    segmentation_->ShallowCopy(input);
  }

  if(!input->GetPointData())
    return -2;

  // scalars
  if (scalarField_.length()) {
     vtkInputScalars_ = input->GetPointData()->GetArray(scalarField_.data());
  } else {
     vtkInputScalars_ = input->GetPointData()->GetArray(FieldId);
  }

  if(!vtkInputScalars_)
    return -2;

  varyingDataValues_=(vtkInputScalars_->GetMTime() > GetMTime());
  if (varyingMesh_ || varyingDataValues_ || !inputScalarsName_->size()) {
    if (input->GetPointData()) {
      int numberOfArrays = input->GetPointData()->GetNumberOfArrays();
      int numberOfScalarArrays{};

      for (int i = 0; i < numberOfArrays; ++i) {
        vtkDataArray* inputArray = input->GetPointData()->GetArray(i);
        if (inputArray) {
          if (inputArray->GetNumberOfTuples() == numberOfVertices_ &&
              inputArray->GetNumberOfComponents() == 1) {
            ++numberOfScalarArrays;
          }
        }
      }

      inputScalars_->resize(numberOfScalarArrays);
      inputScalarsName_->resize(numberOfScalarArrays);

      int k{};
      for (int i = 0; i < numberOfArrays; ++i) {
        vtkDataArray* inputArray = input->GetPointData()->GetArray(i);
        if (inputArray) {
          if (inputArray->GetNumberOfTuples() == numberOfVertices_ &&
              inputArray->GetNumberOfComponents() == 1) {
            (*inputScalars_)[k].resize(numberOfVertices_);
            (*inputScalarsName_)[k] = inputArray->GetName();

            for (unsigned int j = 0; j < numberOfVertices_; ++j) {
              (*inputScalars_)[k][j] = inputArray->GetTuple1(j);
            }

            ++k;
          }
        }
      }
    }
    Modified();
  }

  if (scalarField_.size() == 0) {
     vertexScalars_ = &((*inputScalars_)[FieldId]);
     scalarField_   = (*inputScalarsName_)[FieldId];
  } else {
    for (unsigned int i = 0; i < inputScalarsName_->size(); ++i) {
       if ((*inputScalarsName_)[i] == scalarField_) {
          vertexScalars_ = &((*inputScalars_)[i]);
          FieldId        = i;
       }
    }
  }

  auto   result    = std::minmax_element(vertexScalars_->begin(), vertexScalars_->end());
  double scalarMin = *result.first;
  double scalarMax = *result.second;
  deltaScalar_     = (scalarMax - scalarMin);

  // neighbors
  triangulation_->setWrapper(this);

  // offsets
  if (varyingMesh_ || varyingDataValues_ || !vertexSoSoffsets_->size()) {

    vertexSoSoffsets_->clear();

    if((InputOffsetFieldId != -1)&&(inputOffsetScalarFieldName_.empty())){
      if(input->GetPointData()->GetArray(InputOffsetFieldId)){
        inputOffsetScalarFieldName_ =
          input->GetPointData()->GetArray(InputOffsetFieldId)->GetName();
        useInputOffsetScalarField_ = true;
      }
    }

    if (useInputOffsetScalarField_ and inputOffsetScalarFieldName_.length()) {

      auto offsets =
        input->GetPointData()->GetArray(inputOffsetScalarFieldName_.data());

      if(offsets->GetNumberOfTuples() != numberOfVertices_){
        stringstream msg;
        msg << "[ttkContourForests] Mesh and offset sizes do not match :("
          << endl;
        msg << "[ttkContourForests] Using default offset field instead..."
          << endl;
        dMsg(cerr, msg.str(), Debug::infoMsg);
      }
      else{
        vertexSoSoffsets_->resize(offsets->GetNumberOfTuples());
        for (unsigned int i = 0; i < vertexSoSoffsets_->size(); ++i){
          (*vertexSoSoffsets_)[i] = offsets->GetTuple1(i);
        }
      }
    }
    if(vertexSoSoffsets_->empty()){
      vertexSoSoffsets_->resize(numberOfVertices_);
      for (unsigned int i = 0; i < vertexSoSoffsets_->size(); ++i){
        (*vertexSoSoffsets_)[i] = i;
      }
    }
    toUpdateVertexSoSoffsets_ = false;
  }

  if (varyingMesh_ || varyingDataValues_ || !isLoaded_) {
    stringstream msg;
    msg << "[ttkContourForests] Convenient data storage has been loaded." << endl;
    msg << "[ttkContourForests]   Number of input scalars: " << inputScalars_->size() << endl;
    msg << "[ttkContourForests]   Input scalars name:" << endl;
    for (unsigned int i = 0; i < inputScalarsName_->size(); ++i)
      msg << "[ttkContourForests]     " << (*inputScalarsName_)[i] << endl;
    msg << "[ttkContourForests]   Active scalar name: " << scalarField_ << endl;
    msg << "[ttkContourForests]   Number of tuples: " << vertexScalars_->size() << endl;
    msg << "[ttkContourForests]   [min max]: [" << scalarMin << " " << scalarMax << "]" << endl;
    msg << "[ttkContourForests]   Number of vertices: " << numberOfVertices_ << endl;
    msg << "[ttkContourForests]   Vertex offsets: " << boolalpha << (bool)vertexSoSoffsets_->size()
      << endl;
    dMsg(cout, msg.str(), detailedInfoMsg);
  }

  stringstream msg;
  msg << "[ttkContourForests] Launching computation for field '"
    << scalarField_ << "'..." << endl;
  dMsg(cout, msg.str(), timeMsg);

  isLoaded_ = true;
  return 0;
}

bool ttkContourForests::isCoincident(float p1[], double p2[])
{
  double sPrev[3];
  double sNext[3];
  for (unsigned int k = 0; k < 3; k++) {
    sPrev[k] = p2[k] - p1[k];
    sNext[k] = sPrev[k];
  }

  return (vtkMath::Normalize(sNext) == 0.0);
}

bool ttkContourForests::isCoincident(double p1[], double p2[])
{
  double sPrev[3];
  double sNext[3];
  for (unsigned int k = 0; k < 3; k++) {
    sPrev[k] = p2[k] - p1[k];
    sNext[k] = sPrev[k];
  }

  return (vtkMath::Normalize(sNext) == 0.0);
}

void ttkContourForests::getSkeletonArcs()
{
  vtkSmartPointer<vtkAppendPolyData> app = vtkSmartPointer<vtkAppendPolyData>::New();

  vtkDoubleArray* scalars{};
  vtkIntArray*    identifierScalars{};
  vtkIntArray*    typeScalars{};
  vtkIntArray*    sizeScalars{};
  vtkDoubleArray* spanScalars{};
  int type = static_cast<int>(TreeComponent::Arc);

  float point1[3];
  vector<double> point2(3);
  // get skeleton scalars
  vector<vector<vector<double>>> skeletonScalars(inputScalars_->size());
  for (unsigned int f = 0; f < inputScalars_->size(); ++f)
    getSkeletonScalars((*inputScalars_)[f], skeletonScalars[f]);

  double    inputScalar;
  SuperArc* a;
  int       regionSize;
  double    regionSpan;
  int       currentZone = 0;
  int       regionId;

  for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
    a = tree_->getSuperArc(i);

    if (a->isVisible()) {
      int      upNodeId   = tree_->getSuperArc(i)->getUpNodeId();
      int      upVertex   = tree_->getNode(upNodeId)->getVertexId();
      float    coordUp[3];
      triangulation_->getVertexPoint(
        upVertex, coordUp[0], coordUp[1], coordUp[2]);

      int      downNodeId   = tree_->getSuperArc(i)->getDownNodeId();
      int      downVertex   = tree_->getNode(downNodeId)->getVertexId();
      float    coordDown[3];
      triangulation_->getVertexPoint(
        downVertex, coordDown[0], coordDown[1], coordDown[2]);

      regionSize = tree_->getNumberOfVisibleRegularNode(i);
      regionSpan = Geometry::distance(coordUp, coordDown, 3);
      regionId       = currentZone++;

      // Line //
      if ((*barycenters_)[static_cast<int>(treeType_)][i].size()) {
        // init: min
        int downNodeVId;
        if (treeType_ == TreeType::Split)
          downNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();
        else
          downNodeVId = tree_->getNode(a->getDownNodeId())->getVertexId();

        triangulation_->getVertexPoint(
          downNodeVId, point1[0], point1[1], point1[2]);
        vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();
        line->SetPoint1(point1);

        const auto nbBarycenter = (*barycenters_)[static_cast<int>(treeType_)][i].size();
        for (unsigned int j = 0; j < nbBarycenter; ++j) {
          point2 = (*barycenters_)[static_cast<int>(treeType_)][i][j];
          line->SetPoint2(point2.data());

          if (!isCoincident(point1, point2.data())) {
            line->Update();
            vtkPolyData* lineData = line->GetOutput();

            // Point data //
            for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
              inputScalar = skeletonScalars[f][i][j];

              scalars = vtkDoubleArray::New();
              scalars->SetName((*inputScalarsName_)[f].data());
              for (unsigned int k = 0; k < 2; ++k)
                scalars->InsertTuple1(k, inputScalar);
              lineData->GetPointData()->AddArray(scalars);
              scalars->Delete();
            }

            // Cell data //
            // Identifier
            identifierScalars = vtkIntArray::New();
            identifierScalars->SetName("SegmentationId");
            for (unsigned int k = 0; k < 2; ++k)
              identifierScalars->InsertTuple1(k, regionId);
            lineData->GetCellData()->AddArray(identifierScalars);
            identifierScalars->Delete();
            // Type
            typeScalars = vtkIntArray::New();
            typeScalars->SetName("Type");
            for (unsigned int k = 0; k < 2; ++k)
              typeScalars->InsertTuple1(k, type);
            lineData->GetCellData()->AddArray(typeScalars);
            typeScalars->Delete();
            // Size
            sizeScalars = vtkIntArray::New();
            sizeScalars->SetName("RegionSize");
            for (unsigned int k = 0; k < 2; ++k)
              sizeScalars->InsertTuple1(k, regionSize);
            lineData->GetCellData()->AddArray(sizeScalars);
            sizeScalars->Delete();
            // Span
            spanScalars = vtkDoubleArray::New();
            spanScalars->SetName("RegionSpan");
            for (unsigned int k = 0; k < 2; ++k)
              spanScalars->InsertTuple1(k, regionSpan);
            lineData->GetCellData()->AddArray(spanScalars);
            spanScalars->Delete();

            app->AddInputData(lineData);
          }

          line = vtkSmartPointer<vtkLineSource>::New();
          line->SetPoint1(point2.data());
        }

        // end: max
        int upNodeVId;
        if (treeType_ == TreeType::Split)
          upNodeVId = tree_->getNode(a->getDownNodeId())->getVertexId();
        else
          upNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();

        float pt[3];
        triangulation_->getVertexPoint(upNodeVId, pt[0], pt[1], pt[2]);
        point2[0]=pt[0];
        point2[1]=pt[1];
        point2[2]=pt[2];
        line->SetPoint2(point2.data());

        if (!isCoincident(point1, point2.data())) {
          line->Update();
          vtkPolyData* lineData = line->GetOutput();

          // Point data //
          for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
            inputScalar =
              skeletonScalars[f][i][(*barycenters_)[static_cast<int>(treeType_)][i].size()];

            scalars = vtkDoubleArray::New();
            scalars->SetName((*inputScalarsName_)[f].data());
            for (unsigned int k = 0; k < 2; ++k)
              scalars->InsertTuple1(k, inputScalar);
            lineData->GetPointData()->AddArray(scalars);
            scalars->Delete();
          }

          // Cell data //
          // Identifier
          identifierScalars = vtkIntArray::New();
          identifierScalars->SetName("SegmentationId");
          for (unsigned int k = 0; k < 2; ++k)
            identifierScalars->InsertTuple1(k, regionId);
          lineData->GetCellData()->AddArray(identifierScalars);
          identifierScalars->Delete();
          // Type
          typeScalars = vtkIntArray::New();
          typeScalars->SetName("Type");
          for (unsigned int k = 0; k < 2; ++k)
            typeScalars->InsertTuple1(k, type);
          lineData->GetCellData()->AddArray(typeScalars);
          typeScalars->Delete();
          // Size
          sizeScalars = vtkIntArray::New();
          sizeScalars->SetName("RegionSize");
          for (unsigned int k = 0; k < 2; ++k)
            sizeScalars->InsertTuple1(k, regionSize);
          lineData->GetCellData()->AddArray(sizeScalars);
          sizeScalars->Delete();
          // Span
          spanScalars = vtkDoubleArray::New();
          spanScalars->SetName("RegionSpan");
          for (unsigned int k = 0; k < 2; ++k)
            spanScalars->InsertTuple1(k, regionSpan);
          lineData->GetCellData()->AddArray(spanScalars);
          spanScalars->Delete();

          app->AddInputData(lineData);
        }
      } else {
        vtkSmartPointer<vtkLineSource> line = vtkSmartPointer<vtkLineSource>::New();

        int downNodeVId = tree_->getNode(a->getDownNodeId())->getVertexId();
        triangulation_->getVertexPoint(
          downNodeVId, point1[0], point1[1], point1[2]);
        line->SetPoint1(point1);

        int upNodeVId = tree_->getNode(a->getUpNodeId())->getVertexId();
        float pt[3];
        triangulation_->getVertexPoint(upNodeVId, pt[0], pt[1], pt[2]);
        point2[0]=pt[0];
        point2[1]=pt[1];
        point2[2]=pt[2];
        line->SetPoint2(point2.data());

        if (!isCoincident(point1, point2.data())) {
          line->Update();
          vtkPolyData* lineData = line->GetOutput();

          // Point data //
          for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
            inputScalar = skeletonScalars[f][i][0];

            scalars = vtkDoubleArray::New();
            scalars->SetName((*inputScalarsName_)[f].data());
            for (unsigned int k = 0; k < 2; ++k)
              scalars->InsertTuple1(k, inputScalar);
            lineData->GetPointData()->AddArray(scalars);
            scalars->Delete();
          }

          // Cell data //
          // Identifier
          identifierScalars = vtkIntArray::New();
          identifierScalars->SetName("SegmentationId");
          for (int k = 0; k < 2; ++k)
            identifierScalars->InsertTuple1(k, regionId);
          lineData->GetCellData()->AddArray(identifierScalars);
          identifierScalars->Delete();
          // Type
          typeScalars = vtkIntArray::New();
          typeScalars->SetName("Type");
          for (unsigned int k = 0; k < 2; ++k)
            typeScalars->InsertTuple1(k, type);
          lineData->GetCellData()->AddArray(typeScalars);
          typeScalars->Delete();
          // Size
          sizeScalars = vtkIntArray::New();
          sizeScalars->SetName("RegionSize");
          for (unsigned int k = 0; k < 2; ++k)
            sizeScalars->InsertTuple1(k, regionSize);
          lineData->GetCellData()->AddArray(sizeScalars);
          sizeScalars->Delete();
          // Span
          spanScalars = vtkDoubleArray::New();
          spanScalars->SetName("RegionSpan");
          for (unsigned int k = 0; k < 2; ++k)
            spanScalars->InsertTuple1(k, regionSpan);
          lineData->GetCellData()->AddArray(spanScalars);
          spanScalars->Delete();

          app->AddInputData(lineData);
        }
      }
    } else {
      //cout << " pruned _ :" << tree_->getNode(a->getDownNodeId())->getVertexId() << " -  "
      //<< tree_->getNode(a->getUpNodeId())->getVertexId() << endl;
    }
  }

  app->Update();
  skeletonArcs_->ShallowCopy(app->GetOutput());
}

int ttkContourForests::getSkeletonScalars(const vector<double>& scalars,
    vector<vector<double>>& skeletonScalars) const
{
  skeletonScalars.clear();
  skeletonScalars.resize(tree_->getNumberOfSuperArcs());

  int nodeId;
  int vertexId;

  double f;
  double f0;
  double f1;
  double fmin;
  double fmax;
  int nodeMinId;
  int nodeMaxId;
  int nodeMinVId;
  int nodeMaxVId;
  const SuperArc* a;
  for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
    a = tree_->getSuperArc(i);

    if (!a->isPruned()) {
      if (treeType_ == TreeType::Split) {
        nodeMinId = a->getUpNodeId();
        nodeMaxId = a->getDownNodeId();
      } else {
        nodeMaxId = a->getUpNodeId();
        nodeMinId = a->getDownNodeId();
      }

      nodeMaxVId = tree_->getNode(nodeMaxId)->getVertexId();
      nodeMinVId = tree_->getNode(nodeMinId)->getVertexId();

      fmax = scalars[nodeMaxVId];
      fmin = scalars[nodeMinVId];

      // init: min
      f0 = fmin;

      // iteration
      for (unsigned int j = 0; j < (*samples_)[static_cast<int>(treeType_)][i].size(); ++j) {
        const vector<int>& sample = (*samples_)[static_cast<int>(treeType_)][i][j];

        f = 0;
        for (unsigned int k = 0; k < sample.size(); ++k) {
          nodeId   = sample[k];
          vertexId = nodeId;
          f += scalars[vertexId];
        }
        if (sample.size()) {
          f /= sample.size();

          f1 = f;
          // update the arc
          skeletonScalars[i].push_back((f0 + f1) / 2);
          f0 = f1;
        }
      }

      // end: max
      f1 = fmax;

      // update the arc
      skeletonScalars[i].push_back((f0 + f1) / 2);
    }
  }

  return 0;
}

void ttkContourForests::getSkeletonNodes()
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  float point[3];

  double scalar{};
  vector<vtkDoubleArray*> scalars(inputScalars_->size());
  for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
    scalars[f] = vtkDoubleArray::New();
    scalars[f]->SetName((*inputScalarsName_)[f].data());
  }

  vtkIntArray* nodeIdentifierScalars = vtkIntArray::New();
  nodeIdentifierScalars->SetName("NodeIdentifier");

  vtkIntArray* vertexIdentifierScalars = vtkIntArray::New();
  vertexIdentifierScalars->SetName("VertexIdentifier");

  int type{};
  vtkIntArray* nodeTypeScalars = vtkIntArray::New();
  nodeTypeScalars->SetName("NodeType");

  vtkIntArray* regionSizeScalars = vtkIntArray::New();
  regionSizeScalars->SetName("RegionSize");

  int identifier{};
  for (unsigned i = 0; i < criticalPoints_->size(); ++i) {
    int nodeId        = (*criticalPoints_)[i];
    if(tree_->getNode(nodeId)->isHidden()) continue;
    int vertexId      = tree_->getNode(nodeId)->getVertexId();
    NodeType nodeType = getNodeType(nodeId);

    if ((nodeType == NodeType::Local_minimum and showMin_) or
        (nodeType == NodeType::Local_maximum and showMax_) or
        (nodeType == NodeType::Saddle1 and showSaddle1_) or
        (nodeType == NodeType::Saddle2 and showSaddle2_) or
        ((showSaddle1_ and showSaddle2_) and
         (nodeType == NodeType::Regular or nodeType == NodeType::Degenerate))) {
      // Positions
      triangulation_->getVertexPoint(vertexId, point[0], point[1], point[2]);
      points->InsertPoint(identifier, point);

      // Scalars
      for (unsigned int f = 0; f < inputScalars_->size(); ++f) {
        scalar = (*inputScalars_)[f][vertexId];
        scalars[f]->InsertTuple1(identifier, scalar);
      }

      // NodeIdentifier
      nodeIdentifierScalars->InsertTuple1(identifier, nodeId);

      // VertexIdentifier
      vertexIdentifierScalars->InsertTuple1(identifier, vertexId);

      // Type
      type = static_cast<int>(nodeType);
      nodeTypeScalars->InsertTuple1(identifier, type);

      // RegionSize
      int regionSize=0;
      if(nodeType==NodeType::Local_maximum){
        const int arcId=tree_->getNode(nodeId)->getDownSuperArcId(0);
        regionSize=tree_->getSuperArc(arcId)->getNumberOfRegularNodes()+1;
      }
      else if(nodeType==NodeType::Local_minimum){
        const int arcId=tree_->getNode(nodeId)->getUpSuperArcId(0);
        regionSize=tree_->getSuperArc(arcId)->getNumberOfRegularNodes()+1;
      }
      regionSizeScalars->InsertTuple1(identifier, regionSize);

      ++identifier;
    }
  }
  skeletonNodes_->SetPoints(points);
  for (unsigned int f = 0; f < inputScalars_->size(); ++f)
    skeletonNodes_->GetPointData()->AddArray(scalars[f]);
  skeletonNodes_->GetPointData()->AddArray(nodeIdentifierScalars);
  skeletonNodes_->GetPointData()->AddArray(vertexIdentifierScalars);
  skeletonNodes_->GetPointData()->AddArray(nodeTypeScalars);
  skeletonNodes_->GetPointData()->AddArray(regionSizeScalars);

  for (unsigned int f = 0; f < inputScalars_->size(); ++f)
    scalars[f]->Delete();
  nodeIdentifierScalars->Delete();
  vertexIdentifierScalars->Delete();
  nodeTypeScalars->Delete();
  regionSizeScalars->Delete();
}

NodeType ttkContourForests::getNodeType(int id)
{
  return getNodeType(id, treeType_, tree_);
}

NodeType ttkContourForests::getNodeType(int id, TreeType type, MergeTree* tree)
{
  int upDegree{};
  int downDegree{};
  if (type == TreeType::Join || type == TreeType::Contour) {
    upDegree   = tree->getNode(id)->getUpValence();
    downDegree = tree->getNode(id)->getDownValence();
  } else {
    downDegree   = tree->getNode(id)->getUpValence();
    upDegree = tree->getNode(id)->getDownValence();
  }
  int degree = upDegree + downDegree;

  // saddle point
  if (degree > 1) {
    if (upDegree == 2 && downDegree == 1)
      return NodeType::Saddle2;
    else if (upDegree == 1 && downDegree == 2)
      return NodeType::Saddle1;
    else if(upDegree == 1 && downDegree == 1)
      return NodeType::Regular;
    else
      return NodeType::Degenerate;
  }
  // local extremum
  else {
    if (upDegree)
      return NodeType::Local_minimum;
    else
      return NodeType::Local_maximum;
  }
}

void ttkContourForests::getCriticalPoints()
{
  vector<bool> isCriticalPoint(numberOfVertices_);

  criticalPoints_->clear();
  for (unsigned int i = 0; i < numberOfVertices_; ++i)
    isCriticalPoint[i] = false;

  //const int nbVert = triangulation_->getNumberOfOriginalVertices();

  // looking for critical points
  for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
    auto a = tree_->getSuperArc(i);

    if (!a->isPruned()) {
      int upId   = a->getUpNodeId();
      int up_vId = tree_->getNode(upId)->getVertexId();
      if (!isCriticalPoint[up_vId]) {
        isCriticalPoint[up_vId] = true;
        criticalPoints_->push_back(upId);
      }

      int downId   = a->getDownNodeId();
      int down_vId = tree_->getNode(downId)->getVertexId();
      if (!isCriticalPoint[down_vId]) {
        isCriticalPoint[down_vId] = true;
        criticalPoints_->push_back(downId);
      }
    }
  }
  //{
  //stringstream msg;
  //msg << "[ttkContourForests] List of critical points :" << endl;
  //for (unsigned int it = 0; it < criticalPoints_->size(); ++it)
  //msg << "[ttkContourForests]   NodeId:" << (*criticalPoints_)[it]
  //<< ", VertexId:" << tree_->getNode(it)->getVertexId() << endl;
  //dMsg(cout, msg.str(), advancedInfoMsg);
  //}
}

int ttkContourForests::sample(unsigned int samplingLevel)
{
  samples_->resize(3);
  (*samples_)[static_cast<int>(treeType_)].resize(tree_->getNumberOfSuperArcs());
  vector<vector<int>> sampleList(samplingLevel);

  SuperArc* a;
  for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
    a = tree_->getSuperArc(i);

    if (!a->isPruned()) {
      for (unsigned int j = 0; j < samplingLevel; ++j) {
        sampleList[j].clear();
      }

      double fmax, fmin;
      int nodeMaxId, nodeMinId;
      int nodeMaxVId, nodeMinVId;
      double delta;
      if (a->getNumberOfRegularNodes()) {
        if (treeType_ == TreeType::Split) {
          nodeMaxId = a->getDownNodeId();
          nodeMinId = a->getUpNodeId();
        } else {
          nodeMaxId = a->getUpNodeId();
          nodeMinId = a->getDownNodeId();
        }

        nodeMaxVId = tree_->getNode(nodeMaxId)->getVertexId();
        nodeMinVId = tree_->getNode(nodeMinId)->getVertexId();

        fmax = (*vertexScalars_)[nodeMaxVId];
        fmin = (*vertexScalars_)[nodeMinVId];

        delta = (fmax - fmin) / samplingLevel;

        double f;
        int nodeId;
        int vertexId;
        for (int j = 0; j < a->getNumberOfRegularNodes(); ++j) {
          nodeId   = a->getRegularNodeId(j);
          if(a->isMasqued(j)) continue;
          vertexId = nodeId;
          f        = (*vertexScalars_)[vertexId];

          for (unsigned int k = 0; k < samplingLevel; ++k) {
            if (f <= (k + 1) * delta + fmin) {
              sampleList[k].push_back(nodeId);
              break;
            }
          }
        }

        // update the arc
        for (unsigned int j = 0; j < sampleList.size(); ++j)
          (*samples_)[static_cast<int>(treeType_)][i].push_back(sampleList[j]);
      }
    }
  }

  return 0;
}

int ttkContourForests::computeBarycenters()
{
  barycenters_->resize(3);
  (*barycenters_)[static_cast<int>(treeType_)].resize(tree_->getNumberOfSuperArcs());
  vector<float> barycenter(3);
  int vertexId;

  const SuperArc* a;
  for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
    a = tree_->getSuperArc(i);
    if (!a->isPruned()) {
      for (unsigned int j = 0; j < (*samples_)[static_cast<int>(treeType_)][i].size(); ++j) {
        vector<int>& sample = (*samples_)[static_cast<int>(treeType_)][i][j];

        for (unsigned int k = 0; k < 3; ++k)
          barycenter[k] = 0;

        for (unsigned int k = 0; k < sample.size(); ++k) {
          vertexId = sample[k];

          float pt[3];
          triangulation_->getVertexPoint(vertexId, pt[0], pt[1], pt[2]);
          barycenter[0]+=pt[0];
          barycenter[1]+=pt[1];
          barycenter[2]+=pt[2];
        }
        if (sample.size()) {
          for (unsigned int k = 0; k < 3; ++k)
            barycenter[k] /= sample.size();

          // update the arc
          unsigned int nbBar = (*barycenters_)[static_cast<int>(treeType_)][i].size();
          (*barycenters_)[static_cast<int>(treeType_)][i].resize(nbBar + 1);
          (*barycenters_)[static_cast<int>(treeType_)][i][nbBar].resize(3);

          for (unsigned int k = 0; k < 3; ++k)
            (*barycenters_)[static_cast<int>(treeType_)][i][nbBar][k] = barycenter[k];
        }
      }
    }
  }

  return 0;
}

void ttkContourForests::computeSkeleton(unsigned int arcRes)
{
  sample(arcRes);
  computeBarycenters();
}

void ttkContourForests::smoothSkeleton(unsigned int skeletonSmoothing)
{
  for (unsigned int i = 0; i < skeletonSmoothing; i++) {
    for (unsigned int j = 0; j < tree_->getNumberOfSuperArcs(); j++) {
      if (!tree_->getSuperArc(j)->isPruned()) {
        smooth(j, !(treeType_ == TreeType::Split));
      }
    }
  }
}

void ttkContourForests::smooth(const int idArc, bool order)
{
  int N = (*barycenters_)[static_cast<int>(treeType_)][idArc].size();
  if (N) {
    // init //
    vector<vector<double>> barycenterList(N);
    for (unsigned int i = 0; i < barycenterList.size(); ++i)
      barycenterList[i].resize(3);

    int up_vId;
    int down_vId;
    if (order) {
      up_vId   = tree_->getNode(tree_->getSuperArc(idArc)->getUpNodeId())->getVertexId();
      down_vId = tree_->getNode(tree_->getSuperArc(idArc)->getDownNodeId())->getVertexId();
    } else {
      down_vId = tree_->getNode(tree_->getSuperArc(idArc)->getUpNodeId())->getVertexId();
      up_vId   = tree_->getNode(tree_->getSuperArc(idArc)->getDownNodeId())->getVertexId();
    }

    float p0[3];
    float p1[3];
    triangulation_->getVertexPoint(down_vId, p0[0], p0[1], p0[2]);
    triangulation_->getVertexPoint(up_vId, p1[0], p1[1], p1[2]);

    // filtering //
    if (N > 1) {
      // first
      for (unsigned int k = 0; k < 3; ++k)
        barycenterList[0][k] =
          (p0[k] + (*barycenters_)[static_cast<int>(treeType_)][idArc][1][k]) * 0.5;

      // main
      for (int i = 1; i < N - 1; ++i) {
        for (unsigned int k = 0; k < 3; ++k)
          barycenterList[i][k] =
            ((*barycenters_)[static_cast<int>(treeType_)][idArc][i - 1][k] +
             (*barycenters_)[static_cast<int>(treeType_)][idArc][i + 1][k]) *
            0.5;
      }
      // last
      for (unsigned int k = 0; k < 3; ++k)
        barycenterList[N - 1][k] =
          (p1[k] + (*barycenters_)[static_cast<int>(treeType_)][idArc][N - 1][k]) * 0.5;
    } else {
      for (unsigned int k = 0; k < 3; ++k)
        barycenterList[0][k] = (p0[k] + p1[k]) * 0.5;
    }

    // copy //
    for (int i = 0; i < N; ++i) {
      for (unsigned int k = 0; k < 3; ++k)
        (*barycenters_)[static_cast<int>(treeType_)][idArc][i][k] = barycenterList[i][k];
    }
  }
}

void ttkContourForests::getSkeleton()
{
  Timer t;
  computeSkeleton(arcResolution_);
  smoothSkeleton(skeletonSmoothing_);

  // nodes
  if (showMin_ || showMax_ || showSaddle1_ || showSaddle2_)
    getSkeletonNodes();
  else
    skeletonNodes_->ShallowCopy(voidUnstructuredGrid_);

  // arcs
  if (showArc_)
    getSkeletonArcs();
  else
    skeletonArcs_->ShallowCopy(voidPolyData_);

  // ce qui est fait n'est plus à faire
  toComputeSkeleton_ = false;

  {
    stringstream msg;
    msg << "[ttkContourForests] Topological skeleton built in " << t.getElapsedTime()
      << "s :" << endl;
    msg << "[ttkContourForests]   Arc - Resolution: " << arcResolution_ << endl;
    msg << "[ttkContourForests]   Smoothing: " << boolalpha << skeletonSmoothing_ << endl;
    dMsg(cout, msg.str(), timeMsg);
  }
}

void ttkContourForests::getSegmentation(vtkDataSet* input)
{
  Timer t;

  // field
  int regionId{};
  vtkSmartPointer<vtkIntArray> scalarsRegionId = vtkSmartPointer<vtkIntArray>::New();
  scalarsRegionId->SetName("SegmentationId");
  scalarsRegionId->SetNumberOfTuples(vertexScalars_->size());

  int regionType{};
  vtkSmartPointer<vtkIntArray> scalarsRegionType = vtkSmartPointer<vtkIntArray>::New();
  scalarsRegionType->SetName("RegionType");
  scalarsRegionType->SetNumberOfTuples(vertexScalars_->size());

  int regionSize{};
  vtkSmartPointer<vtkIntArray> scalarsRegionSize = vtkSmartPointer<vtkIntArray>::New();
  scalarsRegionSize->SetName("RegionSize");
  scalarsRegionSize->SetNumberOfTuples(vertexScalars_->size());

  double regionSpan{};
  vtkSmartPointer<vtkDoubleArray> scalarsRegionSpan = vtkSmartPointer<vtkDoubleArray>::New();
  scalarsRegionSpan->SetName("RegionSpan");
  scalarsRegionSpan->SetNumberOfTuples(vertexScalars_->size());

  int currentZone{};

  if (!segmentation_) {
    segmentation_ = input->NewInstance();
    segmentation_->ShallowCopy(input);
  }

  for (unsigned int i = 0; i < numberOfVertices_; i++) {
    scalarsRegionId->SetTuple1(i, -1);
  }

  // nodes
  for (unsigned int it = 0; it < criticalPoints_->size(); ++it) {
    int nodeId   = (*criticalPoints_)[it];
    int vertexId = tree_->getNode(nodeId)->getVertexId();

    // RegionType
    regionType = -1;
    scalarsRegionType->SetTuple1(vertexId, regionType);
  }

  // arcs
  for (unsigned int i = 0; i < tree_->getNumberOfSuperArcs(); ++i) {
    auto a = tree_->getSuperArc(i);
    if (a->isVisible()) {
      int      upNodeId   = tree_->getSuperArc(i)->getUpNodeId();
      NodeType upNodeType = getNodeType(upNodeId);
      int      upVertex   = tree_->getNode(upNodeId)->getVertexId();
      float    coordUp[3];
      triangulation_->getVertexPoint(
        upVertex, coordUp[0], coordUp[1], coordUp[2]);

      int      downNodeId   = tree_->getSuperArc(i)->getDownNodeId();
      NodeType downNodeType = getNodeType(downNodeId);
      int      downVertex   = tree_->getNode(downNodeId)->getVertexId();
      float    coordDown[3];
      triangulation_->getVertexPoint(
        downVertex, coordDown[0], coordDown[1], coordDown[2]);

      regionSize = tree_->getNumberOfVisibleRegularNode(i);
      regionSpan = Geometry::distance(coordUp, coordDown);
      regionId   = currentZone++;
      //regionId = i;

      //cout << "arc : " << tree_->printArc(i);
      //cout << " span : " << regionSpan;
      //cout << " coords : ";
      //cout << coordDown[0] << ",";
      //cout << coordDown[1] << ",";
      //cout << coordDown[2] << " || ";
      //cout << coordUp[0] << ",";
      //cout << coordUp[1] << ",";
      //cout << coordUp[2] << endl;

      scalarsRegionId->SetTuple1(tree_->getNode(downNodeId)->getVertexId(), regionId);
      scalarsRegionId->SetTuple1(tree_->getNode(upNodeId)->getVertexId(), regionId);

      scalarsRegionSize->SetTuple1(tree_->getNode(downNodeId)->getVertexId(), regionSize);
      scalarsRegionSize->SetTuple1(tree_->getNode(upNodeId)->getVertexId(), regionSize);

      scalarsRegionSpan->SetTuple1(tree_->getNode(downNodeId)->getVertexId(), regionSpan);
      scalarsRegionSpan->SetTuple1(tree_->getNode(upNodeId)->getVertexId(), regionSpan);

      for (int j = 0; j < tree_->getSuperArc(i)->getNumberOfRegularNodes(); ++j) {
        int nodeId   = tree_->getSuperArc(i)->getRegularNodeId(j);
        int vertexId = nodeId;
        // cout << vertexId << ", ";
        if(tree_->getSuperArc(i)->isMasqued(j)) {
          // cout << vertexId << ", ";
          continue;
        }

        //cout << vertexId << ", ";
        scalarsRegionId->SetTuple1(vertexId, regionId);
        scalarsRegionSize->SetTuple1(vertexId, regionSize);
        scalarsRegionSpan->SetTuple1(vertexId, regionSpan);
      }
      //cout << endl;

      // RegionType
      if (upNodeType == NodeType::Local_minimum && downNodeType == NodeType::Local_maximum)
        regionType = static_cast<int>(ArcType::Min_arc);
      else if (upNodeType == NodeType::Local_minimum || downNodeType == NodeType::Local_minimum)
        regionType = static_cast<int>(ArcType::Min_arc);
      else if (upNodeType == NodeType::Local_maximum || downNodeType == NodeType::Local_maximum)
        regionType = static_cast<int>(ArcType::Max_arc);
      else if (upNodeType == NodeType::Saddle1 && downNodeType == NodeType::Saddle1)
        regionType = static_cast<int>(ArcType::Saddle1_arc);
      else if (upNodeType == NodeType::Saddle2 && downNodeType == NodeType::Saddle2)
        regionType = static_cast<int>(ArcType::Saddle2_arc);
      else
        regionType = static_cast<int>(ArcType::Saddle1_saddle2_arc);

      for (int j = 0; j < tree_->getSuperArc(i)->getNumberOfRegularNodes(); ++j) {
        int nodeId   = tree_->getSuperArc(i)->getRegularNodeId(j);
        if(tree_->getSuperArc(i)->isMasqued(j)){
          // Ignore masqued ones
          continue;
        }
        int vertexId = nodeId;
        scalarsRegionType->SetTuple1(vertexId, regionType);
      }
    }
  }

  // output
  segmentation_->GetPointData()->AddArray(scalarsRegionId);
  segmentation_->GetPointData()->AddArray(scalarsRegionType);
  segmentation_->GetPointData()->AddArray(scalarsRegionSize);
  segmentation_->GetPointData()->AddArray(scalarsRegionSpan);

  {
    stringstream msg;
    msg << "[ttkContourForests] Topological segmentation built in " << t.getElapsedTime()
      << "s :" << endl;
    msg << "[ttkContourForests]   RegionType: " << boolalpha
      << (bool)scalarsRegionType->GetNumberOfTuples() << endl;
    msg << "[ttkContourForests]   SegmentationId: " << boolalpha
      << (bool)scalarsRegionId->GetNumberOfTuples() << endl;
    dMsg(cout, msg.str(), timeMsg);
  }

  toComputeSegmentation_=false;
}

void ttkContourForests::getTree()
{
  setDebugLevel(debugLevel_);
  // sequential params
  contourTree_->setDebugLevel(debugLevel_);
  contourTree_->setupTriangulation(triangulation_);
  contourTree_->setVertexScalars(
      vtkInputScalars_->GetVoidPointer(0));
  if(vertexSoSoffsets_){
    contourTree_->setVertexSoSoffsets(*vertexSoSoffsets_);
  }
  contourTree_->setTreeType(treeType_);
  // parallel params
  contourTree_->setLessPartition(lessPartition_);
  contourTree_->setThreadNumber(threadNumber_);
  contourTree_->setPartitionNum(partitionNum_);
  // simplification params
  contourTree_->setSimplificationMethod(simplificationType_);
  contourTree_->setSimplificationThreshold(simplificationThreshold_);
  // build
  switch(vtkInputScalars_->GetDataType()){
    vtkTemplateMacro(({contourTree_->build<VTK_TT>();}));
  }

  // ce qui est fait n'est plus à faire
  toComputeContourTree_    = false;
}

void ttkContourForests::updateTree()
{
  // polymorphic tree
  switch (treeType_) {
    case TreeType::Join:
      tree_ = contourTree_->getJoinTree();
      break;
    case TreeType::Split:
      tree_ = contourTree_->getSplitTree();
      break;
    case TreeType::JoinAndSplit:
      tree_ = contourTree_->getJoinTree();
      tree_ = contourTree_->getSplitTree();
      break;
    case TreeType::Contour:
      tree_ = contourTree_;
      break;
  }

  getCriticalPoints();

  toUpdateTree_ = false;
}

int ttkContourForests::doIt(vector<vtkDataSet *> &inputs,
  vector<vtkDataSet *> &outputs){

  Memory m;

  vtkDataSet *input = inputs[0];
  vtkPolyData *outputSkeletonNodes = vtkPolyData::SafeDownCast(outputs[0]);
  vtkPolyData *outputSkeletonArcs = vtkPolyData::SafeDownCast(outputs[1]);
  vtkDataSet *outputSegmentation = outputs[2];

  if(!input->GetNumberOfPoints())
    return 1;

  // conversion
  vtkDataSetToStdVector(input);

  if(simplificationType_ == 0){
    simplificationThreshold_ = simplificationThresholdBuffer_ * deltaScalar_;
  } else if (simplificationType_ == 1) {
     double  coord0[3], coord1[3], spanTotal;
     double* bounds           = input->GetBounds();
     coord0[0]                = bounds[0];
     coord1[0]                = bounds[1];
     coord0[1]                = bounds[2];
     coord1[1]                = bounds[3];
     coord0[2]                = bounds[4];
     coord1[2]                = bounds[5];
     spanTotal                = Geometry::distance(coord0, coord1);
     simplificationThreshold_ = simplificationThresholdBuffer_ * spanTotal;
  } else if (simplificationType_ == 2) {
     simplificationThreshold_ =
         simplificationThresholdBuffer_ * triangulation_->getNumberOfVertices();
  }

  // ContourForestsTree //
  if (varyingMesh_ || varyingDataValues_ ||
      toComputeContourTree_) {
    clearTree();
    getTree();
  }

  // update the trees
  if (varyingMesh_ || varyingDataValues_ || toUpdateTree_)
    updateTree();

  // Skeleton //
  if (varyingMesh_ || varyingDataValues_ || toComputeSkeleton_) {
    clearSkeleton();
    getSkeleton();
  }

  // Segmentation //
  if (varyingMesh_ || varyingDataValues_ || toComputeSegmentation_)
    getSegmentation(input);

  // Output //

  // skeleton
  outputSkeletonNodes->ShallowCopy(skeletonNodes_);
  outputSkeletonArcs->ShallowCopy(skeletonArcs_);
  // segmentation
  outputSegmentation->ShallowCopy(segmentation_);

   {
    stringstream msg;
    msg << "[ttkContourForests] Memory usage: "
      << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

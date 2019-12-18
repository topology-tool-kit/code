#include <regex>
#include <algorithm>
#include <ttkPointDataSelector.h>
#include <vtkSetGet.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPointDataSelector)

  // transmit abort signals
  bool ttkPointDataSelector::needsToAbort() {
  return GetAbortExecute();
}

// transmit progress status
int ttkPointDataSelector::updateProgress(const float &progress) {

  {
    stringstream msg;
    msg << "[ttkPointDataSelector] " << progress * 100 << "% processed...."
        << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  UpdateProgress(progress);
  return 0;
}

int ttkPointDataSelector::doIt(vtkDataSet *input, vtkDataSet *output) {
  Memory m;

  output->ShallowCopy(input);

  vtkPointData *inputPointData = input->GetPointData();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputPointData) {
    cerr << "[ttkPointDataSelector] Error: input has no point data." << endl;
    return -1;
  }
#endif

  vtkSmartPointer<vtkPointData> outputPointData
    = vtkSmartPointer<vtkPointData>::New();
#ifndef TTK_ENABLE_KAMIKAZE
  if(!outputPointData) {
    cerr
      << "[ttkPointDataSelector] Error: vtkPointData memory allocation problem."
      << endl;
    return -1;
  }
#endif

  try {
    const int lastField = std::min((int)ScalarFields.size(), RangeId[1] + 1);
    for(int i = RangeId[0]; i < lastField; i++) {
      auto scalar = ScalarFields[i];
      if(scalar.length() > 0 && regex_match(scalar, regex(RegexpString))) {
        vtkDataArray *arr = inputPointData->GetArray(scalar.data());
        if(arr) {

          if(RenameSelected) {
            if(ScalarFields.size() != 1 && RangeId[1] - RangeId[0] != 0) {
              vtkErrorMacro("Can't rename more than one field.");
              return 0;
            }

            if(localFieldCopy_) {
              localFieldCopy_->Delete();
              localFieldCopy_ = nullptr;
            }

            localFieldCopy_ = arr->NewInstance();

            if(localFieldCopy_) {
              localFieldCopy_->DeepCopy(arr);
              localFieldCopy_->SetName(SelectedFieldName.data());
              arr = localFieldCopy_;
            }
          }

          outputPointData->AddArray(arr);
        }
      }
    }
  } catch(std::regex_error &) {
    vtkWarningMacro("[ttkPointDataSelector]: Bad regexp.");
  }

  output->GetPointData()->ShallowCopy(outputPointData);

  {
    stringstream msg;
    msg << "[ttkPointDataSelector] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 0;
}

int ttkPointDataSelector::RequestInformation(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  NbScalars = input->GetPointData()->GetNumberOfArrays();
  return vtkDataSetAlgorithm::RequestUpdateExtent(request, inputVector, outputVector);
}

int ttkPointDataSelector::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {
  Memory m;

  vtkDataSet *input = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *output = vtkDataSet::GetData(outputVector);

  doIt(input, output);

  {
    stringstream msg;
    msg << "[ttkPointDataSelector] Memory usage: " << m.getElapsedUsage()
        << " MB." << endl;
    dMsg(cout, msg.str(), memoryMsg);
  }

  return 1;
}

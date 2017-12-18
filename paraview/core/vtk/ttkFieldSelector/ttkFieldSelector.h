/// \ingroup vtk
/// \class ttkFieldSelector
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date December 2017
///
/// \brief TTK VTK-filter that selects scalar fields on input with shallow copy.
/// 
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
#pragma once

#ifndef _MSC_VER
// ttk code includes
#include<Wrapper.h>

// VTK includes
#include<vtkCharArray.h>
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>
#else
// VTK includes
#include<vtkCharArray.h>
#include<vtkDataArray.h>
#include<vtkDataSet.h>
#include<vtkDataSetAlgorithm.h>
#include<vtkDoubleArray.h>
#include<vtkFiltersCoreModule.h>
#include<vtkFloatArray.h>
#include<vtkInformation.h>
#include<vtkIntArray.h>
#include<vtkObjectFactory.h>
#include<vtkPointData.h>
#include<vtkSmartPointer.h>

// ttk code includes
#include<Wrapper.h>
#endif

#ifndef TTK_PLUGIN
class VTKFILTERSCORE_EXPORT ttkFieldSelector
#else
class ttkFieldSelector
#endif
: public vtkDataSetAlgorithm, public Wrapper{

  public:

    static ttkFieldSelector* New();
    vtkTypeMacro(ttkFieldSelector, vtkDataSetAlgorithm)

      // default ttk setters
      vtkSetMacro(debugLevel_, int);

    void SetThreads(){
      if(!UseAllCores)
        threadNumber_ = ThreadNumber;
      else{
        threadNumber_ = OsCall::getNumberOfCores();
      }
      Modified();
    }

    void SetThreadNumber(int threadNumber){\
      ThreadNumber = threadNumber;\
        SetThreads();\
    }\
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters

    void SetScalarFields(string s){
      ScalarFields.push_back(s);
      Modified();
    }

    int FillInputPortInformation(int port, vtkInformation *info){
      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
          break;
        default:
          break;
      }

      return 1;
    }

    int FillOutputPortInformation(int port, vtkInformation *info){

      switch(port){
        case 0:
          info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkDataSet");
          break;
        default:
          break;
      }

      return 1;
    }

  protected:

    ttkFieldSelector(){
      UseAllCores = false;

      SetNumberOfInputPorts(1);
      SetNumberOfOutputPorts(1);
    }

    ~ttkFieldSelector(){};

    int RequestData(vtkInformation *request,
        vtkInformationVector **inputVector,
        vtkInformationVector *outputVector);

  private:

    bool UseAllCores;
    int ThreadNumber;
    vector<string> ScalarFields;

    int doIt(vtkDataSet *input, vtkDataSet *output);
    bool needsToAbort();
    int updateProgress(const float &progress);
};

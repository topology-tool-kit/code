/// \ingroup vtkWrappers
/// \class ttkCellDataConverter
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date February 2016
///
/// \brief TTK VTK-filter that converts data types for cell-based scalar fields 
/// (for instance, from double to float).
///
/// \param Input Input cell-based scalar field (vtkDataSet)
/// \param Output Output cell-based scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the 
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a 
/// VTK pipeline.
///
/// \sa vtkPointDataConverter

#ifndef _TTK_CELLDATACONVERTER_H
#define _TTK_CELLDATACONVERTER_H

// ttk code includes
#include                  <Wrapper.h>

// VTK includes -- to adapt
#include                  <vtkUnsignedCharArray.h>
#include                  <vtkUnsignedShortArray.h>
#include                  <vtkDataArray.h>
#include                  <vtkDataSet.h>
#include                  <vtkDataSetAlgorithm.h>
#include                  <vtkDoubleArray.h>
#include                  <vtkFiltersCoreModule.h>
#include                  <vtkFloatArray.h>
#include                  <vtkInformation.h>
#include                  <vtkIntArray.h>
#include                  <vtkObjectFactory.h>
#include                  <vtkCellData.h>
#include                  <vtkSmartPointer.h>

#include<limits>

class VTKFILTERSCORE_EXPORT ttkCellDataConverter 
: public vtkDataSetAlgorithm, public Wrapper{

enum SupportedType{
Double=0,
  Float,
  Int,
  UnsignedShort,
  UnsignedChar,
  };

  public:      
    static ttkCellDataConverter* New();
    
    vtkTypeMacro(ttkCellDataConverter, vtkDataSetAlgorithm);
    
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
    
    void SetThreadNumber(int threadNumber){
      ThreadNumber = threadNumber;
      SetThreads();
    }   
    
    void SetUseAllCores(bool onOff){
      UseAllCores = onOff;
      SetThreads();
    }
    // end of default ttk setters
    
    vtkSetMacro(ScalarField, string);
    vtkGetMacro(ScalarField, string);

    void SetOutputType(int outputType){
      OutputType = outputType;
      Modified();
    }
    vtkGetMacro(OutputType, int);

    void SetUseNormalization(bool onOff){
      UseNormalization = onOff;
      Modified();
    }
    vtkGetMacro(UseNormalization, int);

  protected:
    
    ttkCellDataConverter();    
    ~ttkCellDataConverter();

    int RequestData(vtkInformation *request, 
		      vtkInformationVector **inputVector, vtkInformationVector *outputVector);

    template<typename A, typename B, typename C>
      int convert(vtkDataArray* inputData, vtkDataSet* output);
  private:
    
    bool UseAllCores;
    int ThreadNumber;
    string ScalarField;
    int OutputType;
    bool UseNormalization;
    
    // base code features
    int doIt(vtkDataSet *input, vtkDataSet *output);
    bool needsToAbort();
    int updateProgress(const float &progress);   
};

#endif // _TTK_CELLDATACONVERTER_H

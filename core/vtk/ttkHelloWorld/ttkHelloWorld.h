/// TODO 4: Provide your information
///
/// \ingroup vtk
/// \class ttkHelloWorld
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::HelloWorld module.
///
/// This VTK filter uses the ttk::HelloWorld module to compute the bounding box
/// of a vtkDataSet, which is returned as a vtkUnstructuredGrid.
///
/// \param Input vtkDataSet whose bounding box will be computed.
/// \param Output vtkUnstructuredGrid that corresponds to bounding box of the
/// input.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::HelloWorld
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkHelloWorldModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <ttkTriangulation.h>
#include <ttkScalarWorker.h>

// TTK Base Includes
#include <HelloWorld.h>

class TTKHELLOWORLD_EXPORT ttkHelloWorld
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::HelloWorld // and we inherit from the base class
{
private:
  /**
   * TODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputArrayName{"AveragedScalarField"};

public:
  /**
   * TODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */
  vtkSetMacro(OutputArrayName, std::string);
  vtkGetMacro(OutputArrayName, std::string);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkHelloWorld *New();
  vtkTypeMacro(ttkHelloWorld, ttkAlgorithm);

protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkHelloWorld();
  ~ttkHelloWorld() override;

  /**
   * TODO 8: Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 9: Specify the data object type of each output port
   *         (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  /**
   * Here you call the basecode
   * You need to implement a templated Compute
   * method that recieve the input and output arrays.
   * You can need to use GetRange to transform these vtk
   * arrays into generic range for TTK.
   * You can also access the triangulation here as well
   * as the baseClass.
   */
  struct HelloWorldWorker
    : public ScalarWorker<HelloWorldWorker, ttkHelloWorld> {
    using ScalarWorker::GetRange;
    using ScalarWorker::ScalarWorker;
    using ScalarWorker::operator();

    template <typename ArrayType1, typename ArrayType2>
    void Compute(ArrayType1 *input, ArrayType2 *output) {
      auto inputRange = GetRange(input);
      auto outputRange = GetRange(output);
      baseClass->computeAverages(outputRange, inputRange, triangulation->getData());
    }
  };
};

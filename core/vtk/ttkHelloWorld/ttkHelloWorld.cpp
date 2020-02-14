#include <ttkHelloWorld.h>

#include <vtkArrayDispatch.txx>
#include <vtkDataObject.h> // For port information
#include <vtkObjectFactory.h> // for new macro

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include "HelloWorld.h"
#include "Triangulation.h"

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkHelloWorld);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkHelloWorld::ttkHelloWorld() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkHelloWorld::~ttkHelloWorld() {
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkHelloWorld::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  else
    return 0;

  return 1;
}

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( ttkAlgorithm::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkHelloWorld::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;

  return 1;
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkHelloWorld::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkDataSet as required by FillInputPortInformation
  vtkDataSet *inputDataSet = vtkDataSet::GetData(inputVector[0]);

  // Get input array that will be processed
  //
  // Note: VTK provides abstract functionality to handle array selections, but
  //       this essential functionality is unfortunately not well documented.
  //       Before you read further, please keep in mind the the TTK developer
  //       team is not responsible for the existing VTK Api ;-)
  //
  //       In a nutshell, prior to the RequestData execution one has to call
  //
  //           SetInputArrayToProcess (
  //               int idx,
  //               int port,
  //               int connection,
  //               int fieldAssociation,
  //               const char *name
  //            )
  //
  //       The parameter 'idx' is often missunderstood: lets say the filter
  //       requires n arrays, then idx enumerates them from 0 to n-1.
  //
  //       The 'port' is the input port index at which the object is connected
  //       from which we want to get the array.
  //
  //       The 'connection' is the connection index at that port (we have to
  //       specify this because VTK allows multiple connections at the same
  //       input port).
  //
  //       The 'fieldAssociation' integer specifies if the array should be taken
  //       from 0: point data, 1: cell data, or 2: field data.
  //
  //       The final parameter is the 'name' of the array.
  //
  //       Example: SetInputArrayToProcess(3,1,0,1,"EdgeLength") will store for
  //                the 3rd array the filter needs the cell data array named
  //                "EdgeLength" that it will retrieve from the vtkDataObject
  //                at input port 1 (first connection). During the RequestData
  //                method one can then actually retrieve the 3rd array it
  //                requires for its computation by calling
  //                GetInputArrayToProcess(3, inputVector)
  //
  //       If this filter is run within ParaView, then the UI will automatically
  //       call SetInputArrayToProcess (see HelloWorld.xml file).
  //
  vtkDataArray *inputArray = this->GetInputArrayToProcess(0, inputVector);

  // Create an output array that has the same data type as the input array
  // Note: vtkSmartPointers are well documented
  //       (https://vtk.org/Wiki/VTK/Tutorials/SmartPointers)
  vtkSmartPointer<vtkDataArray> outputArray
    = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputArray->SetName(this->OutputArrayName.data()); // set array name
  outputArray->SetNumberOfComponents(1); // only one component per tuple
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  // Get ttk::triangulation of the input vtkDataSet (will create one if does
  // not exist already)
  ttk::Triangulation *triangulation
    = ttkAlgorithm::GetTriangulation(inputDataSet);

  // Precondition the triangulation (e.g., enable fetching of vertex neighbors)
  this->preconditionTriangulation(triangulation); // implemented in base class

  // Call the base code computation using the right vtk Dispatch
  HelloWorldWorker scalarWorker(this, triangulation);
  if(!vtkArrayDispatch::Dispatch2BySameValueType<
       vtkArrayDispatch::AllTypes>::Execute(inputArray, outputArray,
                                            scalarWorker)) {
    std::cerr << "Bad array type" << std::endl;
    return 0;
  }
  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkDataSet::GetData(outputVector, 0);

  // make a SHALLOW copy of the input
  outputDataSet->ShallowCopy(inputDataSet);

  // add to the output point data the computed output array
  outputDataSet->GetPointData()->AddArray(outputArray);

  return 1;
}

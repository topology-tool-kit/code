#pragma once

#include "ArrayIterator.h"
#include "Triangulation.h"

// vtk array manipulation
#include "ArrayIterator.h"
#include "vtkArrayDispatch.h"
#include "vtkDataArrayRange.h"

template <typename Wrapper, typename Algo>
class ScalarWorker {

public:
  ScalarWorker(Algo *a, ttk::Triangulation *t)
    : baseClass(a), triangulation(t) {
  }

  // Helper to transform VTK input into generic range
  template <template <typename> class ArrayType, typename ElType>
  RangeHandler<ElType> GetRange(ArrayType<ElType> *arr) {
    return RangeHandler<ElType>(vtk::DataArrayValueRange<1>(arr));
  }

  Wrapper* getSelf(){
    return reinterpret_cast<Wrapper*>(this);
  }

  // Define the required operator()
  template <typename... Args>
  auto operator()(Args &&... args) -> void {
    return getSelf()->Compute(std::forward<Args>(args)...);
  }

  // User define method should look like this
  // template<typename ArrayType1, typename ArrayType2>
  // void Compute(ArrayType1 *input, ArrayType2 *output)
  // {
  //   auto inputRange = GetRange(input);
  //   auto outputRange = GetRange(output);
  //   baseClass->compute(inputRange, outputRange, triangulation);
  // }

protected:
  ttk::Triangulation *triangulation;
  Algo *baseClass;
};

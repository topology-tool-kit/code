/// \ingroup base
/// \class ttk::ArrayIterator
/// \author Charles Gueunet charles.gueunet@kitware.com
/// \date 2020-01-08
///
///\brief TTK class to iterate over any array with a generic interface

#pragma once

#include "DataTypes.h"
#include <iterator>

template <typename ElementType>
class GenericIterator
  : public std::iterator<std::random_access_iterator_tag, ElementType> {
  ElementType *p;

public:
  GenericIterator() : p(nullptr) {
  }
  GenericIterator(ElementType *x) : p(x) {
  }
  GenericIterator(const GenericIterator<ElementType> &it) : p(it.p) {
  }
  GenericIterator &operator++() {
    ++p;
    return *this;
  }
  GenericIterator operator++(int) {
    GenericIterator tmp(*this);
    operator++();
    return tmp;
  }
  GenericIterator &operator--() {
    --p;
    return *this;
  }
  GenericIterator operator--(int) {
    GenericIterator tmp(*this);
    operator--();
    return tmp;
  }
  GenericIterator operator+(int offset) const {
    return this->p + offset;
  }
  GenericIterator &operator=(const GenericIterator<ElementType>& it)
  {
    this->p = it->p;
    return *this;
  }
  bool operator==(const GenericIterator<ElementType> &rhs) const {
    return this->p == rhs.p;
  }
  bool operator!=(const GenericIterator<ElementType> &rhs) const {
    return this->p != rhs.p;
  }
  ElementType &operator*() {
    return *this->p;
  }
  const ElementType &operator*() const {
    return *this->p;
  }
};

template <typename ElementType>
class RangeHandler : public GenericIterator<ElementType> {

  using iteratorT = GenericIterator<ElementType>;
  iteratorT beginIt, endIt;

public:
  RangeHandler(iteratorT b, iteratorT e) : beginIt(b), endIt(e) {
  }
  RangeHandler(ElementType *b, ElementType *e) : beginIt(b), endIt(e) {
  }
  template<typename RangeType>
  RangeHandler(RangeType range) : beginIt(std::begin(range)), endIt(std::end(range)) {
  }

  iteratorT begin()
  {
    return this->beginIt;
  }
  iteratorT end()
  {
    return this->endIt;
  }
  const iteratorT& cbegin() const
  {
    return this->beginIt;
  }
  const iteratorT& cend() const
  {
    return this->endIt;
  }
  // TODO: deal with reverse iterators...

  // Access underlying data

  ElementType &operator[](int idx) {
    return *(this->beginIt + idx);
  }

  const ElementType &operator[](int idx) const {
    return *(this->beginIt + idx);
  }
};

/// \ingroup base
/// \class FTMTreeUtils_Template
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEUTILS_TEMPLATE_H
#define _FTMTREEUTILS_TEMPLATE_H

#pragma once

//#include <FTMTreeMT.h>

namespace ttk {
  namespace ftm {
    // --------------------
    // Get
    // --------------------

    template <class dataType>
    std::tuple<dataType, dataType> FTMTree_MT::getBirthDeath(idNode nodeId) {
      idNode originId = this->getNode(nodeId)->getOrigin();
      if(this->isNodeOriginDefined(
           nodeId)) { // Avoid error if origin is not defined
        dataType pers1 = this->getValue<dataType>(nodeId);
        dataType pers2 = this->getValue<dataType>(originId);
        dataType birth = std::min(pers1, pers2);
        dataType death = std::max(pers1, pers2);
        return std::make_tuple(birth, death);
      }
      return std::make_tuple(0.0, 0.0);
    }

    template <class dataType>
    dataType FTMTree_MT::getBirth(idNode nodeId) {
      return std::get<0>(this->getBirthDeath<dataType>(nodeId));
    }

    template <class dataType>
    dataType FTMTree_MT::getNodePersistence(idNode nodeId) {
      std::tuple<dataType, dataType> birthDeath
        = this->getBirthDeath<dataType>(nodeId);
      return std::get<1>(birthDeath) - std::get<0>(birthDeath);
    }

    // --------------------
    // Is
    // --------------------

    template <class dataType>
    bool FTMTree_MT::isJoinTree() {
      auto root = this->getRoot();
      idNode child = this->getChildren(root)[0];
      if(this->isFullMerge()) {
        dataType min = std::numeric_limits<dataType>::max();
        for(unsigned int i = 0; i < this->getNumberOfNodes(); ++i) {
          dataType value = this->getValue<dataType>(i);
          if(not this->isNodeAlone(i) and value < min) {
            min = value;
            child = i;
          }
        }
      }
      return this->getValue<dataType>(root) > this->getValue<dataType>(child);
    }

    template <class dataType>
    bool FTMTree_MT::isImportantPair(idNode nodeId, double threshold) {
      dataType rootPers = this->getNodePersistence<dataType>(this->getRoot());
      if(threshold > 1)
        threshold /= 100.0;
      threshold = rootPers * threshold;
      return this->getNodePersistence<dataType>(nodeId) > threshold;
    }

    // --------------------
    // Utils
    // --------------------

    template <typename type>
    static type myAbs(const type var) {
      return (var >= 0) ? var : -var;
    }

    template <class dataType>
    bool isEqual(dataType first, dataType two, double eps = 1e-6) {
      return myAbs<dataType>(first - two)
             < eps * std::max(myAbs<dataType>(first), myAbs<dataType>(two));
    }

  } // namespace ftm
} // namespace ttk

#endif

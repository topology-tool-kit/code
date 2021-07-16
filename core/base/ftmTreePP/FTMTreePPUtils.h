/// \ingroup base
/// \class FTMTreePPUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _FTMTREEPPUTILS_H
#define _FTMTREEPPUTILS_H

#pragma once

#include <FTMTree.h>
#include <FTMTreePP.h>

namespace ttk {
  namespace ftm {

    template <class dataType>
    void getPersistencePairs(
      FTMTree_MT *tree,
      std::vector<std::tuple<SimplexId, SimplexId, dataType>> &pairs) {
      FTMTreePP pairsCompute;
      pairsCompute.setCustomTree(tree);
      pairsCompute.computePersistencePairs<dataType>(
        pairs, tree->isJoinTree<dataType>());
    }

    template <class dataType>
    std::vector<std::tuple<SimplexId, SimplexId, dataType>>
      computePersistencePairs(FTMTree_MT *tree) {
      std::vector<std::tuple<SimplexId, SimplexId, dataType>> pairs;
      getPersistencePairs<dataType>(tree, pairs);
      for(auto pair : pairs) {
        if(tree->getNode(std::get<0>(pair))->getOrigin() < std::get<0>(pair)
           and tree->getNode(std::get<0>(pair))->getOrigin() >= 0)
          tree->getNode(tree->getNode(std::get<0>(pair))->getOrigin())
            ->setOrigin(std::get<1>(pair));

        tree->getNode(std::get<0>(pair))->setOrigin(std::get<1>(pair));
        tree->getNode(std::get<1>(pair))->setOrigin(std::get<0>(pair));
      }
      return pairs;
    }

  } // namespace ftm
} // namespace ttk

#endif

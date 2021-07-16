/// \ingroup base
/// \class MergeTreeVisualization
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Visualization module for merge trees.
///

#ifndef _MERGETREEVISUALIZATION_H
#define _MERGETREEVISUALIZATION_H

#pragma once

#include <FTMTree.h>

#include <stack>

using namespace ttk;
using namespace ftm;

class MergeTreeVisualization : virtual public Debug {
protected:
  // Visualization parameters
  bool PlanarLayout = false;
  bool BranchDecompositionPlanarLayout = false;
  double BranchSpacing = 1.;
  bool RescaleTreesIndividually = false;
  double ImportantPairs = 50.; // important pairs threshold
  double ImportantPairsSpacing = 1.;
  double NonImportantPairsSpacing = 1.;
  double NonImportantPairsProximity = 0.05;

public:
  MergeTreeVisualization() = default;
  virtual ~MergeTreeVisualization() = default;

  // ==========================================================================
  // Getter / Setter
  // ==========================================================================
  // Visualization parameters
  void setPlanarLayout(bool b) {
    PlanarLayout = b;
  }
  void setBranchDecompositionPlanarLayout(bool b) {
    BranchDecompositionPlanarLayout = b;
  }
  void setBranchSpacing(double d) {
    BranchSpacing = d;
  }
  void setRescaleTreesIndividually(bool b) {
    RescaleTreesIndividually = b;
  }
  void setImportantPairs(double d) {
    ImportantPairs = d;
  }
  void setImportantPairsSpacing(double d) {
    ImportantPairsSpacing = d;
  }
  void setNonImportantPairsSpacing(double d) {
    NonImportantPairsSpacing = d;
  }
  void setNonImportantPairsProximity(double d) {
    NonImportantPairsProximity = d;
  }

  // ==========================================================================
  // Planar Layout
  // ==========================================================================
  // TODO manage multi pers pairs
  template <class dataType>
  void treePlanarLayoutBDImpl(
    FTMTree_MT *tree,
    std::vector<float> &retVec,
    std::vector<LongSimplexId> &treeSimplexId,
    std::vector<idNode> &branching,
    std::vector<std::vector<ftm::idNode>> &nodeBranching) {
    idNode treeRoot = tree->getRoot();
    idNode treeRootOrigin = tree->getNode(treeRoot)->getOrigin();
    float rootY = retVec[treeSimplexId[treeRoot] * 2 + 1];
    float rootOriginY = retVec[treeSimplexId[treeRootOrigin] * 2 + 1];
    float rootYmin = std::min(rootY, rootOriginY);
    float rootYmax = std::max(rootY, rootOriginY);
    dataType rootPers = tree->getNodePersistence<dataType>(treeRoot);
    std::vector<std::tuple<float, float>> allNodeSpanX(
      tree->getNumberOfNodes());
    std::vector<std::tuple<float, float>> allNodeImportantSpanX(
      tree->getNumberOfNodes());

    // Compute gap
    float nonImportantPairsGap
      = (rootYmax - rootYmin) * 0.05 * NonImportantPairsSpacing;
    float importantPairsGap
      = std::max(nonImportantPairsGap, (float)ImportantPairsSpacing);

    // Some functions
    auto compLowerPers = [&](const ftm::idNode a, const ftm::idNode b) {
      return tree->getNodePersistence<dataType>(a)
             < tree->getNodePersistence<dataType>(b);
    };
    auto compEqUpperPers = [&](const ftm::idNode a, const ftm::idNode b) {
      return not compLowerPers(a, b);
    };

    // Go
    std::vector<ftm::idNode> leaves = tree->getLeavesFromTree();
    std::sort(leaves.begin(), leaves.end(), compLowerPers);
    std::queue<ftm::idNode> queue;
    for(auto node : leaves)
      queue.emplace(node);
    while(!queue.empty()) {
      idNode node = queue.front();
      queue.pop();
      idNode nodeOrigin = tree->getNode(node)->getOrigin();

      dataType nodePers = tree->getNodePersistence<dataType>(node);
      retVec[treeSimplexId[nodeOrigin] * 2] = 0;
      retVec[treeSimplexId[nodeOrigin] * 2 + 1]
        = nodePers / rootPers * (rootYmax - rootYmin) + rootYmin;

      // Positioning nodes in the branch
      if(tree->isBranchOrigin(nodeOrigin)) {
        float prevX = 0;
        std::vector<idNode> nodeBranchingVector;
        for(size_t i = 1; i < nodeBranching[nodeOrigin].size(); ++i)
          nodeBranchingVector.push_back(nodeBranching[nodeOrigin][i]);
        std::sort(nodeBranchingVector.begin(), nodeBranchingVector.end(),
                  compEqUpperPers);

        // Iterate through each node of the branch
        int lastIndexImportant = -1;
        for(size_t i = 0; i < nodeBranchingVector.size(); ++i) {
          idNode nodeBranchingI = nodeBranchingVector[i];

          // Get old node span X
          float oldMin = std::get<0>(allNodeSpanX[nodeBranchingI]);
          float oldMax = std::get<1>(allNodeSpanX[nodeBranchingI]);

          // Get x spacing
          float nodeSpacing = 0;
          if(i > 0) {
            if(tree->isImportantPair<dataType>(
                 nodeBranchingVector[i], ImportantPairs)) {
              nodeSpacing = importantPairsGap;
            } else if(tree->isImportantPair<dataType>(
                        nodeBranchingVector[i - 1], ImportantPairs)) {
              nodeSpacing = NonImportantPairsProximity;
              // prevX =
            } else {
              nodeSpacing = nonImportantPairsGap;
            }
          } else if(not tree->isImportantPair<dataType>(
                      nodeBranchingVector[i], ImportantPairs)
                    and tree->isImportantPair<dataType>(
                      nodeOrigin, ImportantPairs))
            nodeSpacing = NonImportantPairsProximity;
          float newMin = prevX + nodeSpacing;
          float shiftX = newMin - oldMin;

          // Set y coordinate according difference in persistence
          dataType nodeBranchingIPers
            = tree->getNodePersistence<dataType>(nodeBranchingI);
          float shiftY = nodeBranchingIPers * BranchSpacing;
          float diffY = retVec[treeSimplexId[nodeBranchingI] * 2 + 1];
          retVec[treeSimplexId[nodeBranchingI] * 2 + 1]
            = retVec[treeSimplexId[nodeOrigin] * 2 + 1] - shiftY;
          diffY = retVec[treeSimplexId[nodeBranchingI] * 2 + 1] - diffY;

          // Shift this branch
          std::queue<idNode> queueBranching;
          queueBranching.emplace(nodeBranchingI);
          while(!queueBranching.empty()) {
            idNode nodeBranchOrigin = queueBranching.front();
            queueBranching.pop();
            retVec[treeSimplexId[nodeBranchOrigin] * 2] += shiftX;
            if(nodeBranchOrigin != nodeBranchingI)
              retVec[treeSimplexId[nodeBranchOrigin] * 2 + 1] += diffY;
            if(tree->isBranchOrigin(nodeBranchOrigin))
              for(auto nodeB : nodeBranching[nodeBranchOrigin])
                queueBranching.emplace(nodeB);
          }

          // Update node span X
          allNodeSpanX[nodeBranchingI]
            = std::make_tuple(oldMin + shiftX, oldMax + shiftX);
          float oldMinImp = std::get<0>(allNodeImportantSpanX[nodeBranchingI]);
          float oldMaxImp = std::get<1>(allNodeImportantSpanX[nodeBranchingI]);
          allNodeImportantSpanX[nodeBranchingI]
            = std::make_tuple(oldMinImp + shiftX, oldMaxImp + shiftX);

          // Update x base for next iteration
          prevX = std::get<1>(allNodeSpanX[nodeBranchingI]);
          if(tree->isImportantPair<dataType>(
               nodeBranchingVector[i], ImportantPairs)) {
            lastIndexImportant = i;
            prevX = std::get<1>(allNodeImportantSpanX[nodeBranchingI]);
            if(i < nodeBranchingVector.size() - 1
               and not tree->isImportantPair<dataType>(
                 nodeBranchingVector[i + 1], ImportantPairs)) {
              float spanMin = std::get<0>(allNodeSpanX[nodeBranchingVector[0]]);
              float spanMax = std::get<1>(
                allNodeImportantSpanX[nodeBranchingVector[lastIndexImportant]]);
              prevX = (spanMin + spanMax) / 2;
            }
          }
        } // end for nodeBranching

        // Update node span X
        float spanMin = std::get<0>(allNodeSpanX[nodeBranchingVector[0]]);
        if(lastIndexImportant != -1) {
          float spanMaxImp = std::get<1>(
            allNodeImportantSpanX[nodeBranchingVector[lastIndexImportant]]);
          allNodeImportantSpanX[nodeOrigin]
            = std::make_tuple(spanMin, spanMaxImp);
        } else {
          allNodeImportantSpanX[nodeOrigin] = std::make_tuple(0, 0);
          spanMin = 0;
        }
        float spanMax = std::get<1>(
          allNodeSpanX[nodeBranchingVector[nodeBranchingVector.size() - 1]]);
        allNodeSpanX[nodeOrigin] = std::make_tuple(spanMin, spanMax);
      } else {
        allNodeSpanX[nodeOrigin] = std::make_tuple(0, 0);
        allNodeImportantSpanX[nodeOrigin] = std::make_tuple(0, 0);
      }

      // Positioning of this node x coordinate
      float spanMin = std::get<0>(allNodeImportantSpanX[nodeOrigin]);
      float spanMax = std::get<1>(allNodeImportantSpanX[nodeOrigin]);
      retVec[treeSimplexId[nodeOrigin] * 2] = (spanMin + spanMax) / 2;
    }

    // Copy coordinates of nodeOrigin
    for(auto node : leaves)
      queue.emplace(node);
    while(!queue.empty()) {
      idNode node = queue.front();
      queue.pop();
      idNode nodeOrigin = tree->getNode(node)->getOrigin();
      retVec[treeSimplexId[node] * 2] = retVec[treeSimplexId[nodeOrigin] * 2];
      retVec[treeSimplexId[node] * 2 + 1]
        = retVec[treeSimplexId[nodeOrigin] * 2 + 1];
    }
  }

  // TODO manage multi pers pairs
  template <class dataType>
  std::vector<float> treePlanarLayoutImpl(
    FTMTree_MT *tree,
    std::tuple<double, double, double, double, double, double> oldBounds,
    double refPersistence) {
    int verbose = 0;

    if(verbose > 0) {
      printMsg(debug::Separator::L1, debug::Priority::VERBOSE);
      printMsg("Planar Layout", debug::Priority::VERBOSE);
    }

    // ----------------------------------------------------
    // Init internal parameters
    // ----------------------------------------------------
    Timer t_init;
    if(verbose > 1)
      printMsg("Init internal parameters", debug::Priority::VERBOSE);

    auto nPoints = tree->getRealNumberOfNodes();
    int outNumberOfPoints = nPoints * 2;
    std::vector<float> retVec(outNumberOfPoints);

    int cptNode = 0;
    std::vector<LongSimplexId> treeSimplexId(tree->getNumberOfNodes());
    std::vector<idNode> branching;
    std::vector<int> branchingID;
    std::vector<std::vector<ftm::idNode>> nodeBranching;
    tree->getTreeBranching(branching, branchingID, nodeBranching);

    if(verbose > 0) {
      std::stringstream ss;
      ss << "INIT            = " << t_init.getElapsedTime();
      printMsg(ss.str(), debug::Priority::VERBOSE);
    }

    // ----------------------------------------------------
    // Iterate through tree
    // ----------------------------------------------------
    Timer t_iterate;
    if(verbose > 1)
      printMsg("Iterate through tree", debug::Priority::VERBOSE);

    std::queue<idNode> queue;
    ftm::idNode treeRoot = tree->getRoot();
    ftm::idNode treeRootOrigin = tree->getNode(treeRoot)->getOrigin();
    queue.emplace(treeRoot);
    while(!queue.empty()) {
      idNode node = queue.front();
      queue.pop();

      // Get and insert point
      treeSimplexId[node] = cptNode;
      ++cptNode;

      // Push children to the queue
      auto children = tree->getChildren(node);
      for(size_t i = 0; i < children.size(); ++i) {
        auto child = children[i];
        queue.emplace(child);
      }
    }

    if(verbose > 0) {
      std::stringstream ss;
      ss << "ITERATE TREE    = " << t_iterate.getElapsedTime();
      printMsg(ss.str(), debug::Priority::VERBOSE);
    }

    // ----------------------------------------------------
    // Prepositioning coordinates
    // ----------------------------------------------------
    std::queue<idNode> queue2;
    queue2.emplace(treeRoot);
    while(!queue2.empty()) {
      idNode node = queue2.front();
      queue2.pop();

      retVec[treeSimplexId[node] * 2]
        = tree->getValue<dataType>(branching[node]);
      retVec[treeSimplexId[node] * 2 + 1] = tree->getValue<dataType>(node);

      // Push children to the queue
      for(auto child : tree->getChildren(node))
        queue2.emplace(child);
    }

    // ----------------------------------------------------
    // Rescale coordinates
    // ----------------------------------------------------
    Timer t_rescale;
    if(verbose > 1)
      printMsg("Rescale coordinates ", debug::Priority::VERBOSE);

    float x_min, y_min, x_max, y_max;
    x_min = std::numeric_limits<float>::max();
    y_min = std::numeric_limits<float>::max();
    x_max = std::numeric_limits<float>::lowest();
    y_max = std::numeric_limits<float>::lowest();
    for(int i = 0; i < outNumberOfPoints; i += 2) {
      x_min = std::min(x_min, retVec[i]);
      x_max = std::max(x_max, retVec[i]);
      y_min = std::min(y_min, retVec[i + 1]);
      y_max = std::max(y_max, retVec[i + 1]);
    }
    auto newBounds = std::make_tuple(x_min, x_max, y_min, y_max, 0, 0);

    // TODO correctly manage diff and offset if RescaleTreesIndividually
    double diff = std::max((std::get<1>(oldBounds) - std::get<0>(oldBounds)),
                           (std::get<3>(oldBounds) - std::get<2>(oldBounds)));
    double offset = std::max(std::get<0>(oldBounds), std::get<2>(oldBounds));
    if(not RescaleTreesIndividually) {
      // diff *= getNodePersistence<dataType>(tree, treeRoot) / refPersistence;
      diff = tree->getNodePersistence<dataType>(treeRoot);
      offset = tree->getBirth<dataType>(treeRoot);
    }

    for(int i = 0; i < outNumberOfPoints; i += 2) {
      auto divisor1
        = std::get<1>(newBounds) - std::get<0>(newBounds); // (x_max - x_min)
      divisor1 = (divisor1 == 0 ? 1 : divisor1);
      auto divisor2
        = std::get<3>(newBounds) - std::get<2>(newBounds); // (y_max - y_min)
      divisor2 = (divisor2 == 0 ? 1 : divisor2);

      // x coordinate
      retVec[i] = (retVec[i] - std::get<0>(newBounds)) / divisor1;
      retVec[i] = retVec[i] * diff / 2 + offset;

      // y coordinate
      retVec[i + 1] = (retVec[i + 1] - std::get<2>(newBounds)) / divisor2;
      retVec[i + 1] = retVec[i + 1] * diff + offset;
    }

    if(verbose > 0) {
      std::stringstream ss;
      ss << "RESCALE COORD.  = " << t_rescale.getElapsedTime();
      printMsg(ss.str(), debug::Priority::VERBOSE);
    }

    // ----------------------------------------------------
    // Call Branch Decomposition Planar Layout if asked
    // ----------------------------------------------------
    if(BranchDecompositionPlanarLayout) {
      treePlanarLayoutBDImpl<dataType>(
        tree, retVec, treeSimplexId, branching, nodeBranching);
      return retVec;
    }

    // ----------------------------------------------------
    // Move nodes given scalars
    // ----------------------------------------------------
    Timer t_move;
    if(verbose > 1)
      printMsg("Move nodes given scalars", debug::Priority::VERBOSE);

    float rootY = retVec[treeSimplexId[treeRoot] * 2 + 1];
    float rootOriginY = retVec[treeSimplexId[treeRootOrigin] * 2 + 1];
    float rootYmin = std::min(rootY, rootOriginY);
    float rootYmax = std::max(rootY, rootOriginY);
    auto rootBirthDeath = tree->getBirthDeath<dataType>(treeRoot);
    dataType rootBirth = std::get<0>(rootBirthDeath);
    dataType rootDeath = std::get<1>(rootBirthDeath);
    for(size_t i = 0; i < tree->getNumberOfNodes(); ++i) {
      retVec[treeSimplexId[i] * 2 + 1]
        = (tree->getValue<dataType>(i) - rootBirth) / (rootDeath - rootBirth);
      retVec[treeSimplexId[i] * 2 + 1]
        = retVec[treeSimplexId[i] * 2 + 1] * (rootYmax - rootYmin) + rootYmin;
    }

    if(verbose > 0) {
      std::stringstream ss;
      ss << "MOVE SCALAR     = " << t_move.getElapsedTime();
      printMsg(ss.str(), debug::Priority::VERBOSE);
    }

    // ----------------------------------------------------
    // Scale pairs given persistence
    // ----------------------------------------------------
    Timer t_scale;
    if(verbose > 1)
      printMsg("Scale pairs given persistence", debug::Priority::VERBOSE);

    dataType rootPers = tree->getNodePersistence<dataType>(treeRoot);

    std::vector<ftm::idNode> leaves = tree->getLeavesFromTree();
    auto compLowerPers = [&](const ftm::idNode a, const ftm::idNode b) {
      return tree->getNodePersistence<dataType>(a)
             < tree->getNodePersistence<dataType>(b);
    };
    std::sort(leaves.begin(), leaves.end(), compLowerPers);
    std::stack<ftm::idNode> stack;
    for(auto node : leaves)
      stack.emplace(node);
    std::vector<bool> nodeDone(tree->getNumberOfNodes(), false);
    while(!stack.empty()) {
      ftm::idNode node = stack.top();
      stack.pop();
      nodeDone[node] = true;

      if(node == treeRoot or node == treeRootOrigin or tree->isNodeAlone(node))
        continue;

      dataType nodePers = tree->getNodePersistence<dataType>(node);
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();

      // Manage leaf
      if(tree->isLeaf(node)) {
        float nodeDiff = (retVec[treeSimplexId[node] * 2]
                          - retVec[treeSimplexId[nodeOrigin] * 2]);
        int sign = nodeDiff / std::abs(nodeDiff);
        auto inc = sign * nodePers / rootPers * (rootYmax - rootYmin) / 2;
        retVec[treeSimplexId[node] * 2]
          = retVec[treeSimplexId[nodeOrigin] * 2] + inc;

        // Push nodes in the branch to the stack
        ftm::idNode nodeParent = tree->getParentSafe(node);
        ftm::idNode oldNodeParent = -1;
        while(nodeParent != nodeOrigin) {
          if(not nodeDone[nodeParent])
            stack.emplace(nodeParent);
          else
            break;
          oldNodeParent = nodeParent;
          nodeParent = tree->getParentSafe(nodeParent);
          if(oldNodeParent == nodeParent) {
            std::stringstream ss;
            ss << "treePlanarLayoutImpl oldNodeParent == nodeParent";
            printMsg(ss.str(), debug::Priority::VERBOSE);
            break;
          }
        }
      }

      // Manage saddle
      if(not tree->isLeaf(node) and not tree->isRoot(node)) {
        float branchY
          = retVec[treeSimplexId[tree->getNode(branching[node])->getOrigin()]
                   * 2];
        retVec[treeSimplexId[node] * 2] = branchY;
      }
    }

    if(verbose > 0) {
      std::stringstream ss;
      ss << "SCALE PERS.     = " << t_scale.getElapsedTime();
      printMsg(ss.str(), debug::Priority::VERBOSE);
    }

    // ----------------------------------------------------
    // Branches positionning and avoid edges crossing
    // ----------------------------------------------------
    Timer t_avoid;
    if(verbose > 1)
      printMsg("Avoid edges crossing", debug::Priority::VERBOSE);

    bool isJT = tree->isJoinTree<dataType>();
    auto compValue = [&](const ftm::idNode a, const ftm::idNode b) {
      return (isJT ? tree->getValue<dataType>(a) < tree->getValue<dataType>(b)
                   : tree->getValue<dataType>(a) > tree->getValue<dataType>(b));
    };
    std::vector<std::tuple<float, float, float, float>> allBranchBounds(
      tree->getNumberOfNodes());
    std::vector<std::vector<ftm::idNode>> allBranchOrigins(
      tree->getNumberOfNodes());
    std::vector<int> allBranchOriginsSize(tree->getNumberOfNodes());
    std::queue<ftm::idNode> queueCrossing;

    // ----- Get important and non-important pairs gap and store saddle nodes of
    // each branch
    int maxSize = std::numeric_limits<int>::lowest();
    for(auto leaf : leaves)
      queueCrossing.emplace(leaf);
    while(!queueCrossing.empty()) {
      ftm::idNode node = queueCrossing.front();
      queueCrossing.pop();
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();

      // Get saddle nodes in the branch
      auto tupBranchOrigins = tree->getBranchOriginsFromThisBranch(nodeOrigin);
      allBranchOrigins[nodeOrigin] = std::get<0>(tupBranchOrigins);
      std::vector<ftm::idNode> nonBranchOrigins = std::get<1>(tupBranchOrigins);
      allBranchOrigins[nodeOrigin].insert(allBranchOrigins[nodeOrigin].end(),
                                          nonBranchOrigins.begin(),
                                          nonBranchOrigins.end());
      std::sort(allBranchOrigins[nodeOrigin].begin(),
                allBranchOrigins[nodeOrigin].end(), compValue);
      allBranchOriginsSize[nodeOrigin] = allBranchOrigins[nodeOrigin].size();

      // Get sizes of sub-branches if they are non-important
      for(size_t i = 0; i < allBranchOrigins[nodeOrigin].size(); ++i) {
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        bool isSubBranchImportant
          = tree->isImportantPair<dataType>(branchNodeOrigin, ImportantPairs);
        if(not isSubBranchImportant)
          allBranchOriginsSize[nodeOrigin]
            += allBranchOriginsSize[branchNodeOrigin];
      }

      if(tree->isImportantPair<dataType>(nodeOrigin, ImportantPairs))
        maxSize = std::max(maxSize, allBranchOriginsSize[nodeOrigin]);
    }
    double nonImportantPairsGap
      = (rootYmax - rootYmin) * 0.005 * NonImportantPairsSpacing;
    double importantPairsGap = (maxSize)*nonImportantPairsGap * 1.05;
    bool customImportantPairsSpacing
      = importantPairsGap < ImportantPairsSpacing;
    if(customImportantPairsSpacing)
      importantPairsGap = ImportantPairsSpacing;

    // ----- Positioning of branches and avoid conflict
    for(auto leaf : leaves)
      queueCrossing.emplace(leaf);
    while(!queueCrossing.empty()) {
      ftm::idNode node = queueCrossing.front();
      queueCrossing.pop();
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();

      // Prepositioning of branches
      // auto restrictedBounds = getBranchBounds(retVec, treeSimplexId,
      // branching, tree, nodeOrigin, true);
      auto restrictedBounds = std::make_tuple(
        retVec[treeSimplexId[node] * 2], retVec[treeSimplexId[node] * 2], 0, 0);
      for(size_t i = 0; i < allBranchOrigins[nodeOrigin].size(); ++i) {
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        ftm::idNode branchNode = tree->getNode(branchNodeOrigin)->getOrigin();

        bool isSubBranchImportant
          = tree->isImportantPair<dataType>(branchNodeOrigin, ImportantPairs);
        bool toLeft = not isSubBranchImportant;

        // float branchNodeOriginXmin =
        // std::get<0>(allBranchBounds[branchNodeOrigin]);
        float branchNodeOriginXmax
          = std::get<1>(allBranchBounds[branchNodeOrigin]);
        float shift
          = toLeft ? std::get<0>(restrictedBounds) - branchNodeOriginXmax :
                   // std::get<1>(restrictedBounds) - branchNodeOriginXmin;
              std::get<1>(restrictedBounds)
                - retVec[treeSimplexId[branchNode] * 2];
        shift += (toLeft ? -1 : 1)
                 * (isSubBranchImportant ? importantPairsGap
                                         : nonImportantPairsGap);
        // shift += (toLeft ? -1 : 1) * nonImportantPairsGap;
        shiftBranchBounds(retVec, treeSimplexId, branching, allBranchBounds,
                          allBranchOrigins[branchNodeOrigin], tree,
                          branchNodeOrigin, shift);
      }

      // Shift a branch if conflict with another one
      for(size_t i = 1; i < allBranchOrigins[nodeOrigin].size(); ++i) {
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        ftm::idNode branchNode = tree->getNode(branchNodeOrigin)->getOrigin();
        for(size_t j = 0; j < i; ++j) {
          auto first = allBranchBounds[branchNodeOrigin];
          ftm::idNode previousBranchNodeOrigin
            = allBranchOrigins[nodeOrigin][j];
          auto second = allBranchBounds[previousBranchNodeOrigin];

          bool branchConflict = isConflictingBranchAndBound(
            first, second, tree, previousBranchNodeOrigin, retVec,
            treeSimplexId);
          if(isConflictingBounds(first, second) or branchConflict) {

            // Get left or right orientation given the branch
            int lastIndex = nodeBranching[branchNodeOrigin].size() - 1;
            bool isLeft
              = (retVec
                   [treeSimplexId[nodeBranching[branchNodeOrigin][lastIndex]]
                    * 2]
                 //< retVec[treeSimplexId[nodeOrigin]*2]);
                 < retVec[treeSimplexId[node] * 2]);

            // Get shift
            float branchNodeOriginXmax = std::get<1>(first);
            float previousbranchNodeOriginXmin = std::get<0>(second);
            // float branchNodeOriginXmin = std::get<0>(first);
            float previousbranchNodeOriginXmax = std::get<1>(second);
            float shift
              = isLeft ? previousbranchNodeOriginXmin - branchNodeOriginXmax :
                       // previousbranchNodeOriginXmax - branchNodeOriginXmin;
                  previousbranchNodeOriginXmax
                    - retVec[treeSimplexId[branchNode] * 2];
            bool isSubBranchImportant = tree->isImportantPair<dataType>(
              branchNodeOrigin, ImportantPairs);
            shift += (isLeft ? -1 : 1)
                     * (isSubBranchImportant ? importantPairsGap
                                             : nonImportantPairsGap);
            // shift += (isLeft ? -1 : 1) * nonImportantPairsGap;

            // Shift bounds
            shiftBranchBounds(retVec, treeSimplexId, branching, allBranchBounds,
                              allBranchOrigins[branchNodeOrigin], tree,
                              branchNodeOrigin, shift);
          }
        }
      } // end for

      // TODO optimize get branch bounds by using results previously computed
      // Get branch x and y bounds
      allBranchBounds[nodeOrigin]
        = getBranchBounds(retVec, treeSimplexId, branching, tree, nodeOrigin);
    } // end while

    // ----- Correction of important/non-important pairs gap
    // TODO the gap between important pairs can be higher than the minimum gap
    // needed to avoid conflict. The gap is computed using the maximum number of
    // non-important pairs attached to an inmportant pairs Unfortunately the
    // real gap can only be computed here, after the conflicts has been avoided.
    // The maximum real gap must be calculated and propagated to all important
    // branches and we also need to manage to avoid conflict with this new gap.
    // Get real gap
    double realImportantPairsGap = std::numeric_limits<double>::lowest();
    /*if(customImportantPairsSpacing)
      realImportantPairsGap = importantPairsGap;
    else{
      for(auto leaf : leaves)
        queueCrossing.emplace(leaf);
      while(!queueCrossing.empty()){
        ftm::idNode node = queueCrossing.front();
        queueCrossing.pop();
        ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();
        //if(isRoot(tree, nodeOrigin)) continue; // TODO manage root gap and gap
    for the others

        bool isBranchImportant = isImportantPair<dataType>(tree, nodeOrigin,
    ImportantPairs); if(not isBranchImportant) continue;

        double gap = retVec[treeSimplexId[node]*2] -
    std::get<0>(allBranchBounds[nodeOrigin]); realImportantPairsGap =
    std::max(gap, realImportantPairsGap);
      }
      realImportantPairsGap *= 1.05;
    }*/
    realImportantPairsGap = importantPairsGap;

    // Shift given real gap
    for(auto leaf : leaves)
      queueCrossing.emplace(leaf);
    while(!queueCrossing.empty()) {
      ftm::idNode node = queueCrossing.front();
      queueCrossing.pop();
      ftm::idNode nodeOrigin = tree->getNode(node)->getOrigin();

      bool isBranchImportant
        = tree->isImportantPair<dataType>(nodeOrigin, ImportantPairs);
      if(not isBranchImportant)
        continue;

      for(size_t i = 0; i < allBranchOrigins[nodeOrigin].size(); ++i) {
        ftm::idNode branchNodeOrigin = allBranchOrigins[nodeOrigin][i];
        bool isSubBranchImportant
          = tree->isImportantPair<dataType>(branchNodeOrigin, ImportantPairs);
        double shift = 0;
        if(not isSubBranchImportant) {
          double gap = retVec[treeSimplexId[node] * 2]
                       - std::get<0>(allBranchBounds[nodeOrigin]);
          shift = -(realImportantPairsGap - gap) * NonImportantPairsProximity;
        } else {
          shift = -(
            importantPairsGap
            - realImportantPairsGap); // TODO multi shift depending on conflict
        }
        shiftBranchBounds(retVec, treeSimplexId, branching, allBranchBounds,
                          allBranchOrigins[branchNodeOrigin], tree,
                          branchNodeOrigin, shift);
      }
    }

    if(verbose > 0) {
      std::stringstream ss;
      ss << "AVOID CROSSING  = " << t_avoid.getElapsedTime();
      printMsg(ss.str(), debug::Priority::VERBOSE);
      printMsg(debug::Separator::L2, debug::Priority::VERBOSE);
    }

    return retVec;
  }

  template <class dataType>
  std::vector<float> treePlanarLayout(
    FTMTree_MT *tree,
    std::tuple<double, double, double, double, double, double> oldBounds,
    double refPersistence) {
    std::vector<float> res;
    res = treePlanarLayoutImpl<dataType>(tree, oldBounds, refPersistence);
    return res;
  }

  // ==========================================================================
  // Bounds Utils
  // ==========================================================================
  void printTuple(std::tuple<float, float, float, float> tup) {
    printMsg(debug::Separator::L2, debug::Priority::VERBOSE);
    std::stringstream ss;
    ss << std::get<0>(tup) << " _ " << std::get<1>(tup) << " _ "
       << std::get<2>(tup) << " _ " << std::get<3>(tup) << " _ ";
    printMsg(ss.str(), debug::Priority::VERBOSE);
  }

  // Branchroot must be a non-leaf node
  std::tuple<float, float, float, float>
    getBranchBounds(std::vector<float> &retVec,
                    std::vector<LongSimplexId> &treeSimplexId,
                    std::vector<ftm::idNode> &branching,
                    FTMTree_MT *tree,
                    ftm::idNode branchRoot,
                    bool restricted = false) {
    float x_min = std::numeric_limits<float>::max();
    float y_min = std::numeric_limits<float>::max();
    float x_max = std::numeric_limits<float>::lowest();
    float y_max = std::numeric_limits<float>::lowest();

    std::queue<ftm::idNode> queue;
    queue.emplace(branchRoot);
    while(!queue.empty()) {
      ftm::idNode node = queue.front();
      queue.pop();

      // Skip if we go in the branch in which is branchRoot
      if(branching[node] != branchRoot
         and tree->getParentSafe(node) == branchRoot and node != branchRoot)
        continue;

      // Skip if restricted
      if(restricted and !tree->isLeaf(node) and branching[node] != branchRoot
         and node != branchRoot)
        continue;

      y_min = std::min(y_min, retVec[treeSimplexId[node] * 2 + 1]);
      y_max = std::max(y_max, retVec[treeSimplexId[node] * 2 + 1]);
      if(node != branchRoot) {
        x_min = std::min(x_min, retVec[treeSimplexId[node] * 2]);
        x_max = std::max(x_max, retVec[treeSimplexId[node] * 2]);
      }

      for(auto child : tree->getChildren(node))
        queue.emplace(child);
    }

    return std::make_tuple(x_min, x_max, y_min, y_max);
  }

  bool
    isConflictingBoundsXOneWay(std::tuple<float, float, float, float> first,
                               std::tuple<float, float, float, float> second) {
    double eps = std::numeric_limits<float>::epsilon();
    return (std::get<0>(first) <= std::get<0>(second)
            and std::get<0>(second) <= std::get<1>(first))
           or (std::get<0>(first) <= std::get<1>(second)
               and std::get<1>(second) <= std::get<1>(first))
           or (isEqual<float>(std::get<0>(first), std::get<0>(second), eps)
               or isEqual<float>(std::get<0>(second), std::get<1>(first), eps)
               or isEqual<float>(std::get<0>(first), std::get<1>(second), eps)
               or isEqual<float>(std::get<1>(second), std::get<1>(first), eps));
  }

  bool isConflictingBoundsX(std::tuple<float, float, float, float> first,
                            std::tuple<float, float, float, float> second) {
    return isConflictingBoundsXOneWay(first, second)
           or isConflictingBoundsXOneWay(second, first);
  }

  bool
    isConflictingBoundsYOneWay(std::tuple<float, float, float, float> first,
                               std::tuple<float, float, float, float> second) {
    double eps = std::numeric_limits<float>::epsilon();
    return (std::get<2>(first) <= std::get<2>(second)
            and std::get<2>(second) <= std::get<3>(first))
           or (std::get<2>(first) <= std::get<3>(second)
               and std::get<3>(second) <= std::get<3>(first))
           or (isEqual<float>(std::get<2>(first), std::get<2>(second), eps)
               or isEqual<float>(std::get<2>(second), std::get<3>(first), eps)
               or isEqual<float>(std::get<2>(first), std::get<3>(second), eps)
               or isEqual<float>(std::get<3>(second), std::get<3>(first), eps));
  }

  bool isConflictingBoundsY(std::tuple<float, float, float, float> first,
                            std::tuple<float, float, float, float> second) {
    return isConflictingBoundsYOneWay(first, second)
           or isConflictingBoundsYOneWay(second, first);
  }

  bool isConflictingBounds(std::tuple<float, float, float, float> first,
                           std::tuple<float, float, float, float> second) {
    return isConflictingBoundsX(first, second)
           and isConflictingBoundsY(first, second);
  }

  bool
    isConflictingBranchAndBound(std::tuple<float, float, float, float> first,
                                std::tuple<float, float, float, float> second,
                                FTMTree_MT *tree,
                                ftm::idNode branchNodeOrigin,
                                std::vector<float> &retVec,
                                std::vector<LongSimplexId> &treeSimplexId) {
    float xBranchNodeOrigin = retVec[treeSimplexId[branchNodeOrigin] * 2];
    float xBranchNode
      = retVec[treeSimplexId[tree->getNode(branchNodeOrigin)->getOrigin()] * 2];
    float myMin = std::min(xBranchNode, xBranchNodeOrigin);
    float myMax = std::max(xBranchNode, xBranchNodeOrigin);
    auto branchBounds = std::make_tuple(myMin, myMax, 0, 0);
    return isConflictingBoundsX(first, branchBounds)
           and isConflictingBoundsY(first, second);
  }

  std::tuple<float, float, float, float>
    shiftBranchBoundsTuple(std::tuple<float, float, float, float> branchBound,
                           float realShift) {
    return std::make_tuple(std::get<0>(branchBound) + realShift,
                           std::get<1>(branchBound) + realShift,
                           std::get<2>(branchBound), std::get<3>(branchBound));
  }

  void shiftBranchBounds(
    std::vector<float> &retVec,
    std::vector<LongSimplexId> &treeSimplexId,
    std::vector<ftm::idNode> &branching,
    std::vector<std::tuple<float, float, float, float>> &allBranchBounds,
    std::vector<ftm::idNode> &branchOrigins,
    FTMTree_MT *tree,
    ftm::idNode branchRoot,
    float shift) {
    std::queue<ftm::idNode> queue;
    queue.emplace(branchRoot);
    while(!queue.empty()) {
      ftm::idNode node = queue.front();
      queue.pop();

      if(branching[node] != branchRoot
         and tree->getParentSafe(node) == branchRoot and node != branchRoot)
        continue;

      if(node != branchRoot)
        retVec[treeSimplexId[node] * 2] += shift;

      for(auto child : tree->getChildren(node))
        queue.emplace(child);
    }
    allBranchBounds[branchRoot]
      = shiftBranchBoundsTuple(allBranchBounds[branchRoot], shift);
    for(auto node : branchOrigins)
      allBranchBounds[node]
        = shiftBranchBoundsTuple(allBranchBounds[node], shift);
  }

  std::vector<double> tupleToVector(
    std::tuple<double, double, double, double, double, double> &tup) {
    std::vector<double> vec{std::get<0>(tup), std::get<1>(tup),
                            std::get<2>(tup), std::get<3>(tup),
                            std::get<4>(tup), std::get<5>(tup)};
    return vec;
  }

  std::tuple<double, double, double, double, double, double>
    vectorToTuple(std::vector<double> &vec) {
    return std::make_tuple(vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]);
  }

  std::tuple<double, double, double, double, double, double> getMaximalBounds(
    std::vector<std::tuple<double, double, double, double, double, double>>
      &allBounds,
    std::vector<int> &clusteringAssignmentT,
    int clusterID) {
    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    double y_max = std::numeric_limits<double>::lowest();
    double z_max = std::numeric_limits<double>::lowest();
    for(size_t i = 0; i < allBounds.size(); ++i)
      if(clusteringAssignmentT[i] == clusterID) {
        x_min = std::min(x_min, std::get<0>(allBounds[i]));
        x_max = std::max(x_max, std::get<1>(allBounds[i]));
        y_min = std::min(y_min, std::get<2>(allBounds[i]));
        y_max = std::max(y_max, std::get<3>(allBounds[i]));
        z_min = std::min(z_min, std::get<4>(allBounds[i]));
        z_max = std::max(z_max, std::get<5>(allBounds[i]));
      }
    return std::make_tuple(x_min, x_max, y_min, y_max, z_min, z_max);
  }
};

#endif

#include <FTMTree.h>

namespace ttk {
  namespace ftm {

    // --------------------
    // Is
    // --------------------

    bool FTMTree_MT::isNodeOriginDefined(idNode nodeId) {
      unsigned int origin = (unsigned int)this->getNode(nodeId)->getOrigin();
      return origin != nullNodes and origin < this->getNumberOfNodes()
             and origin >= 0;
      ;
    }

    bool FTMTree_MT::isRoot(idNode nodeId) {
      return this->getNode(nodeId)->getNumberOfUpSuperArcs() == 0;
    }

    bool FTMTree_MT::isLeaf(idNode nodeId) {
      return this->getNode(nodeId)->getNumberOfDownSuperArcs() == 0;
    }

    bool FTMTree_MT::isNodeAlone(idNode nodeId) {
      return this->isRoot(nodeId) and this->isLeaf(nodeId);
    }

    bool FTMTree_MT::isFullMerge() {
      idNode treeRoot = this->getRoot();
      return (unsigned int)this->getNode(treeRoot)->getOrigin() == treeRoot;
    }

    bool FTMTree_MT::isBranchOrigin(idNode nodeId) {
      return this->getParent(this->getNode(nodeId)->getOrigin()) != nodeId;
    }

    // --------------------
    // Get
    // --------------------

    idNode FTMTree_MT::getRoot() {
      for(idNode node = 0; node < this->getNumberOfNodes(); ++node)
        if(this->isRoot(node) and !this->isLeaf(node))
          return node;
      return nullNodes;
    }

    idNode FTMTree_MT::getParentSafe(idNode nodeId) {
      if(!this->isRoot(nodeId)) {
        // _ Nodes in merge trees should have only one parent
        idSuperArc arcId = this->getNode(nodeId)->getUpSuperArcId(0);
        idNode parentNodeId = this->getSuperArc(arcId)->getUpNodeId();
        return parentNodeId;
      }
      return nodeId;
    }

    std::vector<idNode> FTMTree_MT::getChildren(idNode nodeId) {
      std::vector<idNode> childrens;
      for(idSuperArc i = 0;
          i < this->getNode(nodeId)->getNumberOfDownSuperArcs(); ++i) {
        idSuperArc arcId = this->getNode(nodeId)->getDownSuperArcId(i);
        childrens.push_back(this->getSuperArc(arcId)->getDownNodeId());
      }
      return childrens;
    }

    std::vector<idNode> FTMTree_MT::getLeavesFromTree() {
      std::vector<idNode> treeLeaves;
      for(idNode i = 0; i < this->getNumberOfNodes(); ++i) {
        if(this->isLeaf(i) and !this->isRoot(i))
          treeLeaves.push_back(i);
      }
      return treeLeaves;
    }

    int FTMTree_MT::getNumberOfLeavesFromTree() {
      auto leaves = this->getLeaves();
      return leaves.size();
    }

    int FTMTree_MT::getNumberOfNodeAlone() {
      int cpt = 0;
      for(idNode i = 0; i < this->getNumberOfNodes(); ++i)
        cpt += this->isNodeAlone(i) ? 1 : 0;
      return cpt;
    }

    int FTMTree_MT::getRealNumberOfNodes() {
      return this->getNumberOfNodes() - this->getNumberOfNodeAlone();
    }

    std::tuple<std::vector<idNode>, std::vector<idNode>>
      FTMTree_MT::getBranchOriginsFromThisBranch(idNode node) {
      std::vector<idNode> branchOrigins, nonBranchOrigins;

      idNode nodeOrigin = this->getNode(node)->getOrigin();
      idNode nodeParent = this->getParent(nodeOrigin);
      while(nodeParent != node) {
        if(this->isBranchOrigin(nodeParent))
          branchOrigins.push_back(nodeParent);
        else
          nonBranchOrigins.push_back(nodeParent);
        nodeParent = this->getParent(nodeParent);
      }

      return std::make_tuple(branchOrigins, nonBranchOrigins);
    }

    void FTMTree_MT::getTreeBranching(
      std::vector<idNode> &branching,
      std::vector<int> &branchingID,
      std::vector<std::vector<idNode>> &nodeBranching) {
      branching = std::vector<idNode>(this->getNumberOfNodes());
      branchingID = std::vector<int>(this->getNumberOfNodes(), -1);
      nodeBranching
        = std::vector<std::vector<idNode>>(this->getNumberOfNodes());
      int branchID = 0;
      std::queue<idNode> queue;
      queue.emplace(this->getRoot());
      while(!queue.empty()) {
        idNode node = queue.front();
        queue.pop();
        if(this->isLeaf(node))
          continue;
        auto nodeOrigin = this->getNode(node)->getOrigin();
        idNode parentNodeOrigin = nodeOrigin;
        while(parentNodeOrigin != node) {
          branching[parentNodeOrigin] = node;
          branchingID[parentNodeOrigin] = branchID;
          nodeBranching[node].push_back(parentNodeOrigin);
          parentNodeOrigin = this->getParent(parentNodeOrigin);
        }
        if(this->isRoot(node)) {
          branching[node] = node;
          branchingID[node] = branchID;
        }
        ++branchID;
        for(auto child : this->getChildren(node))
          queue.emplace(child);
      }
    }

    void FTMTree_MT::getTreeBranching(std::vector<idNode> &branching,
                                      std::vector<int> &branchingID) {
      std::vector<std::vector<idNode>> nodeBranching;
      this->getTreeBranching(branching, branchingID, nodeBranching);
    }

  } // namespace ftm
} // namespace ttk

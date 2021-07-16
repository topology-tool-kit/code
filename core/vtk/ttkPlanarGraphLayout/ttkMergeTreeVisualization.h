/// \ingroup vtk
/// \class ttkMergeTreeVisualization
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Visualization module for merge trees.
///

#pragma once

#include <MergeTreeVisualization.h>

#include <FTMTree.h>

#include <Debug.h>
#include <ttkAlgorithm.h>

// VTK Includes
#include <vtkAppendFilter.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

#include <stack>

using namespace ttk;
using namespace ftm;

class ttkMergeTreeVisualization : public MergeTreeVisualization {
private:
  // Visualization parameters
  double DimensionSpacing = 1.;
  int DimensionToShift = 0;
  bool OutputSegmentation = false;
  int ShiftMode = 0; // 0: Star ; 1: Star Barycenter ; 2: Line ; 3: Double Line

  // Offset
  int iSampleOffset = 0;
  int noSampleOffset = 0;

  // Used for critical type and for point coordinates
  std::vector<vtkUnstructuredGrid *> treesNodes;
  std::vector<std::vector<int>>
    treesNodeCorrMesh; // used to acces treesNodes given input trees

  // Segmentation
  std::vector<vtkDataSet *> treesSegmentation;

  // Clustering output
  std::vector<int> clusteringAssignment;
  std::vector<std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>>
    outputMatchingBarycenter;

  // Barycenter output
  std::vector<std::vector<float>> allBaryPercentMatch;

  // Temporal Subsampling Output
  std::vector<bool> interpolatedTrees;

  // Output
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode;
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputArc;
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputSegmentation;

  // Matching output
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputNode1,
    vtkOutputNode2; // input data
  std::vector<std::vector<SimplexId>> nodeCorr1, nodeCorr2;
  vtkSmartPointer<vtkUnstructuredGrid> vtkOutputMatching; // output

  // Filled by the algorithm
  std::vector<std::vector<SimplexId>> nodeCorr;
  std::vector<double> clusterShift;

public:
  ttkMergeTreeVisualization(){};
  ~ttkMergeTreeVisualization(){};

  // ==========================================================================
  // Getter / Setter
  // ==========================================================================
  // Visualization parameters
  void setDimensionSpacing(double d) {
    DimensionSpacing = d;
  }
  void setDimensionToShift(int i) {
    DimensionToShift = i;
  }
  void setOutputSegmentation(bool b) {
    OutputSegmentation = b;
  }
  void setShiftMode(int i) {
    ShiftMode = i;
  }

  // Offset
  void setISampleOffset(int i) {
    iSampleOffset = i;
  }
  void setNoSampleOffset(int i) {
    noSampleOffset = i;
  }

  // Used for critical type and if not planar layout for point coordinates
  void setTreesNodes(std::vector<vtkUnstructuredGrid *> &nodes) {
    treesNodes = nodes;
  }
  void setTreesNodeCorrMesh(std::vector<std::vector<int>> &nodeCorrMesh) {
    treesNodeCorrMesh = nodeCorrMesh;
  }
  void setTreesNodes(vtkUnstructuredGrid *nodes) {
    treesNodes.clear();
    treesNodes.emplace_back(nodes);
  }
  void setTreesNodeCorrMesh(std::vector<int> &nodeCorrMesh) {
    treesNodeCorrMesh.clear();
    treesNodeCorrMesh.emplace_back(nodeCorrMesh);
  }

  // Segmentation
  void setTreesSegmentation(std::vector<vtkDataSet *> &segmentation) {
    treesSegmentation = segmentation;
  }
  void setTreesSegmentation(vtkDataSet *segmentation) {
    treesSegmentation.clear();
    treesSegmentation.emplace_back(segmentation);
  }

  // Clustering output
  void setClusteringAssignment(std::vector<int> &asgn) {
    clusteringAssignment = asgn;
  }
  void setOutputMatchingBarycenter(
    std::vector<std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>>
      &matching) {
    outputMatchingBarycenter = matching;
  }

  // Barycenter output
  std::vector<std::vector<float>> getAllBaryPercentMatch() {
    return allBaryPercentMatch;
  }
  void
    setAllBaryPercentMatch(std::vector<std::vector<float>> &baryPercentMatch) {
    allBaryPercentMatch = baryPercentMatch;
  }

  // Temporal Subsampling Output
  void setInterpolatedTrees(std::vector<bool> &isInterpolatedTrees) {
    interpolatedTrees = isInterpolatedTrees;
  }

  // Output
  void setVtkOutputNode(vtkSmartPointer<vtkUnstructuredGrid> vtkNode) {
    vtkOutputNode = vtkNode;
  }
  void setVtkOutputArc(vtkSmartPointer<vtkUnstructuredGrid> vtkArc) {
    vtkOutputArc = vtkArc;
  }
  void setVtkOutputSegmentation(
    vtkSmartPointer<vtkUnstructuredGrid> vtkSegmentation) {
    vtkOutputSegmentation = vtkSegmentation;
  }

  // Matching output
  void setVtkOutputNode1(vtkSmartPointer<vtkUnstructuredGrid> vtkNode1) {
    vtkOutputNode1 = vtkNode1;
  }
  void setVtkOutputNode2(vtkSmartPointer<vtkUnstructuredGrid> vtkNode2) {
    vtkOutputNode2 = vtkNode2;
  }
  void setNodeCorr1(std::vector<std::vector<SimplexId>> &nodeCorrT) {
    nodeCorr1 = nodeCorrT;
  }
  void setNodeCorr2(std::vector<std::vector<SimplexId>> &nodeCorrT) {
    nodeCorr2 = nodeCorrT;
  }
  void setVtkOutputMatching(vtkSmartPointer<vtkUnstructuredGrid> vtkMatching) {
    vtkOutputMatching = vtkMatching;
  }
  void setOutputMatching(
    std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> &matching) {
    outputMatchingBarycenter = std::vector<
      std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>>(1);
    outputMatchingBarycenter[0]
      = std::vector<std::vector<std::tuple<ftm::idNode, ftm::idNode>>>(1);

    std::vector<std::tuple<ftm::idNode, ftm::idNode>> tMatching;
    for(auto match : matching)
      tMatching.emplace_back(
        std::make_tuple(std::get<0>(match), std::get<1>(match)));

    outputMatchingBarycenter[0][0] = tMatching;
  }

  // Filled by the algorithm
  std::vector<std::vector<SimplexId>> getNodeCorr() {
    return nodeCorr;
  }

  std::vector<double> getClusterShift() {
    return clusterShift;
  }

  // ==========================================================================
  // Matching Visualization
  // ==========================================================================
  template <class dataType>
  void
    makeMatchingOutput(FTMTree_MT *tree1, FTMTree_MT *tree2, int verbose = 0) {
    std::vector<FTMTree_MT *> trees{tree1, tree2};
    std::vector<FTMTree_MT *> barycenters;

    makeMatchingOutput<dataType>(trees, barycenters, verbose);
  }

  template <class dataType>
  void makeMatchingOutput(std::vector<FTMTree_MT *> &trees,
                          std::vector<FTMTree_MT *> &barycenters,
                          int verbose = 0) {
    // verbose = 2;

    int numInputs = trees.size();
    int NumberOfBarycenters = barycenters.size();
    bool clusteringOutput = (NumberOfBarycenters != 0);
    NumberOfBarycenters
      = std::max(NumberOfBarycenters, 1); // to always enter the outer loop
    if(not clusteringOutput)
      numInputs = 1;

    vtkSmartPointer<vtkUnstructuredGrid> vtkMatching
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> pointsM = vtkSmartPointer<vtkPoints>::New();

    // Fields
    vtkNew<vtkIntArray> matchingID{};
    matchingID->SetName("MatchingID");
    vtkNew<vtkIntArray> matchingType{};
    matchingType->SetName("MatchingType");
    vtkNew<vtkDoubleArray> matchPers{};
    matchPers->SetName("MeanMatchedPersistence");
    // vtkNew<vtkDoubleArray> costArray{}; // TODO
    // costArray->SetName("Cost");
    vtkNew<vtkIntArray> tree1NodeIdField{};
    tree1NodeIdField->SetName("tree1NodeId");
    vtkNew<vtkIntArray> tree2NodeIdField{};
    tree2NodeIdField->SetName("tree2NodeId");

    vtkNew<vtkFloatArray> matchingPercentMatch{};
    matchingPercentMatch->SetName("MatchingPercentMatch");

    // Iterate through clusters and trees
    if(verbose > 0)
      printMsg(
        "// Iterate through clusters and trees", debug::Priority::VERBOSE);
    int count = 0;
    for(int c = 0; c < NumberOfBarycenters; ++c) {
      for(int i = 0; i < numInputs; ++i) {
        for(std::tuple<idNode, idNode> match : outputMatchingBarycenter[c][i]) {
          vtkIdType pointIds[2];
          ftm::idNode tree1NodeId = std::get<0>(match);
          ftm::idNode tree2NodeId = std::get<1>(match);
          // double cost = std::get<2>(match);
          FTMTree_MT *tree1;
          FTMTree_MT *tree2 = trees[i];
          if(not clusteringOutput) {
            tree1 = trees[0];
            tree2 = trees[1];
          } else
            tree1 = barycenters[c];

          // Get first point
          if(verbose > 1)
            printMsg("// Get first point", debug::Priority::VERBOSE);
          SimplexId pointToGet1 = clusteringOutput ? nodeCorr2[c][tree1NodeId]
                                                   : nodeCorr1[0][tree1NodeId];
          double *point1 = vtkOutputNode2->GetPoints()->GetPoint(pointToGet1);
          const SimplexId nextPointId1 = pointsM->InsertNextPoint(point1);
          pointIds[0] = nextPointId1;

          // Get second point
          if(verbose > 1)
            printMsg("// Get second point", debug::Priority::VERBOSE);
          SimplexId pointToGet2 = clusteringOutput ? nodeCorr1[i][tree2NodeId]
                                                   : nodeCorr1[1][tree2NodeId];
          double *point2 = vtkOutputNode1->GetPoints()->GetPoint(pointToGet2);
          const SimplexId nextPointId2 = pointsM->InsertNextPoint(point2);
          pointIds[1] = nextPointId2;

          // Add cell
          if(verbose > 1)
            printMsg("// Add cell", debug::Priority::VERBOSE);
          vtkMatching->InsertNextCell(VTK_LINE, 2, pointIds);

          // Add arc matching percentage
          if(verbose > 1)
            printMsg(
              "// Add arc matching percentage", debug::Priority::VERBOSE);
          if(allBaryPercentMatch.size() != 0)
            matchingPercentMatch->InsertNextTuple1(
              allBaryPercentMatch[c][tree1NodeId]);

          // Add tree1 and tree2 node ids
          if(verbose > 1)
            printMsg(
              "// Add tree1 and tree2 node ids", debug::Priority::VERBOSE);
          tree1NodeIdField->InsertNextTuple1(pointToGet1);
          tree2NodeIdField->InsertNextTuple1(pointToGet2);

          // Add matching ID
          matchingID->InsertNextTuple1(count);

          // Add matching type
          if(verbose > 1)
            printMsg("// Add matching type", debug::Priority::VERBOSE);
          int thisType = 0;
          int tree1NodeDown
            = tree1->getNode(tree1NodeId)->getNumberOfDownSuperArcs();
          int tree1NodeUp
            = tree1->getNode(tree1NodeId)->getNumberOfUpSuperArcs();
          int tree2NodeDown
            = tree2->getNode(tree2NodeId)->getNumberOfDownSuperArcs();
          int tree2NodeUp
            = tree2->getNode(tree2NodeId)->getNumberOfUpSuperArcs();
          if(tree1NodeDown != 0 and tree1NodeUp != 0 and tree2NodeDown != 0
             and tree2NodeUp != 0)
            thisType = 1; // Saddle to Saddle
          if(tree1NodeDown == 0 and tree1NodeUp != 0 and tree2NodeDown == 0
             and tree2NodeUp != 0)
            thisType = 2; // Leaf to leaf
          if(tree1NodeDown != 0 and tree1NodeUp == 0 and tree2NodeDown != 0
             and tree2NodeUp == 0)
            thisType = 3; // Root to root
          matchingType->InsertNextTuple1(thisType);

          // Add mean matched persistence
          if(verbose > 1)
            printMsg(
              "// Add mean matched persistence", debug::Priority::VERBOSE);
          double tree1Pers = tree1->getNodePersistence<dataType>(tree1NodeId);
          double tree2Pers = tree2->getNodePersistence<dataType>(tree2NodeId);
          double meanPersistence = (tree1Pers + tree2Pers) / 2;
          matchPers->InsertNextTuple1(meanPersistence);

          // Add cost
          // costArray->InsertNextTuple1(cost);

          count++;
        }
      }
    }
    vtkMatching->SetPoints(pointsM);
    vtkMatching->GetCellData()->AddArray(matchingType);
    vtkMatching->GetCellData()->AddArray(matchPers);
    vtkMatching->GetCellData()->AddArray(matchingID);
    // vtkMatching->GetCellData()->AddArray(costArray);
    vtkMatching->GetCellData()->AddArray(tree1NodeIdField);
    vtkMatching->GetCellData()->AddArray(tree2NodeIdField);
    if(allBaryPercentMatch.size() != 0)
      vtkMatching->GetCellData()->AddArray(matchingPercentMatch);
    vtkOutputMatching->ShallowCopy(vtkMatching);
  }

  // ==========================================================================
  // Trees Visualization
  // ==========================================================================
  template <class dataType>
  void makeTreesOutput(FTMTree_MT *tree1, int verbose = 0) {
    std::vector<FTMTree_MT *> trees{tree1};

    makeTreesOutput<dataType>(trees, verbose);
  }

  template <class dataType>
  void makeTreesOutput(FTMTree_MT *tree1, FTMTree_MT *tree2, int verbose = 0) {
    std::vector<FTMTree_MT *> trees{tree1, tree2};

    makeTreesOutput<dataType>(trees, verbose);
  }

  template <class dataType>
  void makeTreesOutput(std::vector<FTMTree_MT *> &trees, int verbose = 0) {
    std::vector<FTMTree_MT *> barycenters;
    clusteringAssignment = std::vector<int>(trees.size(), 0);

    makeTreesOutput<dataType>(trees, barycenters, verbose);
  }

  template <class dataType>
  void makeTreesOutput(std::vector<FTMTree_MT *> &trees,
                       std::vector<FTMTree_MT *> &barycenters,
                       int verbose = 0) {
    // verbose = 5;

    int numInputs = trees.size();
    int numInputsOri = numInputs;
    int NumberOfBarycenters = barycenters.size();
    bool clusteringOutput = (NumberOfBarycenters != 0);
    NumberOfBarycenters
      = std::max(NumberOfBarycenters, 1); // to always enter the outer loop

    // Bounds
    std::vector<std::tuple<double, double, double, double, double, double>>
      allBounds(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      if(OutputSegmentation) {
        double *tBounds = treesSegmentation[i]->GetBounds();
        allBounds[i] = std::make_tuple(tBounds[0], tBounds[1], tBounds[2],
                                       tBounds[3], tBounds[4], tBounds[5]);
      } else
        allBounds[i]
          = getRealBounds(treesNodes[i], trees[i], treesNodeCorrMesh[i]);
      /*if(PlanarLayout){
        // TODO correctly manage bounds when planar layout
        std::vector<double> tBound = tupleToVector(allBounds[i]);
        int level = getTreeDepth(trees[i]);
        for(int j = 0; j < tBound.size(); ++j)
          tBound[j] = tBound[j]*level/5;
        allBounds[i] = vectorToTuple(tBound);
      }*/
    }
    std::vector<std::tuple<double, double, double, double, double, double>>
      allBaryBounds(barycenters.size());
    std::vector<std::vector<ftm::idNode>> allBaryBranching(barycenters.size());
    std::vector<std::vector<int>> allBaryBranchingID(barycenters.size());
    for(size_t c = 0; c < barycenters.size(); ++c) {
      allBaryBounds[c] = getMaximalBounds(allBounds, clusteringAssignment, c);
      barycenters[c]->getTreeBranching(
        allBaryBranching[c], allBaryBranchingID[c]);
    }
    if(not clusteringOutput)
      allBaryBounds.emplace_back(
        getMaximalBounds(allBounds, clusteringAssignment, 0));

    // ----------------------------------------------------------------------
    // Make Trees Output
    // ----------------------------------------------------------------------
    if(verbose > 0)
      printMsg("--- Make Trees Output", debug::Priority::VERBOSE);
    std::vector<FTMTree_MT *> treesOri(trees);
    if(ShiftMode == 1) { // Star Barycenter
      trees.clear();
      clusteringAssignment.clear();
      for(unsigned int j = 0; j < barycenters.size(); ++j) {
        trees.emplace_back(barycenters[j]);
        clusteringAssignment.emplace_back(j);
      }
      numInputs = trees.size();
    }
    // - Declare VTK arrays
    vtkSmartPointer<vtkUnstructuredGrid> vtkArcs
      = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Node fields
    vtkNew<vtkIntArray> criticalType{};
    criticalType->SetName("CriticalType");
    vtkNew<vtkFloatArray> persistenceNode{};
    persistenceNode->SetName("Persistence");
    vtkNew<vtkIntArray> clusterIDNode{};
    clusterIDNode->SetName("ClusterID");
    vtkNew<vtkIntArray> isDummyNode{};
    isDummyNode->SetName("isDummyNode");
    vtkNew<vtkIntArray> branchNodeID{};
    branchNodeID->SetName("BranchNodeID");
    vtkNew<vtkFloatArray> scalar{};
    scalar->SetName("Scalar");
    vtkNew<vtkIntArray> isImportantPairsNode{};
    isImportantPairsNode->SetName("isImportantPair");
    vtkNew<vtkIntArray> nodeID{};
    nodeID->SetName("NodeId");

    vtkNew<vtkIntArray> treeIDNode{};
    treeIDNode->SetName("TreeID");
    vtkNew<vtkIntArray> branchBaryNodeID{};
    branchBaryNodeID->SetName("BranchBaryNodeID");
    vtkNew<vtkIntArray> isInterpolatedTreeNode{};
    isInterpolatedTreeNode->SetName("isInterpolatedTree");

    vtkNew<vtkFloatArray> percentMatch{};
    percentMatch->SetName("PercentMatchNode");
    vtkNew<vtkFloatArray> persistenceBaryNode{};
    persistenceBaryNode->SetName("PersistenceBarycenter");

    // Arc fields
    vtkNew<vtkFloatArray> persistenceArc{};
    persistenceArc->SetName("Persistence");
    vtkNew<vtkIntArray> clusterIDArc{};
    clusterIDArc->SetName("ClusterID");
    vtkNew<vtkIntArray> isImportantPairsArc{};
    isImportantPairsArc->SetName("isImportantPair");
    vtkNew<vtkIntArray> isDummyArc{};
    isDummyArc->SetName("isDummyArc");
    vtkNew<vtkIntArray> branchID{};
    branchID->SetName("BranchID");
    vtkNew<vtkIntArray> upNodeId{};
    upNodeId->SetName("upNodeId");
    vtkNew<vtkIntArray> downNodeId{};
    downNodeId->SetName("downNodeId");

    vtkNew<vtkIntArray> treeIDArc{};
    treeIDArc->SetName("TreeID");
    vtkNew<vtkIntArray> branchBaryID{};
    branchBaryID->SetName("BranchBaryID");
    vtkNew<vtkIntArray> isInterpolatedTreeArc{};
    isInterpolatedTreeArc->SetName("isInterpolatedTree");

    vtkNew<vtkFloatArray> percentMatchArc{};
    percentMatchArc->SetName("PercentMatchArc");
    vtkNew<vtkFloatArray> persistenceBaryArc{};
    persistenceBaryArc->SetName("PersistenceBarycenter");

    // Segmentation
    // vtkSmartPointer<vtkPoints> segPoints = vtkSmartPointer<vtkPoints>::New();
    // vtkSmartPointer<vtkUnstructuredGrid> segCells=
    // vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkAppendFilter> appendFilter
      = vtkSmartPointer<vtkAppendFilter>::New();

    // Internal data
    int cellCount = 0;
    bool foundOneInterpolatedTree = false;
    nodeCorr = std::vector<std::vector<SimplexId>>(numInputs);
    clusterShift = std::vector<double>(NumberOfBarycenters, 0);
    allBaryPercentMatch = std::vector<std::vector<float>>(NumberOfBarycenters);

    // --------------------------------------------------------
    // Iterate through all clusters
    // --------------------------------------------------------
    if(verbose > 0)
      printMsg("Iterate through all clusters", debug::Priority::VERBOSE);
    for(int c = 0; c < NumberOfBarycenters; ++c) {

      // Get radius
      if(verbose > 1)
        printMsg("// Get radius", debug::Priority::VERBOSE);
      double delta_max = std::numeric_limits<double>::lowest();
      int noSample = 0 + noSampleOffset;
      for(int i = 0; i < numInputsOri; ++i) {
        delta_max = std::max(
          (std::get<3>(allBounds[i]) - std::get<2>(allBounds[i])), delta_max);
        delta_max = std::max(
          (std::get<1>(allBounds[i]) - std::get<0>(allBounds[i])), delta_max);
        if(clusteringAssignment[i] != c)
          continue;
        noSample += 1;
      }
      double radius = delta_max * 2 * DimensionSpacing;
      int iSample = 0 + iSampleOffset;

      // Double line attributes
      double prevXMax = 0;
      std::vector<double> allPrevXMax;
      double prevYMax = std::numeric_limits<double>::lowest();

      // ------------------------------------------
      // Iterate through all trees of this cluster
      // ------------------------------------------
      if(verbose > 0)
        printMsg("Iterate through all trees of this cluster",
                 debug::Priority::VERBOSE);
      for(int i = 0; i < numInputs; ++i) {
        if(clusteringAssignment[i] != c)
          continue;

        // Get is interpolated tree (temporal subsampling)
        bool isInterpolatedTree = false;
        if(interpolatedTrees.size() != 0)
          isInterpolatedTree = interpolatedTrees[i];
        foundOneInterpolatedTree |= isInterpolatedTree;

        // Get branching
        if(verbose > 1)
          printMsg("// Get branching", debug::Priority::VERBOSE);
        std::vector<ftm::idNode> treeBranching;
        std::vector<int> treeBranchingID;
        trees[i]->getTreeBranching(treeBranching, treeBranchingID);

        // Get shift
        if(verbose > 1)
          printMsg("// Get shift", debug::Priority::VERBOSE);
        double angle = 360.0 / noSample * iSample;
        double pi = 3.14159265359;
        double diff_x = 0, diff_y = 0;
        switch(ShiftMode) {
          case 0: // Star
            diff_x = radius * std::cos(-1 * angle * pi / 180) + clusterShift[c];
            diff_y = radius * std::sin(-1 * angle * pi / 180);
            break;
          case 1: // Star Barycenter
            diff_x = clusterShift[c];
            diff_y = 0;
            break;
          case 2: // Line
            diff_x = prevXMax + radius;
            break;
          case 3: // Double Line
            diff_x = prevXMax + radius;
            if(i >= numInputs / 2) {
              diff_y = -(prevYMax + radius / 2);
              diff_x = allPrevXMax[i - int(numInputs / 2)] + radius;
            } else
              allPrevXMax.emplace_back(prevXMax);
            break;
          default:
            break;
        }
        iSample += 1;

        // Get dimension shift
        if(verbose > 1)
          printMsg("// Get dimension shift", debug::Priority::VERBOSE);
        double diff_z = -std::get<4>(allBounds[i]);
        if(/*not clusteringOutput and*/ not PlanarLayout) {
          if(DimensionToShift != 0) { // is not X
            if(DimensionToShift == 2) // is Z
              diff_z = diff_x;
            else if(DimensionToShift == 1) // is Y
              diff_y = diff_x;
            diff_x = -std::get<0>(allBounds[i]);
          }
        }

        // Planar layout
        if(verbose > 1)
          printMsg("// Planar Layout", debug::Priority::VERBOSE);
        std::vector<float> layout;
        if(PlanarLayout) {
          double refPersistence;
          if(clusteringOutput)
            refPersistence = barycenters[0]->getNodePersistence<dataType>(
              barycenters[0]->getRoot());
          else
            refPersistence
              = trees[0]->getNodePersistence<dataType>(trees[0]->getRoot());
          layout = treePlanarLayout<dataType>(
            trees[i], allBaryBounds[c], refPersistence);
        }

        // Internal arrays
        if(verbose > 1)
          printMsg("// Internal arrays", debug::Priority::VERBOSE);
        int cptNode = 0;
        nodeCorr[i] = std::vector<SimplexId>(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> treeSimplexId(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> treeDummySimplexId(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> layoutCorr(trees[i]->getNumberOfNodes());
        std::vector<ftm::idNode> treeMatching(trees[i]->getNumberOfNodes(), -1);
        if(clusteringOutput)
          for(auto match : outputMatchingBarycenter[c][i])
            treeMatching[std::get<1>(match)] = std::get<0>(match);
        // _ m[i][j] contains the node in treesOri[j] matched to the node i in
        // the barycenter
        std::vector<std::vector<ftm::idNode>> baryMatching(
          trees[i]->getNumberOfNodes(),
          std::vector<ftm::idNode>(numInputsOri, -1));
        if(ShiftMode == 1) {
          for(size_t j = 0; j < outputMatchingBarycenter[c].size(); ++j)
            for(auto match : outputMatchingBarycenter[c][j])
              baryMatching[std::get<0>(match)][j] = std::get<1>(match);
          allBaryPercentMatch[c]
            = std::vector<float>(trees[i]->getNumberOfNodes(), 100);
        }

        // ----------------------------
        // Tree traversal
        // ----------------------------
        if(verbose > 0)
          printMsg("// Tree traversal", debug::Priority::VERBOSE);
        std::queue<idNode> queue;
        queue.emplace(trees[i]->getRoot());
        while(!queue.empty()) {
          idNode node = queue.front();
          queue.pop();
          idNode nodeOrigin = trees[i]->getNode(node)->getOrigin();

          // Push children to the queue
          if(verbose > 2)
            printMsg("// Push children to the queue", debug::Priority::VERBOSE);
          for(auto child : trees[i]->getChildren(node))
            queue.emplace(child);

          // --------------
          // Insert point
          // --------------
          if(verbose > 2)
            printMsg("// Get and insert point", debug::Priority::VERBOSE);
          int nodeMesh = -1;
          int nodeMeshTreeIndex = -1;
          double noMatched = 0.0;
          double point[3] = {0, 0, 0};
          if(ShiftMode == 1) { // Star barycenter
            for(int j = 0; j < numInputsOri; ++j) {
              if((int)baryMatching[node][j] != -1) {
                nodeMesh = treesNodeCorrMesh[j][baryMatching[node][j]];
                double *pointTemp
                  = treesNodes[j]->GetPoints()->GetPoint(nodeMesh);
                for(int k = 0; k < 3; ++k)
                  point[k] += pointTemp[k];
                noMatched += 1;
                nodeMeshTreeIndex = j;
              }
            }
            for(int k = 0; k < 3; ++k)
              point[k] /= noMatched;
          } else if(not isInterpolatedTree) {
            nodeMesh = treesNodeCorrMesh[i][node];
            double *pointP = treesNodes[i]->GetPoints()->GetPoint(nodeMesh);
            for(int p = 0; p < 3; ++p)
              point[p] = pointP[p];
          }
          if(PlanarLayout) {
            layoutCorr[node] = cptNode;
            point[0] = layout[cptNode];
            point[1] = layout[cptNode + 1];
            point[2] = 0;
            cptNode += 2;
          } else
            point[2] += diff_z;
          point[0] += diff_x;
          point[1] += diff_y;

          // Get x Max and y Min for next iteration if needed (double line mode)
          prevXMax = std::max(prevXMax, point[0]);
          if(ShiftMode == 3) { // Double line
            if(i < numInputs / 2)
              prevYMax = std::max(prevYMax, point[1]);
            if(i == int(numInputs / 2) - 1)
              prevXMax = 0;
          }

          // TODO too many dummy nodes are created
          bool dummyNode = PlanarLayout and not BranchDecompositionPlanarLayout
                           and !trees[i]->isRoot(node)
            /*and !isLeaf(trees[i], node)
            and isBranchOrigin(trees[i], node)*/
            ;
          if(dummyNode)
            treeDummySimplexId[node] = points->InsertNextPoint(
              point); // will be modified when processing son
          SimplexId nextPointId = points->InsertNextPoint(point);
          treeSimplexId[node] = nextPointId;
          nodeCorr[i][node] = nextPointId;
          if(dummyNode)
            nodeCorr[i][node] = treeDummySimplexId[node];

          // --------------
          // Insert cell connecting parent
          // --------------
          if(verbose > 2)
            printMsg("// Add cell connecting parent", debug::Priority::VERBOSE);
          if(!trees[i]->isRoot(node)) {
            vtkIdType pointIds[2];
            pointIds[0] = treeSimplexId[node];

            ftm::idNode nodeParent = trees[i]->getParentSafe(node);
            // TODO too many dummy cells are created
            bool dummyCell = PlanarLayout
                             and not BranchDecompositionPlanarLayout
                             and treeBranching[node] == nodeParent
                             and !trees[i]->isRoot(nodeParent);
            if(PlanarLayout and BranchDecompositionPlanarLayout) {
              pointIds[1] = treeSimplexId[treeBranching[node]];
            } else if(dummyCell) {
              double dummyPoint[3]
                = {point[0], layout[layoutCorr[nodeParent] + 1] + diff_y, 0.};
              SimplexId dummyPointId = treeDummySimplexId[nodeParent];
              points->SetPoint(dummyPointId, dummyPoint);
              vtkIdType dummyPointIds[2];
              dummyPointIds[0] = dummyPointId;
              dummyPointIds[1] = treeSimplexId[nodeParent];
              vtkArcs->InsertNextCell(VTK_LINE, 2, dummyPointIds);
              pointIds[1] = dummyPointId;
            } else
              pointIds[1] = treeSimplexId[nodeParent];

            vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);

            // --------------
            // Arc field
            // --------------
            int toAdd = (dummyCell ? 2 : 1);
            for(int toAddT = 0; toAddT < toAdd; ++toAddT) {
              // Add arc matching percentage
              if(ShiftMode == 1) // Star Barycenter
                percentMatchArc->InsertNextTuple1(
                  allBaryPercentMatch[c][allBaryBranching[c][node]]);

              // Add branch bary ID
              if(verbose > 2)
                printMsg(
                  "// Push arc bary branch id", debug::Priority::VERBOSE);
              if(clusteringOutput and ShiftMode != 1) {
                int tBranchID = -1;
                if(treeMatching[node] >= 0
                   and treeMatching[node] < allBaryBranchingID[c].size()) {
                  tBranchID = allBaryBranchingID[c][treeMatching[node]];
                  if(not trees[i]->isLeaf(node)
                     and treeMatching[nodeOrigin] >= 0
                     and treeMatching[nodeOrigin]
                           < allBaryBranchingID[c].size())
                    tBranchID = allBaryBranchingID[c][treeMatching[nodeOrigin]];
                }
                branchBaryID->InsertNextTuple1(tBranchID);
              }

              // Add branch ID
              int tBranchID = treeBranchingID[node];
              branchID->InsertNextTuple1(tBranchID);

              // Add up and down nodeId
              upNodeId->InsertNextTuple1(treeSimplexId[nodeParent]);
              downNodeId->InsertNextTuple1(treeSimplexId[node]);

              // Add arc persistence
              if(verbose > 2)
                printMsg("// Push arc persistence", debug::Priority::VERBOSE);
              idNode nodeToGetPers = treeBranching[node];
              if(PlanarLayout and BranchDecompositionPlanarLayout)
                nodeToGetPers = node;
              double persToAdd
                = trees[i]->getNodePersistence<dataType>(nodeToGetPers);
              persistenceArc->InsertNextTuple1(persToAdd);

              // Add arc persistence barycenter
              if(clusteringOutput and ShiftMode != 1) {
                idNode nodeToGet = treeBranching[node];
                if(PlanarLayout and BranchDecompositionPlanarLayout)
                  nodeToGet = node;
                if(treeMatching[nodeToGet] >= 0
                   and treeMatching[nodeToGet] < allBaryBranchingID[c].size())
                  persistenceBaryArc->InsertNextTuple1(
                    barycenters[c]->getNodePersistence<dataType>(
                      treeMatching[nodeToGet]));
              }

              // Add arc cluster ID
              clusterIDArc->InsertNextTuple1(clusteringAssignment[i]);

              // Add arc tree ID
              treeIDArc->InsertNextTuple1(i);

              // Add isImportantPair
              bool isImportant = false;
              idNode nodeToGetImportance = treeBranching[node];
              if(PlanarLayout and BranchDecompositionPlanarLayout)
                nodeToGetImportance = node;
              isImportant = trees[i]->isImportantPair<dataType>(
                nodeToGetImportance, ImportantPairs);
              isImportantPairsArc->InsertNextTuple1(isImportant);

              // Add isDummyArc
              bool isDummy = toAdd == 2 and toAddT == 0;
              isDummyArc->InsertNextTuple1(isDummy);

              // Add isInterpolatedTree
              isInterpolatedTreeArc->InsertNextTuple1(isInterpolatedTree);

              cellCount++;
            }
          }

          // --------------
          // Node field
          // --------------
          int toAdd = (dummyNode ? 2 : 1);
          for(int toAddT = 0; toAddT < toAdd; ++toAddT) {
            // Add node id
            nodeID->InsertNextTuple1(treeSimplexId[node]);

            // Add node scalar
            scalar->InsertNextTuple1(trees[i]->getValue<dataType>(node));

            // Add criticalType
            if(verbose > 2)
              printMsg("// Add criticalType", debug::Priority::VERBOSE);
            int criticalTypeT = -1;
            if(not isInterpolatedTree) {
              if(ShiftMode == 1) {
                if(nodeMeshTreeIndex != -1) {
                  auto array
                    = treesNodes[nodeMeshTreeIndex]->GetPointData()->GetArray(
                      "CriticalType");
                  criticalTypeT = array->GetTuple1(nodeMesh);
                }
              } else {
                auto array
                  = treesNodes[i]->GetPointData()->GetArray("CriticalType");
                criticalTypeT = array->GetTuple1(nodeMesh);
              }
            } else {
              // TODO
            }
            criticalType->InsertNextTuple1(criticalTypeT);

            // Add node matching percentage
            if(ShiftMode == 1) { // Star Barycenter
              float percentMatchT = noMatched * 100 / numInputs;
              percentMatch->InsertNextTuple1(percentMatchT);
              allBaryPercentMatch[c][node] = percentMatchT;
            }

            // Add node branch bary id
            if(verbose > 2)
              printMsg("// Add node bary branch id", debug::Priority::VERBOSE);
            if(clusteringOutput and ShiftMode != 1) {
              int tBranchID = -1;
              if(treeMatching[node] >= 0
                 and treeMatching[node] < allBaryBranchingID[c].size()) {
                tBranchID = allBaryBranchingID[c][treeMatching[node]];
                if(not trees[i]->isLeaf(node) and treeMatching[nodeOrigin] >= 0
                   and treeMatching[nodeOrigin] < allBaryBranchingID[c].size())
                  tBranchID = allBaryBranchingID[c][treeMatching[nodeOrigin]];
              }
              branchBaryNodeID->InsertNextTuple1(tBranchID);
            }

            // Add node branch bary id
            int tBranchID = treeBranchingID[node];
            if(not trees[i]->isLeaf(node))
              tBranchID = treeBranchingID[nodeOrigin];
            branchNodeID->InsertNextTuple1(tBranchID);

            // Add node persistence
            if(verbose > 2)
              printMsg("// Push node persistence", debug::Priority::VERBOSE);
            persistenceNode->InsertNextTuple1(
              trees[i]->getNodePersistence<dataType>(node));

            // Add node persistence barycenter
            if(clusteringOutput and ShiftMode != 1)
              if(treeMatching[node] >= 0
                 and treeMatching[node] < allBaryBranchingID[c].size())
                persistenceBaryNode->InsertNextTuple1(
                  barycenters[c]->getNodePersistence<dataType>(
                    treeMatching[node]));

            // Add node clusterID
            clusterIDNode->InsertNextTuple1(clusteringAssignment[i]);

            // Add arc tree ID
            treeIDNode->InsertNextTuple1(i);

            // Add isDummyNode
            bool isDummy
              = toAdd == 2 and toAddT == 1 and !trees[i]->isRoot(node);
            isDummyNode->InsertNextTuple1(isDummy);

            // Add isInterpolatedTree
            isInterpolatedTreeNode->InsertNextTuple1(isInterpolatedTree);

            // Add isImportantPair
            bool isImportant = false;
            isImportant
              = trees[i]->isImportantPair<dataType>(node, ImportantPairs);
            isImportantPairsNode->InsertNextTuple1(isImportant);
          }

          if(verbose > 2)
            printMsg("end loop", debug::Priority::VERBOSE);
        } // end tree traversal

        // --------------
        // Manage segmentation
        // --------------
        if(verbose > 1)
          printMsg("// Shift segmentation", debug::Priority::VERBOSE);
        if(OutputSegmentation and not PlanarLayout) {
          auto iTreesSegmentationCopy
            = vtkSmartPointer<vtkUnstructuredGrid>::New();
          iTreesSegmentationCopy->DeepCopy(treesSegmentation[i]);
          auto iVkOutputSegmentationTemp
            = vtkUnstructuredGrid::SafeDownCast(iTreesSegmentationCopy);
          for(int p = 0;
              p < iVkOutputSegmentationTemp->GetPoints()->GetNumberOfPoints();
              ++p) {
            double *point = iVkOutputSegmentationTemp->GetPoints()->GetPoint(p);
            point[0] += diff_x;
            point[1] += diff_y;
            point[2] += diff_z;
            iVkOutputSegmentationTemp->GetPoints()->SetPoint(p, point);
          }
          appendFilter->AddInputData(iVkOutputSegmentationTemp);
        }
      }
      if(c < NumberOfBarycenters - 1)
        clusterShift[c + 1] = radius * 4 + clusterShift[c];
    }

    // --- Add VTK arrays to output
    if(verbose > 1)
      printMsg("// Add VTK arrays to output", debug::Priority::VERBOSE);
    // Manage node output
    vtkOutputNode->SetPoints(points);
    vtkOutputNode->GetPointData()->AddArray(criticalType);
    vtkOutputNode->GetPointData()->AddArray(persistenceNode);
    vtkOutputNode->GetPointData()->AddArray(clusterIDNode);
    vtkOutputNode->GetPointData()->AddArray(treeIDNode);
    vtkOutputNode->GetPointData()->AddArray(isDummyNode);
    vtkOutputNode->GetPointData()->AddArray(branchNodeID);
    vtkOutputNode->GetPointData()->AddArray(nodeID);
    vtkOutputNode->GetPointData()->AddArray(isImportantPairsNode);
    if(not BranchDecompositionPlanarLayout)
      vtkOutputNode->GetPointData()->AddArray(scalar);
    if(clusteringOutput and ShiftMode != 1) {
      vtkOutputNode->GetPointData()->AddArray(branchBaryNodeID);
      vtkOutputNode->GetPointData()->AddArray(persistenceBaryNode);
    }
    if(foundOneInterpolatedTree)
      vtkOutputNode->GetPointData()->AddArray(isInterpolatedTreeNode);
    if(ShiftMode == 1) // Star Barycenter
      vtkOutputNode->GetPointData()->AddArray(percentMatch);

    // Manage arc output
    vtkArcs->SetPoints(points);
    vtkArcs->GetCellData()->AddArray(persistenceArc);
    vtkArcs->GetCellData()->AddArray(clusterIDArc);
    vtkArcs->GetCellData()->AddArray(treeIDArc);
    vtkArcs->GetCellData()->AddArray(isImportantPairsArc);
    vtkArcs->GetCellData()->AddArray(isDummyArc);
    vtkArcs->GetCellData()->AddArray(branchID);
    vtkArcs->GetCellData()->AddArray(upNodeId);
    vtkArcs->GetCellData()->AddArray(downNodeId);
    if(clusteringOutput and ShiftMode != 1) {
      vtkArcs->GetCellData()->AddArray(branchBaryID);
      vtkArcs->GetCellData()->AddArray(persistenceBaryArc);
    }
    if(foundOneInterpolatedTree)
      vtkArcs->GetCellData()->AddArray(isInterpolatedTreeArc);
    if(ShiftMode == 1) // Star Barycenter
      vtkArcs->GetCellData()->AddArray(percentMatchArc);
    vtkOutputArc->ShallowCopy(vtkArcs);

    // Manage segmentation output
    if(OutputSegmentation and not PlanarLayout) {
      appendFilter->Update();
      vtkOutputSegmentation->ShallowCopy(appendFilter->GetOutput());
    }

    //
    if(ShiftMode == 1) // Star Barycenter
      trees = treesOri;
  }

  // ==========================================================================
  // Bounds Utils
  // ==========================================================================
  std::tuple<double, double, double, double, double, double>
    getRealBounds(vtkUnstructuredGrid *treeNodes,
                  FTMTree_MT *tree,
                  std::vector<int> &nodeCorrT) {
    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    double y_max = std::numeric_limits<double>::lowest();
    double z_max = std::numeric_limits<double>::lowest();
    std::queue<idNode> queue;
    queue.emplace(tree->getRoot());
    while(!queue.empty()) {
      idNode node = queue.front();
      queue.pop();
      double *point = treeNodes->GetPoints()->GetPoint(nodeCorrT[node]);
      x_min = std::min(x_min, point[0]);
      x_max = std::max(x_max, point[0]);
      y_min = std::min(y_min, point[1]);
      y_max = std::max(y_max, point[1]);
      z_min = std::min(z_min, point[2]);
      z_max = std::max(z_max, point[2]);
      for(auto child : tree->getChildren(node))
        queue.emplace(child);
    }
    return std::make_tuple(x_min, x_max, y_min, y_max, z_min, z_max);
  }

  std::tuple<double, double, double, double, double, double>
    getRealBounds(vtkUnstructuredGrid *treeNodes, FTMTree_MT *tree) {
    std::vector<int> nodeCorrT(tree->getNumberOfNodes());
    for(size_t i = 0; i < nodeCorrT.size(); ++i)
      nodeCorrT[i] = i;

    return getRealBounds(treeNodes, tree, nodeCorrT);
  }
};

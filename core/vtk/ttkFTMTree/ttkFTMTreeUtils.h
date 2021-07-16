/// \ingroup vtk
/// \class ttk::ttkFTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#ifndef _TTKFTMTREEUTILS_H
#define _TTKFTMTREEUTILS_H

#pragma once

#include <FTMTree.h>

#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

namespace ttk {
  namespace ftm {

    template <class dataType>
    MergeTree makeTree(vtkUnstructuredGrid *treeNodes,
                       vtkUnstructuredGrid *treeArcs) {
      // Init Scalars
      Scalars scalars;
      vtkSmartPointer<vtkDataArray> nodesScalar
        = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar
      scalars.size = nodesScalar->GetNumberOfTuples();
      scalars.values = ttkUtils::GetVoidPointer(nodesScalar);

      // Init Tree
      Params params;
      params.treeType = Join_Split;
      MergeTree mergeTree(scalars, params);

      // Add Nodes
      vtkSmartPointer<vtkDataArray> nodesId
        = treeNodes->GetPointData()->GetArray("NodeId"); // 0: NodeId
      vtkIdType nodesNumTuples = nodesId->GetNumberOfTuples();
      for(vtkIdType i = 0; i < nodesNumTuples; ++i) {
        mergeTree.tree.makeNode(i);
      }

      // Add Arcs
      vtkSmartPointer<vtkDataArray> arcsUp
        = treeArcs->GetCellData()->GetArray("upNodeId"); // 1: upNodeId
      vtkSmartPointer<vtkDataArray> arcsDown
        = treeArcs->GetCellData()->GetArray("downNodeId"); // 2: downNodeId
      vtkIdType arcsNumTuples = arcsUp->GetNumberOfTuples();
      std::set<std::tuple<double, double>> added_arcs; // Avoid duplicates
      for(vtkIdType i = 0; i < arcsNumTuples; ++i) {
        double downId = arcsDown->GetTuple1(i);
        double upId = arcsUp->GetTuple1(i);
        auto it = added_arcs.find(std::make_tuple(downId, upId));
        if(it == added_arcs.end()) { // arc not added yet
          mergeTree.tree.makeSuperArc(downId, upId); // (down, Up)
          added_arcs.insert(std::make_tuple(downId, upId));
        }
      }

      // Manage inconsistent arcs
      // manageInconsistentArcsMultiParent(tree);;

      // Remove self link
      // removeSelfLink(tree);

      // tree->printTree2();

      return mergeTree;
    }

  } // namespace ftm
} // namespace ttk

#endif

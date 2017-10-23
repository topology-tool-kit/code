/// \ingroup base
/// \class ttk::ContourForestsTree
/// \author Charles Gueunet <charles.gueunet@lip6.fr>
/// \date June 2016.
///
///\brief TTK processing package that efficiently computes the contour tree of
/// scalar data and more (data segmentation, topological simplification,
/// persistence diagrams, persistence curves, etc.).
///
///\param dataType Data type of the input scalar field (char, float,
/// etc.).
///
/// \b Related \b publication \n
/// "Contour Forests: Fast Multi-threaded Augmented Contour Trees" \n
/// Charles Gueunet, Pierre Fortin, Julien Jomier, Julien Tierny \n
/// Proc. of IEEE LDAV 2016.
///
/// \sa vtkContourForests.cpp %for a usage example.

#ifndef _CONTOURFORESTSTEMPLATE_H
#define _CONTOURFORESTSTEMPLATE_H

#include "ContourForests.h"

// TODO
// Template wrapper

// ------------------- Contour Forests

// Process
// {

template <typename scalarType>
int ContourForests::build()
{

#ifdef TTK_ENABLE_OPENMP
        omp_set_num_threads(parallelParams_.nbThreads);
#endif

   DebugTimer timerTOTAL;

   // -----------
   // Paramemters
   // -----------
   initTreeType();
   initNbScalars();
   initNbPartitions();
   initSoS();

   if (params_->debugLevel >= 2) {
      // print params:
      cout << "threads :" << static_cast<unsigned>(parallelParams_.nbThreads) << endl;
      cout << "partitions : " <<static_cast<unsigned>(parallelParams_.nbPartitions) << endl;
      if(params_->simplifyThreshold){
         cout << "simplify method : " << params_->simplifyMethod << endl;
         cout << "simplify thresh.: " << params_->simplifyThreshold << endl;
      }
   }
   printDebug(timerTOTAL, "Initialization                   ");

   // ---------
   // Sort Step
   // ---------

   DebugTimer timerSort;
   sortInput<scalarType>();
   printDebug(timerSort, "Sort scalars (+mirror)           ");

   // -------------------
   // Interface & Overlap
   // -------------------

   DebugTimer timerInitOverlap;
   initInterfaces();
   initOverlap();
   if(params_->debugLevel > 3){
        for (idInterface i = 0; i < parallelParams_.nbInterfaces; i++) {
            cout << "interface : " << static_cast<unsigned>(i);
            cout << " seed : " << parallelData_.interfaces[i].getSeed();
            cout << endl;
        }
   }
   printDebug(timerInitOverlap, "Interface and overlap init.      ");

   // -----------------------
   // Allocate parallel trees
   // -----------------------

   DebugTimer timerAllocPara;
   // Union find vector for each partition
   vector<vector<ExtendedUnionFind *>> vect_baseUF_JT(parallelParams_.nbPartitions),
                                       vect_baseUF_ST(parallelParams_.nbPartitions);
   const idVertex &resSize = (scalars_->size / parallelParams_.nbPartitions) / 10;

   parallelData_.trees.clear();
   parallelData_.trees.reserve(parallelParams_.nbPartitions);

   for (idPartition tree = 0; tree < parallelParams_.nbPartitions; ++tree) {
      // Tree array initialization
      parallelData_.trees.emplace_back(params_,mesh_,scalars_, tree);
   }

#pragma omp parallel for num_threads(parallelParams_.nbPartitions) schedule(static)
   for (idPartition tree = 0; tree < parallelParams_.nbPartitions; ++tree) {
      // Tree array initialization
      parallelData_.trees[tree].flush();

      // UF-array reserve
      vect_baseUF_JT[tree].resize(scalars_->size);
      vect_baseUF_ST[tree].resize(scalars_->size);

      // Statistical reserve
      parallelData_.trees[tree].jt_->treeData_.nodes.reserve(resSize);
      parallelData_.trees[tree].jt_->treeData_.superArcs.reserve(resSize);
      parallelData_.trees[tree].st_->treeData_.nodes.reserve(resSize);
      parallelData_.trees[tree].st_->treeData_.superArcs.reserve(resSize);

   }
   printDebug(timerAllocPara, "Parallel allocations             ");

   // -------------------------
   // Build trees in partitions
   // -------------------------

   DebugTimer timerbuild;
   parallelBuild<scalarType>(vect_baseUF_JT, vect_baseUF_ST);

   if (params_->debugLevel >= 4) {
      if (params_->treeType == TreeType::Contour) {
         for (idPartition i = 0; i < parallelParams_.nbPartitions; i++) {
            cout << i << " :" << endl;
            parallelData_.trees[i].printTree2();
            cout << "-----" << endl;
         }
       } else {
          for (idPartition i = 0; i < parallelParams_.nbPartitions; i++) {
             cout << i << " jt:" << endl;
             parallelData_.trees[i].jt_->printTree2();
             cout << i << " st:" << endl;
             parallelData_.trees[i].st_->printTree2();
             cout << "-----" << endl;
         }
       }
   }

   printDebug(timerbuild, "ParallelBuild                    ");

   // --------------------
   // Stitching partitions
   // --------------------

   DebugTimer timerZip;
   if (parallelParams_.partitionNum == -1 && parallelParams_.nbPartitions > 1 ) {
      stitch();
      for (idPartition p = 0; p < parallelParams_.nbPartitions; ++p) {
         parallelData_.trees[p].parallelInitNodeValence(parallelParams_.nbThreads);
      }
   }

   if (params_->debugLevel >= 4) {
      printVectCT();
   }

   printDebug(timerZip, "Stitch                           ");

   // -------------------------------------------------
   // Unification : create one tree from stitched trees
   // -------------------------------------------------

   DebugTimer timerUnify;
   if (params_->treeType == TreeType::Contour) {
       if(parallelParams_.partitionNum >= 0){
          if (parallelParams_.partitionNum > parallelParams_.nbInterfaces) {
             clone(&parallelData_.trees[parallelParams_.nbPartitions - 1]);
          } else {
             clone(&parallelData_.trees[parallelParams_.partitionNum]);
          }
       } else if (parallelParams_.nbPartitions == 1) {
          clone(&parallelData_.trees[0]);
       } else {
          unify();
          // for global simlify
          parallelInitNodeValence(parallelParams_.nbThreads);
       }
   } else {
       if(parallelParams_.partitionNum >= 0){
           if(parallelParams_.partitionNum > parallelParams_.nbInterfaces){
              jt_->clone(parallelData_.trees[parallelParams_.nbInterfaces].jt_);
              st_->clone(parallelData_.trees[parallelParams_.nbInterfaces].st_);
           } else {
              jt_->clone(parallelData_.trees[parallelParams_.partitionNum].jt_);
              st_->clone(parallelData_.trees[parallelParams_.partitionNum].st_);
           }
       } else if (parallelParams_.nbPartitions == 1) {
          jt_->clone(parallelData_.trees[0].jt_);
          st_->clone(parallelData_.trees[0].st_);
       } else {
          unify();
          jt_->parallelInitNodeValence(parallelParams_.nbThreads);
          st_->parallelInitNodeValence(parallelParams_.nbThreads);
       }
   }

   printDebug(timerUnify, "Create Contour tree              ");

   // -------------------
   // Simplification step
   // -------------------

   if (params_->treeType == TreeType::Contour && parallelParams_.partitionNum == -1 &&
       params_->simplifyThreshold) {
      DebugTimer timerGlobalSimplify;
      idEdge simplifed = globalSimplify<scalarType>(-1, nullVertex);
      if(params_->debugLevel >=1){
         printDebug(timerGlobalSimplify, "Simplify Contour tree            ");
         cout << " ( " << simplifed << " pairs merged )" << endl;
      }
   }

   printDebug(timerTOTAL, "TOTAL                            ");

   // ------------------------------
   // Debug print and memory reclaim
   // ------------------------------

    if (params_->debugLevel >= 5)  {
        if(params_->treeType == TreeType::Contour)
           printTree2();
        else {
           cout << "JT :" << endl;
           jt_->printTree2();
           cout << "ST :" << endl;
           st_->printTree2();
        }
    } else if (params_->debugLevel > 2) {
        if(params_->treeType == TreeType::Contour)
           cout << "max node : " << getNumberOfNodes() << endl;
        else {
            cout << "JT max node : " << jt_->getNumberOfNodes() << endl;
            cout << "ST max node : " << st_->getNumberOfNodes() << endl;
        }
    }

   if (params_->treeType == TreeType::Contour) {
      updateSegmentation();
   } else {
      jt_->updateSegmentation();
      st_->updateSegmentation();
   }

   // reclaim memory
   {
      for (idPartition tree = 0; tree < parallelParams_.nbPartitions; ++tree) {
         parallelData_.trees[tree].jt_->treeData_.nodes.shrink_to_fit();
         parallelData_.trees[tree].jt_->treeData_.superArcs.shrink_to_fit();
         parallelData_.trees[tree].st_->treeData_.nodes.shrink_to_fit();
         parallelData_.trees[tree].st_->treeData_.superArcs.shrink_to_fit();
      }
      // Not while arc segmentation depends on vector in partitions
      //parallelData_.interfaces.clear();
      //parallelData_.trees.clear();
   }

   cout << "Tree computed ..." << endl;

   return 0;
}

template <typename scalarType>
int ContourForests::parallelBuild(vector<vector<ExtendedUnionFind *>> &vect_baseUF_JT,
                                  vector<vector<ExtendedUnionFind *>> &vect_baseUF_ST)
{
   vector<float> timeSimplify(parallelParams_.nbPartitions, 0);
   vector<float> speedProcess(parallelParams_.nbPartitions*2, 0);
#ifdef TTK_CONTOUR_FORESTS_ENABLE_PARALLEL_SIMPLIFY
   idEdge nbPairMerged = 0;
#endif

#ifdef TTK_ENABLE_OPENMP
   omp_set_nested(1);
#endif

//cout << "NO PARALLEL DEBUG MODE" << endl;
#pragma omp parallel for num_threads(parallelParams_.nbPartitions) schedule(static)
   for (idPartition i = 0; i < parallelParams_.nbPartitions; ++i) {
      DebugTimer timerMergeTree;

      // ------------------------------------------------------
      // Skip partition that are not asked to compute if needed
      // ------------------------------------------------------

      if (parallelParams_.partitionNum != -1 && parallelParams_.partitionNum != i)
         continue;

      // ------------------------------------------------------
      // Retrieve boundary & overlap list for current partition
      // ------------------------------------------------------

      tuple<idVertex, idVertex> rangeJT = getJTRange(i);
      tuple<idVertex, idVertex> rangeST = getSTRange(i);
      tuple<idVertex, idVertex> seedsPos = getSeedsPos(i);
      tuple<vector<idVertex>, vector<idVertex>> overlaps = getOverlaps(i);
      const idVertex &partitionSize = abs(get<0>(rangeJT) - get<1>(rangeJT)) +
                                      get<0>(overlaps).size() + get<1>(overlaps).size();

      // ---------------
      // Build JT and ST
      // ---------------

#pragma omp parallel sections num_threads(2) if (parallelParams_.lessPartition)
      {

          // if less partition : we built JT and ST in parallel

#pragma omp section
        {
            if(params_->treeType == TreeType::Join
              || params_->treeType == TreeType::Contour
              || params_->treeType == TreeType::JoinAndSplit)
            {
                DebugTimer timerSimplify;
                DebugTimer timerBuild;
                parallelData_.trees[i].getJoinTree()->build(vect_baseUF_JT[i],
                        get<0>(overlaps), get<1>(overlaps),
                        get<0>(rangeJT), get<1>(rangeJT),
                        get<0>(seedsPos), get<1>(seedsPos)
                        );
                speedProcess[i] = partitionSize / timerBuild.getElapsedTime();

#ifdef TTK_CONTOUR_FORESTS_ENABLE_PARALLEL_SIMPLIFY
                timerSimplify.reStart();
                const idEdge tmpMerge =
                    parallelData_.trees[i].getJoinTree()->localSimplify<scalarType>(
                            get<0>(seedsPos), get<1>(seedsPos));
#pragma omp atomic update
                timeSimplify[i] += timerSimplify.getElapsedTime();
#pragma omp atomic update
                nbPairMerged += tmpMerge;
#endif
            }
        }


#pragma omp section
        {
            if(params_->treeType == TreeType::Split
              || params_->treeType == TreeType::Contour
              || params_->treeType == TreeType::JoinAndSplit)
            {
                DebugTimer timerSimplify;
                DebugTimer timerBuild;
                parallelData_.trees[i].getSplitTree()->build(vect_baseUF_ST[i],
                        get<1>(overlaps), get<0>(overlaps),
                        get<0>(rangeST), get<1>(rangeST),
                        get<0>(seedsPos), get<1>(seedsPos)
                        );
                speedProcess[parallelParams_.nbPartitions + i] =
                    partitionSize / timerBuild.getElapsedTime();

#ifdef TTK_CONTOUR_FORESTS_ENABLE_PARALLEL_SIMPLIFY
                timerSimplify.reStart();
                const idEdge tmpMerge =
                    parallelData_.trees[i].getSplitTree()->localSimplify<scalarType>(
                            get<0>(seedsPos), get<1>(seedsPos));
#pragma omp atomic update
                timeSimplify[i] += timerSimplify.getElapsedTime();
#pragma omp atomic update
                nbPairMerged += tmpMerge;
#endif
            }
        }


      }

      {
         stringstream mt;
         mt << "[ParallelBuild] Merge Tree " << static_cast<unsigned>(i)
            << " constructed in : " << timerMergeTree.getElapsedTime() << endl;
         dMsg(cout, mt.str(), infoMsg);
      }

      // Update segmentation of each arc if needed
      if ( params_->simplifyThreshold || params_->treeType != TreeType::Contour){
         DebugTimer timerUpdateSegm;
         parallelData_.trees[i].getJoinTree()->updateSegmentation();
         parallelData_.trees[i].getSplitTree()->updateSegmentation();

         if (params_->debugLevel >= 3) {
            cout << "Local MT : updated in " << timerUpdateSegm.getElapsedTime() << endl;
         }
      }

      // ---------------
      // Combine JT & ST
      // ---------------

      if (params_->treeType == TreeType::Contour) {
         DebugTimer timerCombine;

         // clone here if we do not want to destry original merge trees!
         auto *jt = parallelData_.trees[i].getJoinTree();
         auto *st = parallelData_.trees[i].getSplitTree();

         // Copy missing nodes of a tree to the other one
         // Maintain this traversal order for good insertion
         for (idNode t = 0; t < st->getNumberOfNodes(); ++t) {
            if (!st->getNode(t)->isHidden()) {
               // cout << "insert in jt : " << st->getNode(t)->getVertexId() << endl;
               jt->insertNode(st->getNode(t), true);
            }
         }
         // and vice versa
         for (idNode t = 0; t < jt->getNumberOfNodes(); ++t) {
            if (!jt->getNode(t)->isHidden()) {
               // cout << "insert in st : " << jt->getNode(t)->getVertexId() << endl;
               st->insertNode(jt->getNode(t), true);
            }
         }

         // debug print current JT / ST
         if (params_->debugLevel >= 6) {
            cout << "Local JT :" << endl;
            parallelData_.trees[i].getJoinTree()->printTree2();
            cout << "Local ST :" << endl;
            parallelData_.trees[i].getSplitTree()->printTree2();
            cout << "combine" << endl;
         }

         // Combine, destroy JT and ST to compute CT
         parallelData_.trees[i].combine(get<0>(seedsPos), get<1>(seedsPos));
         parallelData_.trees[i].updateSegmentation();

         if (params_->debugLevel > 2) {
            printDebug(timerCombine, "Trees combined   in    ");
         }

         // debug print CT
         if (params_->debugLevel >= 4) {
            parallelData_.trees[i].printTree2();
         }
      } else {
         if (params_->debugLevel >= 6) {
            cout << "Local JT :" << endl;
            parallelData_.trees[i].getJoinTree()->printTree2();
            cout << "Local ST :" << endl;
            parallelData_.trees[i].getSplitTree()->printTree2();
            cout << "combine" << endl;
         }
      }
   }

   // -------------------------------------
   // Print process speed and simplify info
   // -------------------------------------

   if(params_->debugLevel > 2) {
#ifdef TTK_CONTOUR_FORESTS_ENABLE_PARALLEL_SIMPLIFY
      if (params_->simplifyThreshold) {
         auto  maxSimplifIt = max_element(timeSimplify.cbegin(), timeSimplify.cend());
         float maxSimplif   = *maxSimplifIt;
         cout << "Local simplification maximum time :                         " << maxSimplif;
         cout << " ( " << nbPairMerged << " pairs merged )" << endl;
      }
#endif
      auto maxProcSpeed = max_element(speedProcess.cbegin(), speedProcess.cend());
      auto minProcSpeed = min_element(speedProcess.cbegin(), speedProcess.cend());
      cout << "process speed : ";
      cout << " min is " << *minProcSpeed << " vert/sec";
      cout << " max is " << *maxProcSpeed << " vert/sec";
      cout << endl;
   }

   return 0;
}

//}

#endif /* end of include guard: _CONTOURFORESTSTEMPLATE_H */

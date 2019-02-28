#include <HarmonicFieldComputation.h>

ttk::HarmonicFieldComputation::HarmonicFieldComputation()
    : vertexNumber_{}, constraintNumber_{}, useCotanWeights_{false},
      triangulation_{}, sources_{}, constraints_{},
      outputScalarFieldPointer_{}, solverType_{Cholesky} {
#if defined(TTK_ENABLE_EIGEN) && defined(TTK_ENABLE_OPENMP)
  Eigen::setNbThreads(threadNumber_);
#endif // TTK_ENABLE_OPENMP && TTK_ENABLE_EIGEN
}

ttk::SolverType ttk::HarmonicFieldComputation::findBestSolver() const {
  // TODO
  return ttk::SolverType::Cholesky;
}

// main routine
template <typename scalarFieldType>
int ttk::HarmonicFieldComputation::execute() const {

  using std::cout;
  using std::endl;
  using std::stringstream;

#ifdef TTK_ENABLE_EIGEN
  using SpMat = Eigen::SparseMatrix<scalarFieldType>;
  using SpVec = Eigen::SparseVector<scalarFieldType>;
  using TripletType = Eigen::Triplet<scalarFieldType>;

  Timer t;

  // scalar field constraints vertices
  auto identifiers = static_cast<SimplexId *>(sources_);
  // scalar field: 0 everywhere except on constraint vertices
  auto sf = static_cast<scalarFieldType *>(constraints_);

  {
    stringstream msg;
    msg << "[HarmonicFieldComputation] Beginnning computation" << endl;
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  // get unique constraint vertices
  std::set<SimplexId> identifiersSet;
  for (SimplexId i = 0; i < constraintNumber_; ++i) {
    identifiersSet.insert(identifiers[i]);
  }
  // contains vertices with constraints
  std::vector<SimplexId> identifiersVec(identifiersSet.begin(),
                                        identifiersSet.end());

  // graph laplacian of current mesh
  SpMat lap;
  if (useCotanWeights_) {
    lap = compute_laplacian_with_cotan_weights<SpMat, TripletType,
                                               scalarFieldType>();
  } else {
    lap = compute_laplacian<SpMat, TripletType, scalarFieldType>();
  }

  // constraints vector
  SpVec constraints(vertexNumber_);
  for (size_t i = 0; i < identifiersVec.size(); i++) {
    // put constraint at identifier index
    constraints.coeffRef(identifiersVec[i]) = sf[i];
  }

  // penalty matrix
  SpMat penalty(vertexNumber_, vertexNumber_);
  const scalarFieldType alpha = 1.0e8;
  std::vector<TripletType> triplets;
  triplets.reserve(identifiersVec.size());
  for (auto i : identifiersVec) {
    triplets.emplace_back(TripletType(i, i, alpha));
  }
  penalty.setFromTriplets(triplets.begin(), triplets.end());

  SpMat sol;
  int res;

  SolverType sl = solverType_;
  if (solverType_ == Auto) {
    sl = findBestSolver();
  }

  switch (sl) {
  case Cholesky:
    res = solve<SpMat, SpVec, Eigen::SimplicialCholesky<SpMat>>(
        lap, penalty, constraints, sol);
    break;
  case Iterative:
    res = solve<SpMat, SpVec,
                Eigen::ConjugateGradient<SpMat, Eigen::Upper | Eigen::Lower>>(
        lap, penalty, constraints, sol);
    break;
  case Auto:
    res = Eigen::ComputationInfo::InvalidInput;
    break;
  }

  {
    stringstream msg;
    auto info = static_cast<Eigen::ComputationInfo>(res);
    switch (info) {
    case Eigen::ComputationInfo::Success:
      msg << "[HarmonicFieldComputation] Success!" << endl;
      break;
    case Eigen::ComputationInfo::NumericalIssue:
      msg << "[HarmonicFieldComputation] Numerical Issue!" << endl;
      break;
    case Eigen::ComputationInfo::NoConvergence:
      msg << "[HarmonicFieldComputation] No Convergence!" << endl;
      break;
    case Eigen::ComputationInfo::InvalidInput:
      msg << "[HarmonicFieldComputation] Invalid Input!" << endl;
      break;
    }
    dMsg(cout, msg.str(), advancedInfoMsg);
  }

  auto outputScalarField =
      static_cast<scalarFieldType *>(outputScalarFieldPointer_);
  for (SimplexId i = 0; i < vertexNumber_; ++i) {
    // TODO avoid copy here
    outputScalarField[i] = std::abs(sol.coeffRef(i, 0));
  }

  {
    stringstream msg;
    msg << "[HarmonicFieldComputation] Ending computation after "
        << t.getElapsedTime() << "s (" << threadNumber_ << " thread(s))"
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

#else
  {
    stringstream msg;
    msg << "[HarmonicFieldComputation] Eigen support disabled, computation "
           "skipped "
        << endl;
    dMsg(cout, msg.str(), infoMsg);
  }
#endif // TTK_ENABLE_EIGEN

  return 0;
}

// explicit instantiations for double and float types
template int ttk::HarmonicFieldComputation::execute<double>() const;
template int ttk::HarmonicFieldComputation::execute<float>() const;

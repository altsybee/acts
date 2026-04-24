#include "ExampleAlgorithmBindings.hpp"
#include <algorithm>

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(AliceActsPythonBindings, mex)
{
  mex.doc() = "Alice Acts Examples";

  // Import acts.examples to register IAlgorithm base class
  py::module_::import("acts.examples");

  // Call all algorithm binding functions here
  AliceActsPython::addExampleAlgorithm(mex);
  
  AliceActsPython::addTrackTruthMatcher(mex);

  AliceActsPython::addHitRemoverAlgorithm(mex);

  AliceActsPython::addTrackMergerAlgorithm(mex);

  AliceActsPython::addLutMakerAlgorithm(mex);

  // Validation
  //AliceActsPython::addResPlotTool(mex);
  //AliceActsPython::addRootTrackFitterPerformanceWriter(mex);
  // AliceActsPython::addRootTrackFinderPerformanceWriter(mex);
}

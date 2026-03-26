#pragma once

#include "ExampleAlgorithm.hpp"
#include "TrackTruthMatcher.hpp"
#include "HitRemoverAlgorithm.hpp"
#include "TrackMergerAlgorithm.hpp"

#include "ResPlotTool.hpp"
#include "RootTrackFitterPerformanceWriter.hpp"

#include <algorithm>
#include "ActsPython/Utilities/Macros.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace AliceActsPython
{

// Inline ensures the symbol is available in the .so
inline void addExampleAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::ExampleAlgorithm, mex, "ExampleAlgorithm", testString);
}

inline void addTrackTruthMatcher(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::TrackTruthMatcher, mex, "TrackTruthMatcher",
    inputTracks, inputParticles,
    inputMeasurementParticlesMap, outputTrackParticleMatching,
    outputParticleTrackMatching, matchingRatio, doubleMatching,
    looperProtection, loop_absEta, loop_maxPt, loop_maxParticleHits);
}

inline void addHitRemoverAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::HitRemoverAlgorithm, mex, "HitRemoverAlgorithm",
    inputMeasurements, inputTracks, inputMeasurementParticlesMap,
    sortByOldIndex,
    usedIndices, outputMeasurements,
    outputMeasurementParticlesMap, outputParticleMeasurementsMap,
    outputIndexingMap);
}

inline void addTrackMergerAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::TrackMergerAlgorithm, mex, "TrackMergerAlgorithm",
    inputTrackCollections, inputIndexingMaps,
    outputTrackCollection);
}

  //inline void addResPlotTool(pybind11::module& mex)
  //{
  //py::class_<AliceActsTrk::ResPlotTool::Config>(mex, "ResPlotToolConfig")
  //  .def(py::init<>())
  //  .def_readwrite("varBinning", &AliceActsTrk::ResPlotTool::Config::varBinning);
  //}

  //inline void addRootTrackFitterPerformanceWriter(pybind11::module& mex)
  //{
  //ACTS_PYTHON_DECLARE_WRITER(AliceActsTrk::RootTrackFitterPerformanceWriter, mex,
  //                           "RootTrackFitterPerformanceWriter", inputTracks,
  //                           inputParticles, inputTrackParticleMatching,
  //                           filePath, resPlotToolConfig, effPlotToolConfig,
  //                           trackSummaryPlotToolConfig);
  //}

// inline void addRootTrackFinderPerformanceWriter(pybind11::module& mex) {
//   ACTS_PYTHON_DECLARE_WRITER(AliceActsTrk::RootTrackFinderPerformanceWriter, mex,
//                              "RootTrackFinderPerformanceWriter", inputTracks,
//                              inputParticles, inputTrackParticleMatching,
//                              inputParticleTrackMatching, inputParticleMeasurementsMap,
//                              filePath, fileMode, effPlotToolConfig, fakePlotToolConfig,
//                              duplicationPlotToolConfig, trackSummaryPlotToolConfig,
//                              trackQualityPlotToolConfig, subDetectorTrackSummaryVolumes,
//                              writeMatchingDetails);
// }

} // namespace AliceActsPython

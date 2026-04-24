#pragma once

#include "ExampleAlgorithm.hpp"
#include "TrackTruthMatcher.hpp"
#include "HitRemoverAlgorithm.hpp"
#include "TrackMergerAlgorithm.hpp"
#include "LutMaker.hpp"

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

inline void addLutMakerAlgorithm(pybind11::module& mex)
{
  ACTS_PYTHON_DECLARE_ALGORITHM(
    AliceActsTrk::LutMaker, mex, "LutMakerAlgorithm",
    inputTracks,
    inputParticles,
    // inputMeasurementParticlesMap, 
    inputTrackParticleMatching,
    inputParticleTrackMatching
    // , matchingRatio, doubleMatching,
    // looperProtection, loop_absEta, loop_maxPt, loop_maxParticleHits
  );
}

} // namespace AliceActsPython

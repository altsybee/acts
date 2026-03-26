#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include <memory>
#include <vector>
#include <string>

using namespace ActsExamples;

namespace AliceActsTrk
{

class TrackMergerAlgorithm final : public IAlgorithm
{

 public:
  struct Config {
    std::vector<std::string> inputTrackCollections{};

    /// This will have size = inputTrackCollections.size()-1,
    /// because we assume that the first track collection
    /// used the original measurement vector

    std::vector<std::string> inputIndexingMaps{};
    std::string outputTrackCollection = "";
  };

  TrackMergerAlgorithm(Config cfg, std::unique_ptr<const Acts::Logger> _logger = nullptr);
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::vector<ReadDataHandle<ConstTrackContainer>> m_inputTrackCollections;
  std::vector<ReadDataHandle<std::vector<size_t>>> m_inputIndexingMaps;
  WriteDataHandle<ConstTrackContainer> m_outputTrackCollection{this, "outputTracks"};

}; // TrackMergerAlgorithm
} // namespace AliceActsTrk

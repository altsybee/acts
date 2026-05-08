#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include <memory>
#include <vector>
#include <string>

namespace AliceActsTrk
{

class TrackMergerAlgorithm final : public ActsExamples::IAlgorithm
{

 public:
  struct Config {
    std::vector<std::string> inputTrackCollections{};
    std::string outputTrackCollection = "";
  };

  TrackMergerAlgorithm(Config cfg, std::unique_ptr<const Acts::Logger> _logger = nullptr);
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::vector<ActsExamples::ReadDataHandle<ActsExamples::ConstTrackContainer>>m_inputTrackCollections;
  ActsExamples::WriteDataHandle<ActsExamples::ConstTrackContainer>m_outputTrackCollection{this, "outputTracks"};

}; // TrackMergerAlgorithm
} // namespace AliceActsTrk

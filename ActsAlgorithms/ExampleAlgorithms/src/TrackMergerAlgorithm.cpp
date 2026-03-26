#include "TrackMergerAlgorithm.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/TrackStateProxy.hpp"

#include "ActsExamples/EventData/IndexSourceLink.hpp"

using namespace ActsExamples;

namespace AliceActsTrk
{

TrackMergerAlgorithm::TrackMergerAlgorithm(TrackMergerAlgorithm::Config cfg,
                                           std::unique_ptr<const Acts::Logger> _logger)
  : IAlgorithm("TrackMergerAlgorithm", std::move(_logger)), m_cfg(std::move(cfg))
{

  m_inputTrackCollections.reserve(m_cfg.inputTrackCollections.size());

  for (size_t itc = 0; itc < m_cfg.inputTrackCollections.size(); itc++) {

    auto tc = m_cfg.inputTrackCollections[itc];

    ACTS_DEBUG("Adding trackCollection " << tc << " to the merged output");

    m_inputTrackCollections.emplace_back(this, tc);
    m_inputTrackCollections[itc].initialize(tc);
  }

  m_inputIndexingMaps.reserve(m_cfg.inputIndexingMaps.size());

  for (size_t imap = 0; imap < m_cfg.inputIndexingMaps.size(); imap++) {

    auto map = m_cfg.inputIndexingMaps[imap];
    ACTS_DEBUG("Loading measurement mapping " << map);
    m_inputIndexingMaps.emplace_back(this, map);
    m_inputIndexingMaps[imap].initialize(map);
  }

  m_outputTrackCollection.initialize(m_cfg.outputTrackCollection);
}

ProcessCode TrackMergerAlgorithm::execute(const AlgorithmContext& ctx) const
{

  // Read input data
  std::vector<ConstTrackContainer> trackCollections;

  for (size_t itc = 0; itc < m_cfg.inputTrackCollections.size(); itc++) {
    trackCollections.push_back(m_inputTrackCollections[itc](ctx));
  }

  ACTS_DEBUG("Retrieved and merging " << trackCollections.size() << " track collections");

  // Read input index mapping
  std::vector<std::vector<size_t>> indexingMaps;

  for (size_t imap = 0; imap < m_cfg.inputIndexingMaps.size(); imap++) {
    indexingMaps.push_back(m_inputIndexingMaps[imap](ctx));
  }

  ACTS_DEBUG("Retrieved " << indexingMaps.size() << " indexing maps");

  // Mutable Output Collection
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer actsTracksContainer(trackContainer, trackStateContainer);

  // Make sure the containers have the same dynamic content. I assume the track Collection size is not empty.
  // And all tracks have the same dynamic columns

  actsTracksContainer.ensureDynamicColumns(trackCollections[0]);

  for (size_t itc = 0; itc < trackCollections.size(); itc++) {

    ACTS_DEBUG("Checking track collection index " << itc);

    // Find out what type of TrackStates this collection has.
    // Only use the first track.
    size_t track_counter = 0;
    bool filledTrackStates = true;
    for (const auto& track : trackCollections[itc]) {

      // Check if the track states are filled or not
      if (track_counter == 0) {
        for (const auto& ts : track.trackStatesReversed()) {
          if (ts.getMask() == Acts::TrackStatePropMask::None) {
            filledTrackStates = false;
            break;
          }
        }
        track_counter++;
      }

      auto actsDestProxy = actsTracksContainer.makeTrack();

      if (filledTrackStates) {
        ACTS_DEBUG("Track has filled Track States from collection " << itc);
        actsDestProxy.copyFrom(track);

        // tracks.addColumn<BranchStopper::BranchState>("MyBranchState");
        // tracksTemp.addColumn<BranchStopper::BranchState>("MyBranchState");

        // tracks.addColumn<unsigned int>("trackGroup");
        // tracksTemp.addColumn<unsigned int>("trackGroup");

      } else {
        ACTS_DEBUG("Copying track");
        // actsDestProxy.copyFromShallow(track);
        actsDestProxy.copyFromWithoutStates(track);
        // append track states (cheap), but they're in the wrong order
        for (const auto& srcTrackState : track.trackStatesReversed()) {

          auto destTrackState = actsDestProxy.appendTrackState(srcTrackState.getMask());
          destTrackState.copyFrom(srcTrackState, Acts::TrackStatePropMask::None,
                                  true);

          // Check the track State uncalibrated source link index
          IndexSourceLink sl = destTrackState.getUncalibratedSourceLink().template get<IndexSourceLink>();
          auto hitIndex = sl.index();

          // For collections after the first get the proper indexing
          if (itc > 0) {
            auto og_index = indexingMaps[itc - 1][hitIndex];
            IndexSourceLink remapped_sl(sl.geometryId(), og_index);
            destTrackState.setUncalibratedSourceLink(Acts::SourceLink{remapped_sl});
          }

          // Check!
          IndexSourceLink sl2 = destTrackState.getUncalibratedSourceLink().template get<IndexSourceLink>();
          auto hitIndex2 = sl2.index();

          actsDestProxy.reverseTrackStates();
        }
      }

      // detail::ExpectedLayerPatternHelper::set(actsDestProxy, expectedLayerPattern);
    }
  }

  ACTS_DEBUG("Merged outputTrackContainer size " << actsTracksContainer.size());

  auto constTrackStateContainer =
    std::make_shared<Acts::ConstVectorMultiTrajectory>(
      std::move(*trackStateContainer));

  auto constTrackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(
    std::move(*trackContainer));

  ConstTrackContainer constTracks{constTrackContainer,
                                  constTrackStateContainer};

  m_outputTrackCollection(ctx, std::move(constTracks));

  return ProcessCode::SUCCESS;
}

} // namespace AliceActsTrk

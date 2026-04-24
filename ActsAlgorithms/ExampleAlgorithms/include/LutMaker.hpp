// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "FlatLutEntry.h"

#include <map>
#include <string>
#include <array>
#include <cstdint>
#include <vector>
class TH2F;

#if 1
#define LOGF(level, ...) printf(__VA_ARGS__);

namespace framework
{
class runtime_error_f : public std::runtime_error
{
 public:
  template <typename... Args>
  runtime_error_f(const char* format, Args... args) : std::runtime_error(format)
  {
    size_t size = std::snprintf(nullptr, 0, format, args...) + 1; // Extra space for null terminator
    std::vector<char> buffer(size);
    std::snprintf(buffer.data(), size, format, args...);
    this->message = std::string(buffer.data(), buffer.size() - 1); // Exclude null terminator
  }
  const char* what() const noexcept override
  {
    return message.c_str();
  }

 private:
  std::string message;
};

} // namespace framework

#include "FlatLutEntry.h"

#include <cstring>

namespace o2::delphes
{

float map_t::fracPositionWithinBin(float val) const
{
  float width = (max - min) / nbins;
  int bin;
  float returnVal = 0.5f;
  if (log) {
    bin = static_cast<int>((std::log10(val) - min) / width);
    returnVal = ((std::log10(val) - min) / width) - bin;
  } else {
    bin = static_cast<int>((val - min) / width);
    returnVal = val / width - bin;
  }
  return returnVal;
}

int map_t::find(float val) const
{
  float width = (max - min) / nbins;
  int bin;
  if (log) {
    bin = static_cast<int>((std::log10(val) - min) / width);
  } else {
    bin = static_cast<int>((val - min) / width);
  }
  if (bin < 0) {
    return 0;
  }
  if (bin > nbins - 1) {
    return nbins - 1;
  }
  return bin;
}

void map_t::print() const
{
  LOGF(info, "nbins = %d, min = %f, max = %f, log = %s \n", nbins, min, max, log ? "on" : "off");
}

bool lutHeader_t::check_version() const
{
  return (version == LUTCOVM_VERSION);
}

void lutHeader_t::print() const
{
  LOGF(info, " version: %d \n", version);
  LOGF(info, "     pdg: %d \n", pdg);
  LOGF(info, "   field: %f \n", field);
  LOGF(info, "  nchmap: ");
  nchmap.print();
  LOGF(info, "  radmap: ");
  radmap.print();
  LOGF(info, "  etamap: ");
  etamap.print();
  LOGF(info, "   ptmap: ");
  ptmap.print();
}

void FlatLutData::initialize(const lutHeader_t& header)
{
  mNchBins = header.nchmap.nbins;
  mRadBins = header.radmap.nbins;
  mEtaBins = header.etamap.nbins;
  mPtBins = header.ptmap.nbins;

  size_t headerSize = sizeof(lutHeader_t);
  size_t numEntries = static_cast<size_t>(mNchBins) * mRadBins * mEtaBins * mPtBins;
  size_t entriesSize = numEntries * sizeof(lutEntry_t);
  size_t totalSize = headerSize + entriesSize;

  mData.resize(totalSize);
  // Write header at the beginning
  std::memcpy(mData.data(), &header, headerSize);
  updateRef();
}

size_t FlatLutData::getEntryOffset(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const
{
  size_t headerSize = sizeof(lutHeader_t);

  // Linear index: nch varies slowest, pt varies fastest
  // idx = nch * (rad*eta*pt) + rad * (eta*pt) + eta * pt + pt
  size_t linearIdx = static_cast<size_t>(nch_bin) * (mRadBins * mEtaBins * mPtBins) + static_cast<size_t>(rad_bin) * (mEtaBins * mPtBins) + static_cast<size_t>(eta_bin) * mPtBins + static_cast<size_t>(pt_bin);

  return headerSize + linearIdx * sizeof(lutEntry_t);
}

const lutEntry_t* FlatLutData::getEntryRef(int nch_bin, int rad_bin, int eta_bin, int pt_bin) const
{
  size_t offset = getEntryOffset(nch_bin, rad_bin, eta_bin, pt_bin);
  return reinterpret_cast<const lutEntry_t*>(mDataRef.data() + offset);
}

lutEntry_t* FlatLutData::getEntry(int nch_bin, int rad_bin, int eta_bin, int pt_bin)
{
  size_t offset = getEntryOffset(nch_bin, rad_bin, eta_bin, pt_bin);
  return reinterpret_cast<lutEntry_t*>(mData.data() + offset);
}

const lutHeader_t& FlatLutData::getHeaderRef() const
{
  return *reinterpret_cast<const lutHeader_t*>(mDataRef.data());
}

lutHeader_t& FlatLutData::getHeader()
{
  return *reinterpret_cast<lutHeader_t*>(mData.data());
}

void FlatLutData::updateRef()
{
  mDataRef = std::span{mData.data(), mData.size()};
}

void FlatLutData::cacheDimensions()
{
  auto const& header = getHeaderRef();
  mNchBins = header.nchmap.nbins;
  mRadBins = header.radmap.nbins;
  mEtaBins = header.etamap.nbins;
  mPtBins = header.ptmap.nbins;
}

void FlatLutData::resetDimensions()
{
  mNchBins = 0;
  mRadBins = 0;
  mEtaBins = 0;
  mPtBins = 0;
}

void FlatLutData::adopt(const uint8_t* buffer, size_t size)
{
  mData.resize(size);
  std::memcpy(mData.data(), buffer, size);
  updateRef();
  cacheDimensions();
}

void FlatLutData::view(const uint8_t* buffer, size_t size)
{
  mData.clear();
  mDataRef = std::span{buffer, size};
  cacheDimensions();
}

void FlatLutData::validateBuffer(const uint8_t* buffer, size_t size)
{
  auto header = PreviewHeader(buffer, size);
  auto mNchBins = header.nchmap.nbins;
  auto mRadBins = header.radmap.nbins;
  auto mEtaBins = header.etamap.nbins;
  auto mPtBins = header.ptmap.nbins;

  size_t expectedSize = sizeof(lutHeader_t) + static_cast<size_t>(mNchBins) * mRadBins * mEtaBins * mPtBins * sizeof(lutEntry_t);

  if (size < expectedSize) {
    throw framework::runtime_error_f("Buffer size mismatch: expected %zu, got %zu", expectedSize, size);
  }
}

lutHeader_t FlatLutData::PreviewHeader(const uint8_t* buffer, size_t size)
{
  if (size < sizeof(lutHeader_t)) {
    throw framework::runtime_error_f("Buffer too small for LUT header: expected at least %zu, got %zu", sizeof(lutHeader_t), size);
  }
  const auto* header = reinterpret_cast<const lutHeader_t*>(buffer);
  if (!header->check_version()) {
    throw framework::runtime_error_f("LUT header version mismatch: expected %d, got %d", LUTCOVM_VERSION, header->version);
  }
  return *header;
}

FlatLutData FlatLutData::AdoptFromBuffer(const uint8_t* buffer, size_t size)
{
  validateBuffer(buffer, size);
  FlatLutData data;

  // Copy buffer
  data.adopt(buffer, size);
  return data;
}

FlatLutData FlatLutData::ViewFromBuffer(const uint8_t* buffer, size_t size)
{
  validateBuffer(buffer, size);
  FlatLutData data;

  // Store reference to external buffer
  // WARNING: Caller must ensure buffer lifetime exceeds FlatLutData usage
  data.view(buffer, size);
  return data;
}

FlatLutData FlatLutData::ViewFromBuffer(std::span<std::byte> const& span)
{
  return ViewFromBuffer(reinterpret_cast<const uint8_t*>(span.data()), span.size_bytes());
}

bool FlatLutData::isLoaded() const
{
  return ((!mData.empty()) || (!mDataRef.empty()));
}

lutHeader_t FlatLutData::PreviewHeader(std::ifstream& file, const char* filename)
{
  lutHeader_t tempHeader;
  file.read(reinterpret_cast<char*>(&tempHeader), sizeof(lutHeader_t));
  if (file.gcount() != static_cast<std::streamsize>(sizeof(lutHeader_t))) {
    throw framework::runtime_error_f("Failed to read LUT header from %s", filename);
  }
  if (!tempHeader.check_version()) {
    throw framework::runtime_error_f("LUT header version mismatch: expected %d, got %d", LUTCOVM_VERSION, tempHeader.version);
  }
  return tempHeader;
}

FlatLutData FlatLutData::loadFromFile(std::ifstream& file, const char* filename)
{
  // Read header first
  lutHeader_t tempHeader = PreviewHeader(file, filename);

  FlatLutData data;

  // Initialize flat data structure
  data.initialize(tempHeader);

  // Read all entries sequentially into flat buffer
  size_t headerSize = sizeof(lutHeader_t);
  size_t numEntries = static_cast<size_t>(data.mNchBins) * data.mRadBins * data.mEtaBins * data.mPtBins;
  size_t entriesSize = numEntries * sizeof(lutEntry_t);

  file.read(reinterpret_cast<char*>(data.data() + headerSize), entriesSize);
  if (file.gcount() != static_cast<std::streamsize>(entriesSize)) {
    throw framework::runtime_error_f("Failed to read LUT entries from %s: expected %zu bytes, got %zu", filename, entriesSize, static_cast<size_t>(file.gcount()));
  }

  LOGF(info, "Successfully loaded LUT from %s: %zu entries", filename, numEntries);
  return data;
}

void FlatLutData::reset()
{
  mData.clear();
  updateRef();
  resetDimensions();
}

} // namespace o2::delphes

#endif

using namespace ActsExamples;

namespace AliceActsTrk
{

struct LutContainer {
  o2::delphes::lutHeader_t header;
  std::vector<o2::delphes::lutEntry_t> entries;
};

/// Matches tracks to truth particles and vice versa
class LutMaker final : public IAlgorithm
{
 public:
  struct Config {
    std::string inputTracks;                /// Input (fitted) tracks collection
    std::string inputParticles;             /// Input particles collection.
    std::string inputTrackParticleMatching; /// Input track-particle matching.
    std::string inputParticleTrackMatching; /// Input particle-track matching.
    std::string lutTag;                     /// Tag to identify the LUT (e.g. "lut_v1")
    float magneticField;                    /// Magnetic field value to store in LUT header (in Tesla)
  };

  LutMaker(const Config& config, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ~LutMaker();

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<TrackParametersContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{this, "InputMeasurementParticlesMap"};
  ReadDataHandle<TrackParticleMatching> m_inputTrackParticleMatching{this, "InputTrackParticleMatching"};
  ReadDataHandle<ParticleTrackMatching> m_inputParticleTrackMatching{this, "InputParticleTrackMatching"};

  static std::map<int, LutContainer> mLutMap;                         ///< Map of LUTs indexed by PDG code
  static std::map<int, std::vector<std::uint64_t>> mPtTruthCounts;    ///< Truth particle counts per PDG and pT bin
  static std::map<int, std::vector<std::uint64_t>> mPtMatchedCounts;  ///< Matched track counts per PDG and pT bin
  static std::map<int, std::vector<std::array<double, 15>>> mCovSums; ///< Running sum of covariance lower-triangle elements per PDG and entry index
  static std::map<int, std::vector<std::uint64_t>> mEntryCounts;      ///< Number of accumulated tracks per PDG and entry index
  static std::map<int, std::vector<double>> mQopTResSums;             ///< Running sum of q/pT residuals per PDG and entry index
  static std::map<int, std::vector<double>> mQopTRes2Sums;            ///< Running sum of squared q/pT residuals per PDG and entry index

  static std::map<int, TH2F*> qopt_vs_pT; ///< Resolution graphs for q/pT vs pT per PDG
};

} // namespace AliceActsTrk

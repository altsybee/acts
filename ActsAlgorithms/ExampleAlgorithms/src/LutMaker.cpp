// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "LutMaker.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsExamples/Validation/TrackClassification.hpp"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace ActsExamples;

namespace AliceActsTrk
{

std::map<int, LutContainer> LutMaker::mLutMap;
std::map<int, std::vector<std::uint64_t>> LutMaker::mPtTruthCounts;
std::map<int, std::vector<std::uint64_t>> LutMaker::mPtMatchedCounts;
std::map<int, std::vector<std::array<double, 15>>> LutMaker::mCovSums;
std::map<int, std::vector<std::uint64_t>> LutMaker::mEntryCounts;
std::map<int, std::vector<double>> LutMaker::mQopTResSums;
std::map<int, std::vector<double>> LutMaker::mQopTRes2Sums;
std::map<int, TH2F*> LutMaker::qopt_vs_pT;

LutMaker::LutMaker(const Config& config,
                   std::unique_ptr<const Acts::Logger> logger) : IAlgorithm("LutMaker", std::move(logger)),
                                                                 m_cfg(config)
{
  if (m_cfg.inputTracks.empty()) {
    throw std::invalid_argument("Missing input tracks");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles");
  }

  std::vector<int> pdgCodes = {11, 13, 211, 321, 2212}; // example PDG codes for electron, muon, pion, kaon, proton
  std::map<int, float> pdgMasses = {
    {11, 0.000511f}, // electron mass in GeV/c^2
    {13, 0.10566f},  // muon mass in GeV/c^2
    {211, 0.13957f}, // pion mass in GeV/c^2
    {321, 0.49367f}, // kaon mass in GeV/c^2
    {2212, 0.93827f} // proton mass in GeV/c^2
  };
  for (int pdg : pdgCodes) {
    printf("Initializing LUT for PDG %d\n", pdg);
    o2::delphes::lutHeader_t header;
    header.version = LUTCOVM_VERSION;
    header.pdg = pdg;
    header.mass = pdgMasses.at(pdg);
    header.field = 0.5f; // example magnetic field value in Tesla
    header.nchmap = {1, 0.f, 1000.f, false};
    header.radmap = {1, 0.f, 100.0f, false};
    header.etamap = {10, -4.0f, 4.0f, false};
    // Log10 pT binning from 1e-2 to 1e2 GeV/c.
    header.ptmap = {200, -2.0f, 2.0f, true};
    mLutMap[pdg].header = header;

    printf("Initalizing lutmap\n");
    o2::delphes::lutEntry_t entry;
    // Now filling with empty LUT entries
    for (int nch_bin = 0; nch_bin < header.nchmap.nbins; ++nch_bin) {
      for (int rad_bin = 0; rad_bin < header.radmap.nbins; ++rad_bin) {
        for (int eta_bin = 0; eta_bin < header.etamap.nbins; ++eta_bin) {
          for (int pt_bin = 0; pt_bin < header.ptmap.nbins; ++pt_bin) {
            entry.nch = header.nchmap.eval(nch_bin);
            entry.eta = header.etamap.eval(eta_bin);
            entry.pt = header.ptmap.eval(pt_bin);
            mLutMap[pdg].entries.push_back(entry);
          }
        }
      }
    }

    mPtTruthCounts[pdg] = std::vector<std::uint64_t>(header.ptmap.nbins, 0u);
    mPtMatchedCounts[pdg] = std::vector<std::uint64_t>(header.ptmap.nbins, 0u);

    const std::size_t nEntries = static_cast<std::size_t>(header.nchmap.nbins) *
                                 header.radmap.nbins * header.etamap.nbins * header.ptmap.nbins;
    mCovSums[pdg] = std::vector<std::array<double, 15>>(nEntries, {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.});
    mEntryCounts[pdg] = std::vector<std::uint64_t>(nEntries, 0u);
    mQopTResSums[pdg] = std::vector<double>(nEntries, 0.0);
    mQopTRes2Sums[pdg] = std::vector<double>(nEntries, 0.0);

    qopt_vs_pT[pdg] = new TH2F(Form("qopT_res_vs_pT_pdg%d", pdg),
                               Form("q/pT resolution vs pT for PDG %d; pT [GeV/c]; q/pT residual", pdg),
                               40, 0.f, 100.f, 1000, -0.1f, 0.1f);
  }

  m_inputTracks.initialize(m_cfg.inputTracks);
  m_inputParticles.initialize(m_cfg.inputParticles);
  if (!m_cfg.inputTrackParticleMatching.empty()) {
    m_inputTrackParticleMatching.initialize(m_cfg.inputTrackParticleMatching);
  }
}

LutMaker::~LutMaker()
{
  printf("Called LutMaker destructor, with lutMap of size: %zu\n", mLutMap.size());

  for (auto& [pdg, lut] : mLutMap) {
    // First fill the efficiency values into the LUT
    const auto& truthCounts = mPtTruthCounts[pdg];
    const auto& matchedCounts = mPtMatchedCounts[pdg];

    const int nnch = lut.header.nchmap.nbins;
    const int nrad = lut.header.radmap.nbins;
    const int neta = lut.header.etamap.nbins;
    const int npt = lut.header.ptmap.nbins;

    // Fill efficiency and average covariance + eigen decomposition per entry
    const auto& covSums = mCovSums[pdg];
    const auto& entryCounts = mEntryCounts[pdg];
    const auto& qopTResSums = mQopTResSums[pdg];
    const auto& qopTRes2Sums = mQopTRes2Sums[pdg];

    for (int nch_bin = 0; nch_bin < nnch; ++nch_bin) {
      for (int rad_bin = 0; rad_bin < nrad; ++rad_bin) {
        for (int eta_bin = 0; eta_bin < neta; ++eta_bin) {
          for (int pt_bin = 0; pt_bin < npt; ++pt_bin) {
            const std::size_t entryIdx = static_cast<std::size_t>(nch_bin) * (nrad * neta * npt) + rad_bin * (neta * npt) + eta_bin * npt + pt_bin;
            auto& entry = lut.entries[entryIdx];

            // Efficiency
            const auto denom = truthCounts[pt_bin];
            const auto numer = matchedCounts[pt_bin];
            entry.eff = denom > 0u ? static_cast<float>(numer) / static_cast<float>(denom) : 0.f;

            // Average covariance and eigen decomposition
            const auto count = entryCounts[entryIdx];
            if (count > 0u) {
              Eigen::Matrix<float, 5, 5> cov5 = Eigen::Matrix<float, 5, 5>::Zero();
              int covIdx = 0;
              for (int i = 0; i < 5; ++i) {
                for (int j = 0; j <= i; ++j) {
                  const float avg = static_cast<float>(covSums[entryIdx][covIdx] / static_cast<double>(count));
                  entry.covm[covIdx] = avg;
                  cov5(i, j) = avg;
                  cov5(j, i) = avg;
                  ++covIdx;
                }
              }

              // Calibrate q/pT variance to observed residual width in this bin.
              if (count > 1u) {
                const double meanRes = qopTResSums[entryIdx] / static_cast<double>(count);
                const double varRes = qopTRes2Sums[entryIdx] / static_cast<double>(count) - meanRes * meanRes;
                if (varRes > 0.0) {
                  const float qopTVar = static_cast<float>(varRes);
                  cov5(4, 4) = qopTVar;
                  entry.covm[14] = qopTVar;
                }
              }

              Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float, 5, 5>> eigSolver(cov5);
              if (eigSolver.info() == Eigen::Success) {
                const auto& eigValues = eigSolver.eigenvalues();
                const auto& eigVectors = eigSolver.eigenvectors();
                const auto eigVectorsInv = eigVectors.inverse();
                for (int i = 0; i < 5; ++i) {
                  entry.eigval[i] = eigValues(i);
                  for (int j = 0; j < 5; ++j) {
                    entry.eigvec[i][j] = eigVectors(i, j);
                    entry.eiginv[i][j] = eigVectorsInv(i, j);
                  }
                }
              } else {
                for (int i = 0; i < 5; ++i) {
                  entry.eigval[i] = cov5(i, i);
                  for (int j = 0; j < 5; ++j) {
                    entry.eigvec[i][j] = (i == j) ? 1.f : 0.f;
                    entry.eiginv[i][j] = (i == j) ? 1.f : 0.f;
                  }
                }
              }
              entry.valid = true;
            }
          }
        }
      }
    }

    const std::string lutFileTag = Form("%s%s_%.2fT", m_cfg.lutTag.empty() ? "" : "_", lut.header.pdg, m_cfg.magneticField);
    const std::string filename = std::string("lut_pdg") + std::to_string(pdg) + lutFileTag + ".dat";
    printf("Saving LUT for PDG %d to file %s\n", pdg, filename.c_str());
    // output file
    std::ofstream lutFile(filename, std::ofstream::binary);
    if (!lutFile.is_open()) {
      LOGF(info, "Did not manage to open output file!!");
      return;
    }

    o2::delphes::lutHeader_t lutHeader = lut.header;
    lutFile.write(reinterpret_cast<char*>(&lutHeader), sizeof(o2::delphes::lutHeader_t));

    printf("LUT dimensions: nch_bins = %d, rad_bins = %d, eta_bins = %d, pt_bins = %d\n", nnch, nrad, neta, npt);

    for (int nch_bin = 0; nch_bin < nnch; ++nch_bin) {
      for (int rad_bin = 0; rad_bin < nrad; ++rad_bin) {
        for (int eta_bin = 0; eta_bin < neta; ++eta_bin) {
          for (int pt_bin = 0; pt_bin < npt; ++pt_bin) {
            auto* entry = &lut.entries[nch_bin * (nrad * neta * npt) + rad_bin * (neta * npt) + eta_bin * npt + pt_bin];
            if (entry->valid) {
              entry->print();
            }
            lutFile.write(reinterpret_cast<const char*>(entry), sizeof(o2::delphes::lutEntry_t));
          }
        }
      }
    }

    lutFile.close();
  }

  TFile output("/tmp/lutMakerQA.root", "RECREATE");
  for (const auto& [pdg, hist] : qopt_vs_pT) {
    hist->Write();
    TH1D* reswidth_qopt_vs_pT = hist->ProjectionX(Form("reswidth_qopT_res_vs_pT_pdg%d", pdg));
    TH1D* resmean_qopt_vs_pT = hist->ProjectionX(Form("resmean_qopT_res_vs_pT_pdg%d", pdg));
    TProfile* profileX = hist->ProfileX(Form("qopT_res_vs_pT_pdg%d_profileX", pdg), 1, -1, "s");
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
      reswidth_qopt_vs_pT->SetBinContent(i, profileX->GetBinError(i));
      resmean_qopt_vs_pT->SetBinContent(i, profileX->GetBinContent(i));
    }
    reswidth_qopt_vs_pT->Write();
    resmean_qopt_vs_pT->Write();
  }
}

ActsExamples::ProcessCode LutMaker::execute(const ActsExamples::AlgorithmContext& ctx) const
{
  // Read input tracks
  const auto& tracks = m_inputTracks(ctx);
  printf("Number of input tracks: %zu\n", tracks.size());

  // for (const auto& track : tracks) {
  //   printf("Track pT: %f, eta: %f, phi: %f\n", track.pt(), track.eta(), track.phi());
  // }

  const auto& particles = m_inputParticles(ctx);
  const TrackParticleMatching* trackParticleMatching = nullptr;
  if (!m_cfg.inputTrackParticleMatching.empty()) {
    trackParticleMatching = &m_inputTrackParticleMatching(ctx);
  }

  std::vector<const SimParticle*> particlesByIndex;
  particlesByIndex.reserve(particles.size());
  for (const auto& particle : particles) {
    printf(" - Particle index: %zu, PDG: %d\n", particlesByIndex.size(), particle.pdg());
    particlesByIndex.push_back(&particle);

    const int pdg = std::abs(particle.pdg());
    auto lutIt = mLutMap.find(pdg);
    if (lutIt == mLutMap.end()) {
      continue;
    }

    const auto& p = particle.momentum();
    const double ptTrue = std::hypot(p.x(), p.y());
    const int ptBin = lutIt->second.header.ptmap.find(static_cast<float>(ptTrue));
    ++mPtTruthCounts[pdg][ptBin];
  }

  for (std::size_t itrk = 0; itrk < tracks.size(); ++itrk) {
    printf("Processing track index %zu\n", itrk);
    const auto& track = tracks[itrk];

    const SimParticle* particlePtr = nullptr;
    if (trackParticleMatching != nullptr) {
      auto it = trackParticleMatching->find(itrk);
      if (it == trackParticleMatching->end()) {
        continue;
      }

      if (!it->second.particle.has_value()) {
        continue;
      }
      const auto& particleId = it->second.particle.value();
      auto particleIt = particles.find(particleId);
      if (particleIt == particles.end()) {
        continue;
      }
      particlePtr = &(*particleIt);
    } else {
      if (itrk >= particlesByIndex.size()) {
        continue;
      }
      particlePtr = particlesByIndex[itrk];
    }

    if (particlePtr == nullptr) {
      continue;
    }

    const int pdg = std::abs(particlePtr->pdg());
    if (mLutMap.find(pdg) == mLutMap.end()) {
      printf("Skipping particle with unsupported PDG %d\n", particlePtr->pdg());
      continue;
    } else {
      printf(" + Processing matched particle with PDG %d\n", particlePtr->pdg());
    }

    const double ptReco = track.transverseMomentum();
    const auto& p = particlePtr->momentum();
    const double ptTrue = std::hypot(p.x(), p.y());

    const double relRes = (ptReco - ptTrue) / ptTrue;
    printf("matched pion: ptTrue=%f, ptReco=%f, relRes=%f\n", ptTrue, ptReco, relRes);

    printf("  pxReco = %f, pxTrue = %f\n", track.momentum().x(), p.x());
    printf("  pyReco = %f, pyTrue = %f\n", track.momentum().y(), p.y());
    printf("  pzReco = %f, pzTrue = %f\n", track.momentum().z(), p.z());
    const float mod = std::sqrt(p.x() * p.x() + p.y() * p.y() + p.z() * p.z());
    const float eta = 0.5f * std::log((mod + p.z()) / (mod - p.z()));

    // Fill the LUT information
    auto& lut = mLutMap[pdg];
    const int nch_bin = 0;
    const int rad_bin = 0; // example with only one radial bin
    const int eta_bin = lut.header.etamap.find(eta);
    const int pt_bin = lut.header.ptmap.find(ptTrue);

    ++mPtMatchedCounts[pdg][pt_bin];

    const std::size_t entryIdx = static_cast<std::size_t>(nch_bin) * (lut.header.radmap.nbins * lut.header.etamap.nbins * lut.header.ptmap.nbins) + rad_bin * (lut.header.etamap.nbins * lut.header.ptmap.nbins) + eta_bin * lut.header.ptmap.nbins + pt_bin;

    // Residual definition aligned with q/pT resolution plots.
    const auto& mom = track.momentum();
    const float pAbs = std::sqrt(mom.x() * mom.x() + mom.y() * mom.y() + mom.z() * mom.z());
    const float ptAbs = std::sqrt(mom.x() * mom.x() + mom.y() * mom.y());
    const float sinTheta = (pAbs > 0.f) ? (ptAbs / pAbs) : 0.f;

    const int pdgSigned = particlePtr->pdg();
    const int pdgAbs = std::abs(pdgSigned);
    int baseCharge = 0;
    if (pdgAbs == 11 || pdgAbs == 13) {
      baseCharge = -1;
    } else if (pdgAbs == 211 || pdgAbs == 321 || pdgAbs == 2212) {
      baseCharge = 1;
    }
    const int charge = (pdgSigned >= 0) ? baseCharge : -baseCharge;

    constexpr float kMinSinTheta = 1.e-6f;
    if (charge != 0 && ptTrue > 0.0 && std::abs(sinTheta) > kMinSinTheta) {
      const double qOverPtReco = static_cast<double>(track.qOverP()) / static_cast<double>(sinTheta);
      const double qOverPtTrue = static_cast<double>(charge) / ptTrue;
      const double qOverPtResidual = qOverPtReco - qOverPtTrue;
      mQopTResSums[pdg][entryIdx] += qOverPtResidual;
      mQopTRes2Sums[pdg][entryIdx] += qOverPtResidual * qOverPtResidual;
      qopt_vs_pT[pdg]->Fill(ptTrue, qOverPtResidual);
    }

    // Accumulate covariance lower-triangle elements for later averaging
    if (track.covariance().has_value()) {
      const auto& trackCov = track.covariance().value();

      // ACTS bound covariance uses (loc0, loc1, phi, theta, q/p, time).
      // Convert to a transverse-momentum basis by replacing q/p -> q/pT.
      Eigen::Matrix<float, 5, 5> covAct = Eigen::Matrix<float, 5, 5>::Zero();
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
          covAct(i, j) = static_cast<float>(trackCov(i, j));
        }
      }

      Eigen::Matrix<float, 5, 5> jac = Eigen::Matrix<float, 5, 5>::Identity();
      if (std::abs(sinTheta) > kMinSinTheta) {
        const float cosThetaSafe = std::clamp(static_cast<float>(mom.z() / pAbs), -1.f, 1.f);
        const float theta = std::acos(cosThetaSafe);
        const float cosTheta = std::cos(theta);
        const float qOverP = track.qOverP();

        // y4 = q/pT = (q/p) / sin(theta)
        jac(4, 4) = 1.f / sinTheta;
        jac(4, 3) = -qOverP * cosTheta / (sinTheta * sinTheta);
      }

      const Eigen::Matrix<float, 5, 5> covTrk = jac * covAct * jac.transpose();
      auto& covSum = mCovSums[pdg][entryIdx];

      int covIdx = 0;
      for (int i = 0; i < 5; ++i) {
        for (int j = 0; j <= i; ++j) {
          covSum[covIdx] += static_cast<double>(covTrk(i, j));
          ++covIdx;
        }
      }
      ++mEntryCounts[pdg][entryIdx];
    }
  }

  return ProcessCode::SUCCESS;
}

} // namespace AliceActsTrk

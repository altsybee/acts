// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Kernel/detail/SimulationActor.hpp"

#include <string>

namespace ActsFatras::detail {

void IA_hit_checker(const Acts::GeometryIdentifier& geoId,
                    const Particle& before, const Particle& after,
                    bool& isInKineCuts, float& probToDetect, bool& isInROF,
                    bool& isHitTOF) {
  // if (1 )
  // std::cout << ">> we are in IA_hit_checker() - test" << std::endl;

  float eta =
      Acts::clampValue<float>(Acts::VectorHelpers::eta(after.direction()));
  float pt = after.transverseMomentum();

  if (0 && fabs(eta) < 0.5) {
    std::uint32_t vertPrimary = after.particleId().vertexPrimary();
    std::cout << "vertPrimary = " << vertPrimary << ", afterBC = "
              << after.BC()
              // << ", beforeBC = " << before.BC()
              << std::endl;
  }

  // float etaCut = 0.9;//5;  // 4.5;//1;
  float etaCut = 4.2;  // 0.9;//5;  // 4.5;//1;
  float ptCutMin = 0.01;
  float ptCutMax = 100;
  isInKineCuts = true;

  if (1) {
    isInKineCuts = fabs(eta) < etaCut && (pt > ptCutMin && pt < ptCutMax);

    // Nov 2025 - check low-energy electrons
    if (0 && isInKineCuts) {
      const Acts::Vector4 mom4 = before.fourMomentum();
      float e = mom4[Acts::eEnergy];
      float px = mom4[Acts::eMom0];
      float py = mom4[Acts::eMom1];
      float pz = mom4[Acts::eMom2];
      float m = sqrt(e * e - px * px - py * py - pz * pz);

      if (fabs(m - 0.000511) < 0.00001 && e < 0.002)
        std::cout 
        << " >> e = " << e 
        << " >> px = " << px
        << " >> py = " << py
        << " >> pz = " << pz
        << " >> m = " << m << std::endl;

      // suppress hits from low-energy electrons
      if ( fabs(m - 0.000511) < 0.00001 && e < 0.002)
        isInKineCuts = false;
    }
  }

  probToDetect = 1;
  // 1 ROF
  // if (1) {
  //   probToDetect = (after.BC() >= 0 && after.BC() < 20)
  //                      ? arr_hit_prob_vs_bc[after.BC()]
  //                      : 1.0;
  // }
  // // 2 ROFs
  // else if (0) {
  //   probToDetect = (after.BC() >= 0 && after.BC() < 40)
  //                      ? arr_hit_prob_vs_bc[after.BC()]
  //                      : 1.0;
  // }
  // // for multiROF:
  // else if (0) {
  //   if (after.BC() < 20)
  //     probToDetect = 0.1;  // "noise from the past"
  //   else
  //     probToDetect = 1.0;  // this and next event
  // }
  // vertPrimary == 1 ? 1 : ((double)rand() / (RAND_MAX)) < 0.2;

  // ### IA: layer-by-layer check - to emulate short vs long ROFs
  isInROF = true;
  if (1 && geoId.sensitive() > 0) {
    // check if we are in IRIS
    // if (geoId.volume() == 12 || geoId.volume() == 13 || geoId.volume() == 14)
    // { // IB, old (before Feb 2025?) if (geoId.volume() == 16 ||
    // geoId.volume() == 17 || geoId.volume() == 18) { // IB, Sept-Oct 2025
    if (geoId.volume() == 22) {  // ML
      // check if we are in short ROF for inner layers
      // isInROF = after.BC() < 20 ? true : false;
      // isInROF = after.BC() < 5 ? true : false;
      isInROF = after.BC() < 34 ? true : false;  // 500 ns ROFs, max nColl
      // isInROF = after.BC() < 50 ? true : false; // 1 us ROFs, max nColl
      // isInROF = after.BC() >= 50 ? true : false; // second 1 us ROFs // Sept
      // 26, 2025

      if (0 && fabs(eta) < 0.5)
        std::cout << geoId.sensitive() << " >> " << geoId.volume() << " "
                  << geoId.layer() << " isInROF = " << isInROF
                  << " BC =  " << after.BC() << ".   " << std::endl;
    }

    if (0)
      std::cout << "    xyz = " << before.fourPosition().x() << " "
                << before.fourPosition().y() << " " << before.fourPosition().z()
                << std::endl;
  }

  // MAY 2025: reduced prob to have hits in ML due to occupancy
  if (0 && geoId.volume() == 22)  // ML
  {
    if (geoId.layer() == 2)
      // probToDetect = 1 - 0.015;
      probToDetect = 1 - 0.15;
    else if (geoId.layer() == 4)
      // probToDetect = 1 - 0.005;
      probToDetect = 1 - 0.05;
    else if (geoId.layer() == 6)
      // probToDetect = 1 - 0.002;
      probToDetect = 1 - 0.02;

    if (0)
      std::cout << "geoId.volume() = " << geoId.volume()
                << "geoId.layer() = " << geoId.layer()
                << " probToDetect = " << probToDetect << " eta = " << eta
                << std::endl;
  }

  // ### Feb 2025: try removal of hits from iTOF layer
  isHitTOF = true;
  // std::cout << ">> lala IA TEST in SimActor" << std::endl;
  if (0 && geoId.sensitive() > 0) {
    // check if we are in iTOF layer
    if (geoId.volume() == 22 && geoId.layer() == 8) {
      double tHit = after.fourPosition().w();
      // int bc = 0; // CONSIDER ONLY THE 1st BC!!! //after.BC();
      int bc = after.BC();  // CONSIDER ONLY THE 1st BC!!!
      isHitTOF = (tHit / 3e11 > (bc * 25e-9) - 5e-9) &&
                 (tHit / 3e11 < (bc * 25e-9) + 5e-9);

      if (0 && fabs(eta) < 0.5)
        std::cout << geoId.sensitive() << " >> " << geoId.volume() << " "
                  << geoId.layer() << " tHit = " << tHit
                  << " BC =  " << after.BC() << " isHitTOF = " << isHitTOF
                  << ".   " << std::endl;
    }
  }
}

}  // namespace ActsFatras::detail

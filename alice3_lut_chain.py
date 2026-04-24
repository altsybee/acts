#!/usr/bin/env python3

"""
Chain to create LUTs for the ALICE3 tracking.
It reads the generated particles and the reconstructed tracks from the simulation output, matches them, and produces LUTs for track parameters as a function of particle kinematics.
"""

from typing import Optional
from job_configs import ChainConfig
import argparse


from acts.examples.root import (
    RootParticleReader,
    RootTrackSummaryReader
)

import pathlib
import acts

from AliceActsPythonBindings import LutMakerAlgorithm
print("LutMakerAlgorithm imported successfully!!!")


def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_path",
                        help="Path to input files",
                        type=str)
    return parser


def addLutMakerAlgorithm(s: acts.examples.Sequencer,
                         inputTracks: str,
                         inputParticles: str,
                         inputTrackParticleMatching: str = "",
                         magneticField: float,
                         # inputMeasurementParticlesMap: str,
                         # outputTrackParticleMatching: str,
                         # outputParticleTrackMatching: str,
                         # looperProtection: str = True,
                         # loop_absEta: str = 1.5,
                         # loop_maxPt: str = 1.,
                         # loop_maxParticleHits: int = 11,
                         logLevel: Optional[acts.logging.Level] = None,
                         lutTag: str = ""):

    lutMakerCfg = LutMakerAlgorithm.Config()
    lutMakerCfg.inputTracks = inputTracks
    lutMakerCfg.inputParticles = inputParticles
    lutMakerCfg.inputTrackParticleMatching = inputTrackParticleMatching
    lutMakerCfg.lutTag = lutTag if lutTag is not None else "lut"
    lutMakerCfg.magneticField = magneticField if magneticField is not None else 0.0
    # lutMakerCfg.inputMeasurementParticlesMap = inputMeasurementParticlesMap
    # lutMakerCfg.outputTrackParticleMatching = outputTrackParticleMatching
    # lutMakerCfg.outputParticleTrackMatching = outputParticleTrackMatching
    # lutMakerCfg.matchingRatio = 1.0
    # lutMakerCfg.doubleMatching = False
    # lutMakerCfg.looperProtection = True
    # lutMakerCfg.loop_absEta = 1.5
    # lutMakerCfg.loop_maxPt = 1.
    # lutMakerCfg.loop_maxParticleHits = 11

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    lutAlg = LutMakerAlgorithm(
        config=lutMakerCfg,
        level=customLogLevel())

    s.addAlgorithm(lutAlg)


def runFullChain(cfg, args):
    """
    Run the full ALICE3 simulation and reconstruction chain.

    Parameters
    ----------
    cfg : ChainConfig.Config
        Configuration object for the chain
    args : argparse.Namespace
        Command line arguments (same format as getArgumentParser().parse_args())
    """

    s = acts.examples.Sequencer(events=10000000, skip=False, numThreads=1)

    # Load the simulation files

    simdir = pathlib.Path(args.input_path)
    s.addReader(RootParticleReader(level=acts.logging.WARNING,
                                   outputParticles="particles_generated_selected",
                                   filePath=simdir / "particles.root"))

    s.addReader(RootTrackSummaryReader(level=acts.logging.WARNING,
                                       outputTracks="ckf_tracks",
                                       outputParticles="summary_particles",
                                       filePath=simdir / "tracksummary_ambi.root"))

    addLutMakerAlgorithm(s,
                         inputTracks="ckf_tracks",
                         inputParticles="particles_generated_selected",
                         inputTrackParticleMatching="",
                         magneticField=cfg.general.MF,
                         lutTag=cfg.geo_dir.geodir.split("/")[-2],
                         # inputMeasurementParticlesMap="measurement_particles_map",
                         # outputTrackParticleMatching="seed_merged_particle_matching",
                         # outputParticleTrackMatching="particle_seed_merged_matching",
                         )

    s.run()


if __name__ == "__main__":
    parser = getArgumentParser()
    args = parser.parse_args()
    runFullChain(cfg=ChainConfig.Config(), args=args)

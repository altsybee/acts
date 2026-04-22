#!/usr/bin/env python3

from typing import Optional
from zipfile import ZipFile
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
    parser = argparse.ArgumentParser(
        description="Command line arguments for full chain setup"
    )

    # General parameters
    parser.add_argument(
        "--nThreads",
        "-nthr",
        dest="nThreads",
        help="Number of threads",
        type=int,
        default=-1
    )
    parser.add_argument(
        "--nEv",
        "-n",
        dest="nEvents",
        help="nEvents",
        type=int,
        default=100000,
    )
    parser.add_argument(
        "--skip",
        "-s",
        dest="skip",
        help="Skip events",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--out_dir_prefix",
        dest="out_dir_prefix",
        help="Dir for output",
        type=str,
        default="test",
    )
    return parser


def addLutMakerAlgorithm(
    s: acts.examples.Sequencer,
    inputTracks: str,
    inputParticles: str,
    inputTrackParticleMatching: str = "",
    # inputMeasurementParticlesMap: str,
    # outputTrackParticleMatching: str,
    # outputParticleTrackMatching: str,
    # looperProtection: str = True,
    # loop_absEta: str = 1.5,
    # loop_maxPt: str = 1.,
    # loop_maxParticleHits: int = 11,
    logLevel: Optional[acts.logging.Level] = None,
):

    lutMakerCfg = LutMakerAlgorithm.Config()
    lutMakerCfg.inputTracks = inputTracks
    lutMakerCfg.inputParticles = inputParticles
    lutMakerCfg.inputTrackParticleMatching = inputTrackParticleMatching
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


def runFullChain(cfg=None, args=None):
    """
    Run the full ALICE3 simulation and reconstruction chain.

    Parameters
    ----------
    cfg : ChainConfig.Config
        Configuration object for the chain
    args : argparse.Namespace
        Command line arguments (same format as getArgumentParser().parse_args())
    """
    # Geometry configuration in case cfg not set
    if (cfg is None):
        cfg = ChainConfig.Config()

    if (args is None):
        parser = getArgumentParser()
        args = parser.parse_args()

    # output directory
    IA_outputDirName = (
        "output/"
        + args.out_dir_prefix
        + "_nEv"
        + str(args.nEvents)
        + "_PID"
        + str(cfg.particleGun.gunPID)
        + "_nMeasMin"
        + str(cfg.tracking.nMeasurementsMin)
        + "_ckfChi2Meas"
        + str(cfg.tracking.ckfChi2Measurement)
        + "_ckfMeasPerSurf"
        + str(cfg.tracking.ckfMeasPerSurf) + "_pion"
    )

    s = acts.examples.Sequencer(events=args.nEvents,
                                skip=args.skip,
                                numThreads=args.nThreads)

    # Load the simulation files

    simdir = pathlib.Path(IA_outputDirName)
    s.addReader(
        RootParticleReader(
            level=acts.logging.WARNING,
            outputParticles="particles_generated_selected",
            filePath=simdir / "particles.root",
        )
    )

    s.addReader(
        RootTrackSummaryReader(
            level=acts.logging.WARNING,
            outputTracks="ckf_tracks",
            outputParticles="summary_particles",
            filePath=simdir / "tracksummary_ambi.root",
        )
    )

    # s.addReader(
    #     RootVertexReader(
    #         level=acts.logging.WARNING,
    #         outputVertices="vertices_truth",
    #         filePath=simdir / "vertices.root",
    #     )
    # )

    addLutMakerAlgorithm(
        s,
        inputTracks="ckf_tracks",
        inputParticles="particles_generated_selected",
        inputTrackParticleMatching="",
        # inputMeasurementParticlesMap="measurement_particles_map",
        # outputTrackParticleMatching="seed_merged_particle_matching",
        # outputParticleTrackMatching="particle_seed_merged_matching",
    )

    s.run()


if __name__ == "__main__":
    parser = getArgumentParser()
    args = parser.parse_args()
    runFullChain(cfg=ChainConfig.Config(), args=args)

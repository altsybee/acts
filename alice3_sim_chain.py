#!/usr/bin/env python3

import os
from datetime import datetime
import shutil
import argparse
import sys
from acts.examples.reconstruction import (
    addSeeding,
    # TruthSeedRanges,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    # ParticleSmearingSigmas,
    addCKFTracks,
    addTruthTrackingGsf,
    #    CKFPerformanceConfig,
    TrackSelectorConfig,
    addKalmanTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    CkfConfig,
    addVertexFitting,
    VertexFinder,
    addSpacePointsMaking,  # May 2025
    addHoughVertexFinding,  # May 2025
)

from acts.examples.simulation import (
    addPythia8,
    addGenParticleSelection,
    addSimParticleSelection,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
    addDigiParticleSelection,
)

from acts.examples.root import (
    RootParticleWriter
    )

import pathlib
import acts
import acts.examples
import acts.examples.geant4

# To automatically unzip stuff
from zipfile import ZipFile

# Move it to an utility tool
def unzipFile(zipfile : pathlib.Path):
    if zipfile.exists():
        # unzip
        with ZipFile(zipfile, 'r') as zip_ref:
            zip_ref.extractall(zipfile.parent)
            print(f"Extracted {zipfile}...")
    else:
        raise FileNotFoundError(f"{zipfile} doesn't exist!")




import alice3.performance.generator as alice3_generator
import alice3.performance.writers as alice3_writers
import alice3.performance.seeding as alice3_seeding
import alice3.performance.plotting as alice3_plotting
import alice3.performance.simulation as alice3_simulation

from job_configs import ChainConfig

u = acts.UnitConstants


def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(
        description="Command line arguments for full chain setup"
    )

    ##### General parameters
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
        default=1,
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

    parser.add_argument(
        "--seed",
        "--randomSeed",
        dest="seed",
        help="seed to use in the random number generator",
        type=int,
        default=42)
    #####
    return parser


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
    ## Geometry configuration in case cfg not set
    if (cfg is None):
        cfg = ChainConfig.Config()

    if (args is None):
        parser = getArgumentParser()
        args = parser.parse_args()

    geo_dir = pathlib.Path(os.path.abspath(os.path.dirname(__file__))) / cfg.general.geo_dir
    print("Loading geometry from ... ",str(geo_dir))
    sys.path.insert(0,str(pathlib.Path(geo_dir).resolve()))

    import buildALICE3Geometry

    ##########################
    ##### SOME OTHER PARAMS
    ##########################

    IA_collisionRegion_forSeeds = (
        250 if cfg.seeding.impParForSeeds < 2.0 else 1000
    )  # mm; large values - for V0 daughter reconstruction

    ### output directory
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
        + str(cfg.tracking.ckfMeasPerSurf)
    )

    outputDir = pathlib.Path.cwd() / IA_outputDirName

    if not outputDir.exists():
        outputDir.mkdir(mode=0o777, parents=True, exist_ok=True)

    detector = buildALICE3Geometry.buildALICE3Geometry(
        geo_dir, cfg.general.enableMaterial, False, acts.logging.INFO
    )
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    if cfg.general.fieldMap != "":
        print(">>> !", cfg.general.fieldMap, "field map is used in this sim+rec run.")
        with open(cfg.general.fieldMap) as magFile:  # Checking if it's RZ or XYZ coordinates
            for l in magFile:
                l = l.strip()
                if l.startswith("#"):
                    continue
                while "  " in l:
                    l = l.replace("  ", " ")
                l = l.split(" ")
                if len(l) == 4:
                    field = acts.MagneticFieldMapRz(cfg.general.fieldMap)
                else:
                    field = acts.MagneticFieldMapXyz(cfg.general.fieldMap)
                break
    else:
        print(">>> !Using Constant B-Field")
        field = acts.ConstantBField(acts.Vector3(0.0, 0.0, cfg.general.MF * u.T))

    rnd = acts.examples.RandomNumbers(seed=args.seed)

    # s = acts.examples.Sequencer(events=nEvents, numThreads=-1)
    s = acts.examples.Sequencer(events=args.nEvents,
                                skip=args.skip,
                                numThreads=args.nThreads)


    if not cfg.pythia.usePythia:
        alice3_generator.addAlice3ParticleGun(s, outputDir=outputDir)
    else:
        s = alice3_generator.addAlice3Pythia8(s, rnd=rnd,
                                              outputDir=outputDir)
            
    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            eta=(-cfg.detSim.particleSelectionEta, cfg.detSim.particleSelectionEta),
            pt=(0.02 * u.MeV, None),
            removeNeutral=False,
            rho=(0.0, 24 * u.mm),
            # absZ=(0.0, 1.0 * u.m),
        ),
    )

    if cfg.detSim.simulation == "Fatras":

        addFatras(
            s,
            trackingGeometry,
            field,
            enableInteractions=True,
            rnd=rnd,
            pMin=cfg.detSim.minFatrasPt,  # GeV, May 2025
            outputDirRoot=outputDir,
            logLevel=acts.logging.INFO,
        )
        
    elif cfg.detSim.simulation == "Geant4":


        alice3_gdml = geo_dir / "o2sim_geometry.gdml"
        
        #Check if we need to unzip
        if not alice3_gdml.exists():
            unzipFile(geo_dir / "o2sim_geometry.gdml.zip")
            
            
        gdml_detector = acts.examples.geant4.GdmlDetector(path=str(alice3_gdml))

        addGeant4(
            s,
            gdml_detector,
            trackingGeometry,
            field,
            volumeMappings=["TRKSensor", "FT3Sensor","ITOFSensor"], 
            outputDirRoot=outputDir,
            rnd=rnd,
            logLevel=acts.logging.INFO,
            killVolume=trackingGeometry.highestTrackingVolume,
            killAfterTime=40. * u.ns,  # 25 * u.ns,
            recordHitsOfSecondaries=True,  # False,#True, ### IA
            keepParticlesWithoutHits=True,  # True, ### IA
            killSecondaries=False,
        )
        
    addSimParticleSelection(
        s,
        ParticleSelectorConfig(
            eta=(-cfg.detSim.particleSelectionEta, cfg.detSim.particleSelectionEta),
            pt = (20 * u.MeV, None),
            removeNeutral=True, # May 2025
            hits=(7, None),  # Nov 2025
            #rho=(0.0, 24 * u.mm),
            #absZ=(0.0, 1.0 * u.m),
        ),
    )
    
    s = addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile= cfg.general.digi_file,
        outputDirRoot=outputDir,
        rnd=rnd,
        logLevel=acts.logging.INFO,
        # doMerge=True,   ##!! 
    )

    alice3_simulation.addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            eta=(-cfg.detSim.particleSelectionEta, cfg.detSim.particleSelectionEta),
            pt=(20 * u.MeV, None),
            removeNeutral=True,
            hits=(7, None),
            rho  = tuple(x * u.mm for x in cfg.detSim.digi_particleSelectionRho),
            absZ = tuple(x * u.mm for x in cfg.detSim.digi_particleSelectionZ),
        ),
        measurementLayers=cfg.detSim.particleSelectionLayers,
    )

    s.run()


if __name__ == "__main__":
    parser = getArgumentParser()
    args = parser.parse_args()
    runFullChain(cfg=ChainConfig.Config(), args=args)


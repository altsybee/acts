from pathlib import Path
from typing import Optional, Union, List
from enum import Enum
from collections import namedtuple
import os
import array as arr

import acts
from acts import UnitConstants as u
import acts.examples
import acts.examples.reconstruction as acts_reco

import alice3.performance.plotting as alice3_plotting

from acts.examples.reconstruction import (
    addSeeding,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    # addSeedFilterML,
    # SeedFilterMLDBScanConfig,
)

from enum import Enum, auto

class SeedingMode(Enum):
    HIGH_PT = auto()
    LOW_PT = auto()
    LOOPERS = auto()
    STRANGENESS = auto()


#############################################
#### IGOR's SEEDING: SOME OTHER PARAMS   ####
#############################################

collisionRegion_forSeeds = 250 # 1000 # mm; large values - for V0 daughter reconstruction
seedingAlg = SeedingAlgorithm.GridTriplet
maxSeedsPerSpMConf = 1 
maxQualitySeedsPerSpMConf = 1
enableMaterial = True
seedingGeoSelection_dir = Path(os.environ["MAIN_DIR"])
print(seedingGeoSelection_dir)

##Parameters to be passed
maxSeedsPerSpM   = 2     ## ok for pp, for PbPb higher
sigmaScattering  = 5.0   ## Tuned with particle gun, combine radLengthPerSeed = 0.05
radLengthPerSeed = 0.05
minSeedPt        = 0.2  ## to be changed per iteration
impParForSeeds   = 1.0   #mm
bFieldInZ        = 2.0   #T


##############################################
#### PAVEL's SEEDING: SOME OTHER PARAMS   ####
##############################################
seedParamOption = 0
deltaRmin = arr.array("f", [1., 1.0592])
deltaRmax = arr.array("f", [60., 61.587])
maxSeedsPerMiddleSp = arr.array("I", [1, 8])
maxSigmaScattering = arr.array("f", [5., 3.163])
radiationLengthPerSeed = arr.array("f", [0.1, 0.0785])
maxImpact = arr.array("f", [3., 22.29])
maxCotTheta = arr.array("f", [27.2899, 32.57108])
enableSeedconfirmation = False
initVarInflFactor = 1.0
PavelInitialVarInflation  = [initVarInflFactor] * 6
enableMat = True



def get_seed_finder_config(iteration: int = 0, **overrides) -> SeedFinderConfigArg:
    """
    Get seeding configuration for a given iteration.
    Args:
    iteration: Iteration number (0,1,2)
    overriders: kwargs to override iteration-specific defaults
    """
    
    iteration_params = {
        0 : {
            'minPt' : 0.4 * u.GeV, 
            'impactMax' : 1.0 * u.mm,
            'z': (-1000 * u.mm, 1000 * u.mm),
        },
        
        1 : {
            'minPt' : 0.15 * u.GeV, 
            'impactMax' : 1.0 * u.mm,
            'z': (-1000 * u.mm, 1000 * u.mm),
        },

        2 : {
            'minPt' : 0.08 * u.GeV,
            'impactMax' : 1.0 * u.mm,
            'z': (-1000 * u.mm, 1000 * u.mm),
            },
        
        3 : {
            'minPt' : 0.05 * u.GeV,
            'impactMax' : 1.0 * u.mm,
            'z': (-1000 * u.mm, 1000 * u.mm),
            },
        
    }
    
    params = {
        # --- default parameters ---
        "r"               : (None, 210 * u.mm),
        "deltaR"          : (1 * u.mm, 200 * u.mm),
        "deltaRBottomSP"  : (3 * u.mm, 50 * u.mm),
        "collisionRegion" : (
            -collisionRegion_forSeeds * u.mm,
            collisionRegion_forSeeds * u.mm,
        ),
        "maxSeedsPerSpM" : maxSeedsPerSpM,
        "sigmaScattering": sigmaScattering,
        "radLengthPerSeed": radLengthPerSeed,
        "cotThetaMax": 27.2899,
        "seedConfirmation": True,
        "centralSeedConfirmationRange": acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-620 * u.mm,
            zMaxSeedConf=620 * u.mm,
            rMaxSeedConf=4.9
            * u.mm,  # 36 * u.mm,  # IA: dramatically affects acceptance at eta ~4. <5 * u.mm  gives best results
            nTopForLargeR=1,  # number of top space points that confirm my seed at larger R, 1 - no confirmation
            nTopForSmallR=2,
        ),
        "forwardSeedConfirmationRange": acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-1220 * u.mm,
            zMaxSeedConf=1220 * u.mm,
            rMaxSeedConf=26 * u.mm,  # 15 * u.mm,  #36 * u.mm,
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),
        "useVariableMiddleSPRange": True,
        "deltaRMiddleSPRange": (0.2 * u.mm, 1.0 * u.mm),

        # Iteration 0 defaults that will be overridden later
        **iteration_params[0]
        
        }

    print("PF::params before update ", params)
    
    # Layer on iteration-specific values
    params.update(iteration_params.get(iteration, {}))

    # Layer on any explicit call-site overrides (highest priority)
    params.update(overrides)


    print("PF::params updated ",params)

    return SeedFinderConfigArg(**params)


DefaultSeedFinderConfigArg = SeedFinderConfigArg(
        r=(None, 210 * u.mm),  # iTOF is at 190 mm! if we want it for seeding
        # r=(None, 150 * u.mm),
        # r=(None, 30 * u.mm),
        deltaR=(1 * u.mm, 200 * u.mm),  # deltaR=(1. * u.mm, 60 * u.mm),
        collisionRegion=(
            -collisionRegion_forSeeds * u.mm,
            collisionRegion_forSeeds * u.mm,
        ),
        z=(-1000 * u.mm, 1000 * u.mm),
        maxSeedsPerSpM=maxSeedsPerSpM,  # 2 is minimum, >2 is better for Pb-Pb
        sigmaScattering=sigmaScattering,
        radLengthPerSeed=radLengthPerSeed,  
        minPt=minSeedPt * u.GeV,
        impactMax=impParForSeeds * u.mm,  # important! IB vs ML seeds (e.g. 1 mm is ok for IB seeds, 5 mm - for ML seeds)
        cotThetaMax=27.2899,
        seedConfirmation=True,
        centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-620 * u.mm,
            zMaxSeedConf=620 * u.mm,
            rMaxSeedConf=4.9
            * u.mm,  # 36 * u.mm,  # IA: dramatically affects acceptance at eta ~4. <5 * u.mm  gives best results
            nTopForLargeR=1,  # number of top space points that confirm my seed at larger R, 1 - no confirmation
            nTopForSmallR=2,
        ),
        forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-1220 * u.mm,
            zMaxSeedConf=1220 * u.mm,
            rMaxSeedConf=26 * u.mm,  # 15 * u.mm,  #36 * u.mm,
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),
        # skipPreviousTopSP=True,
        useVariableMiddleSPRange=True,
        deltaRMiddleSPRange=(0.2 * u.mm, 1.0 * u.mm),
    )

DefaultSeedFinderOptionsArg = SeedFinderOptionsArg(bFieldInZ=bFieldInZ * u.T, beamPos=(0 * u.mm, 0 * u.mm))

DefaultSeedFilterConfigArg = SeedFilterConfigArg(
        seedConfirmation=True if impParForSeeds < 2.0 else False,  # mm
        # If seedConfirmation is true we classify seeds as "high-quality" seeds.
        # Seeds that are not confirmed as "high-quality" are only selected if no
        # other "high-quality" seed has been found for that inner-middle doublet
        # Maximum number of normal seeds (not classified as "high-quality" seeds)
        # in seed confirmation
        maxSeedsPerSpMConf=maxSeedsPerSpMConf,  # 1 - USED_FOR_AUG_2025,#3, # CRUCIAL!!!!!!
        maxQualitySeedsPerSpMConf=maxQualitySeedsPerSpMConf,  # 1 - USED_FOR_AUG_2025,   # Core/include/Acts/Seeding/SeedFilterConfig.hpp
        # Maximum number of "high-quality" seeds for each inner-middle SP-dublet in
        # seed confirmation. If the limit is reached we check if there is a normal
        # quality seed to be replaced
    )


DefaultSpacePointGridConfigArg = SpacePointGridConfigArg(
        # zBinEdges=[
        # -4000.0,
        # -2500.0,
        # -2000.0,
        # -1320.0,
        # -625.0,
        # -350.0,
        # -250.0,
        # 250.0,
        # 350.0,
        # 625.0,
        # 1320.0,
        # 2000.0,
        # 2500.0,
        # 4000.0,
        # ],
        impactMax=1.0 * u.mm,
        phiBinDeflectionCoverage=3, #3
    )
DefaultSeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(
    # zBinNeighborsTop=[
    # [0, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 0],
    # ],
    # zBinNeighborsBottom=[
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 1],
    # [0, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # [-1, 0],
    # ],
    numPhiNeighbors=1,
)


def addSpacePointsMaking(
        sequence : acts.examples.Sequencer,
        trackingGeometry : acts.TrackingGeometry,
        geoSelectionConfigFile: Union[Path, str],
        stripGeoSelectionConfigFile: Union[Path, str],
        inputMeasurements : str,
        outputSpacePoints : str,
        logLevel: acts.logging.Level = None,
        ):
    """adds space points making
    For parameters description see addSeeding
    """
    logLevel = acts.examples.defaultLogging(sequence, logLevel)()
    #logLevel=acts.logging.DEBUG
    spAlg = acts.examples.SpacePointMaker(
        level=logLevel,
        inputMeasurements=inputMeasurements,
        outputSpacePoints=outputSpacePoints,
        trackingGeometry=trackingGeometry,
        geometrySelection=acts.examples.json.readJsonGeometryList(
            str(geoSelectionConfigFile)
        ),
        stripGeometrySelection=(
            acts.examples.json.readJsonGeometryList(str(stripGeoSelectionConfigFile))
            if stripGeoSelectionConfigFile
            else []
        ),
    )
    sequence.addAlgorithm(spAlg)
    return spAlg.config.outputSpacePoints

def _make_prefix(prefix: str = "", iterationIndex: int = 0) -> str:
    """
    ACTS v46.4.0 multi-pass convention.

    For prefix="iter_1_", addSeeding/addCKFTracks should use:
      iter_1_measurement_subset
      iter_1_spacepoints
      iter_1_seeds
      iter_1_estimatedparameters
      iter_1_estimatedseeds
      iter_1_seed-prototracks
      iter_1_seed-tracks

    iterationIndex is kept only for backward compatibility.
    """
    if prefix:
        return prefix

    if iterationIndex > 0:
        return f"iter_{iterationIndex}_"

    return ""

def addSeeding(
        s: acts.examples.Sequencer,
        trackingGeometry: acts.TrackingGeometry,
        field: acts.MagneticFieldProvider,
        geoSelectionConfigFile: Optional[Union[Path, str]] = None,
        stripGeoSelectionConfigFile: Optional[Union[Path, str]] = None,
        layerMappingConfigFile: Optional[Union[Path, str]] = None,
        connectorInputConfigFile: Optional[Union[Path, str]] = None,
        lutInputConfigFile: Optional[Union[Path, str]] = None,
        seedingAlgorithm: SeedingAlgorithm = SeedingAlgorithm.GridTriplet,
        trackSmearingSigmas: acts.examples.reconstruction.TrackSmearingSigmas = acts.examples.reconstruction.TrackSmearingSigmas(),
        initialSigmas: Optional[list] = None,
        initialSigmaQoverPt: Optional[float] = None,
        initialSigmaPtRel: Optional[float] = None,
        initialVarInflation: Optional[list] = None,
        seedFinderConfigArg: SeedFinderConfigArg = SeedFinderConfigArg(),
        seedFinderOptionsArg: SeedFinderOptionsArg = SeedFinderOptionsArg(),
        seedFilterConfigArg: SeedFilterConfigArg = SeedFilterConfigArg(),
        spacePointGridConfigArg: SpacePointGridConfigArg = SpacePointGridConfigArg(),
        seedingAlgorithmConfigArg: SeedingAlgorithmConfigArg = SeedingAlgorithmConfigArg(),
        houghTransformConfig: acts.examples.HoughTransformSeeder.Config = acts.examples.HoughTransformSeeder.Config(),
        adaptiveHoughTransformConfig: Optional[
            acts.examples.AdaptiveHoughTransformSeeder.Config
        ] = None,
        hashingTrainingConfigArg: Optional[
            acts.examples.reconstruction.HashingTrainingConfigArg
        ] = acts.examples.reconstruction.HashingTrainingConfigArg(),
        hashingAlgorithmConfigArg: Optional[
            acts.examples.reconstruction.HashingAlgorithmConfigArg
        ] = acts.examples.reconstruction.HashingAlgorithmConfigArg(),
        truthEstimatedSeedingAlgorithmConfigArg: acts.examples.reconstruction.TruthEstimatedSeedingAlgorithmConfigArg = acts.examples.reconstruction.TruthEstimatedSeedingAlgorithmConfigArg(),
        particleHypothesis: Optional[
            acts.ParticleHypothesis
        ] = acts.ParticleHypothesis.pion,
        inputParticles: str = "particles",
        selectedParticles: str = "particles_selected",
        outputDirRoot: Optional[Union[Path, str]] = None,
        outputDirCsv: Optional[Union[Path, str]] = None,
        logLevel: Optional[acts.logging.Level] = None,
        rnd: Optional[acts.examples.RandomNumbers] = None,
        inputMeasurements: str = "measurement_subset",
        outputSpacePoints: str = "spacepoints",
        outputSeeds: str = "seeds",
        iterationIndex: int = 0,
        prefix: str = "",
        
        
) -> None:
    """This function steers the seeding
    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    geoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection. Not required for SeedingAlgorithm.TruthSmeared.
    stripGeoSelectionConfigFile : Path|str, path, None
        Json file for space point geometry selection in strips. Needed for SpacePoint making.
    seedingAlgorithm : SeedingAlgorithm, Default
        seeding algorithm to use: one of Default (no truth information used), TruthSmeared, TruthEstimated
    trackSmearingSigmas : TrackSmearingSigmas(loc0, loc0PtA, loc0PtB, loc1, loc1PtA, loc1PtB, time, phi, theta, ptRel)
        TrackSmearing configuration.
        Defaults specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackParameterSmearing.hpp
    initialSigmas : list
        Sets the initial covariance matrix diagonal. This is ignored in case of TruthSmearing.
        Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp
    initialVarInflation : list
        List of 6 scale factors to inflate the initial covariance matrix
        Defaults (all 1) specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackParameterSmearing.hpp
    seedFinderConfigArg : SeedFinderConfigArg(maxSeedsPerSpM, cotThetaMax, sigmaScattering, radLengthPerSeed, minPt, impactMax, deltaPhiMax, interactionPointCut, deltaZMax, maxPtScattering, zBinEdges, zBinsCustomLooping, rRangeMiddleSP, useVariableMiddleSPRange, binSizeR, seedConfirmation, centralSeedConfirmationRange, forwardSeedConfirmationRange, deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z)
        SeedFinderConfig settings. deltaR, deltaRBottomSP, deltaRTopSP, deltaRMiddleSPRange, collisionRegion, r, z.
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFinderOptionsArg :  SeedFinderOptionsArg(bFieldInZ, beamPos)
        Defaults specified in Core/include/Acts/Seeding/SeedFinderConfig.hpp
    seedFilterConfigArg : SeedFilterConfigArg(compatSeedWeight, compatSeedLimit, numSeedIncrement, seedWeightIncrement, seedConfirmation, maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf, useDeltaRorTopRadius)
                                Defaults specified in Core/include/Acts/Seeding/SeedFilterConfig.hpp
    spacePointGridConfigArg : SpacePointGridConfigArg(rMax, zBinEdges, phiBinDeflectionCoverage, phi, maxPhiBins, impactMax)
                                SpacePointGridConfigArg settings. phi is specified as a tuple of (min,max).
        Defaults specified in Core/include/Acts/Seeding/SpacePointGrid.hpp
    seedingAlgorithmConfigArg : SeedingAlgorithmConfigArg(allowSeparateRMax, zBinNeighborsTop, zBinNeighborsBottom, numPhiNeighbors, useExtraCuts)
                                Defaults specified in Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/SeedingAlgorithm.hpp
    hashingTrainingConfigArg : HashingTrainingConfigArg(annoySeed, f)
                                Defaults specified in Plugins/Hashing/include/ActsPlugins/Hashing/HashingTrainingConfig.hpp
    hashingAlgorithmConfigArg : HashingAlgorithmConfigArg(bucketSize, zBins, phiBins)
                                Defaults specified in Plugins/Hashing/include/ActsPlugins/Hashing/HashingAlgorithmConfig.hpp
    truthEstimatedSeedingAlgorithmConfigArg : TruthEstimatedSeedingAlgorithmConfigArg(deltaR)
        Currently only deltaR=(min,max) range specified here.
    particleHypothesis : Optional[acts.ParticleHypothesis]
        The hypothesis used for track finding. Defaults to pion.
    inputParticles : str, "particles"
        input particles name in the WhiteBoard
    selectedParticles : str, "particles_selected"
        selected particles name in the WhiteBoard
    outputDirRoot : Path|str, path, None
        the output folder for ROOT output, None triggers no output
    logLevel : acts.logging.Level, None
        logging level to override setting given in `s`
    rnd : RandomNumbers, None
        random number generator. Only used by SeedingAlgorithm.TruthSmeared.
    """

    collection_prefix = _make_prefix(prefix=prefix, iterationIndex=iterationIndex)

    outputSpacePointsWithPrefix = collection_prefix + outputSpacePoints
    outputSeedsWithPrefix = collection_prefix + outputSeeds

    outputEstimatedParameters = collection_prefix + "estimatedparameters"
    outputEstimatedSeeds = collection_prefix + "estimatedseeds"
    outputProtoTracks = collection_prefix + "seed-prototracks"
    outputSeedTracks = collection_prefix + "seed-tracks"

    outputSeedParticleMatching = collection_prefix + "seed_particle_matching"
    outputParticleSeedMatching = collection_prefix + "particle_seed_matching"

    trackFinderWriterOutName = collection_prefix + "performance_seeding.root"
    trackParamsWriterOutName = collection_prefix + "estimatedparameters.root"
    
    
    logLevel = acts.examples.defaultLogging(s, logLevel)()
    #logLevel = acts.logging.VERBOSE
    
    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if seedingAlgorithm == SeedingAlgorithm.TruthSmeared:
        acts_reco.addTruthSmearedSeeding(
            s=s,
            rnd=rnd,
            selectedParticles=selectedParticles,
            trackSmearingSigmas=trackSmearingSigmas,
            initialSigmas=initialSigmas,
            initialSigmaQoverPt=initialSigmaQoverPt,
            initialSigmaPtRel=initialSigmaPtRel,
            initialVarInflation=initialVarInflation,
            particleHypothesis=particleHypothesis,
            logLevel=logLevel,
        )
    else:

        spacePoints = addSpacePointsMaking(
            sequence = s,
            trackingGeometry = trackingGeometry,
            geoSelectionConfigFile = geoSelectionConfigFile,
            stripGeoSelectionConfigFile = stripGeoSelectionConfigFile,
            inputMeasurements = inputMeasurements,
            outputSpacePoints=outputSpacePointsWithPrefix,
            logLevel=logLevel,
            
        )
        seeds = None
        perSeedParticleHypothesis = None
        # Run either: truth track finding or seeding
        if seedingAlgorithm == SeedingAlgorithm.TruthEstimated:
            seeds, perSeedParticleHypothesis = acts_reco.addTruthEstimatedSeeding(
                s,
                spacePoints,
                selectedParticles,
                truthEstimatedSeedingAlgorithmConfigArg,
                particleHypothesis=particleHypothesis,
                logLevel=logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.HoughTransform:
            houghTransformConfig.inputSpacePoints = [spacePoints]
            houghTransformConfig.inputMeasurements = inputMeasurements
            houghTransformConfig.outputProtoTracks = collection_prefix + "prototracks"
            houghTransformConfig.outputSeeds = outputSeedsWithPrefix
            houghTransformConfig.trackingGeometry = trackingGeometry
            seeds = acts_reco.addHoughTransformSeeding(s, houghTransformConfig, logLevel)

        elif seedingAlgorithm == SeedingAlgorithm.AdaptiveHoughTransform:
            if adaptiveHoughTransformConfig is None:
                adaptiveHoughTransformConfig = (
                    acts.examples.AdaptiveHoughTransformSeeder.Config()
                )

            adaptiveHoughTransformConfig.inputSpacePoints = [spacePoints]
            adaptiveHoughTransformConfig.outputProtoTracks = (
                collection_prefix + "prototracks"
            )
            adaptiveHoughTransformConfig.outputSeeds = outputSeedsWithPrefix
            adaptiveHoughTransformConfig.trackingGeometry = trackingGeometry
            adaptiveHoughTransformConfig.threshold = 4
            adaptiveHoughTransformConfig.noiseThreshold = 12
            adaptiveHoughTransformConfig.phiMinBinSize = 3.14 / (2.0 * 257.0)
            adaptiveHoughTransformConfig.qOverPtMinBinSize = 1.1 / (2.0 * 257.0)
            adaptiveHoughTransformConfig.qOverPtMin = 1.1
            adaptiveHoughTransformConfig.doSecondPhase = True
            adaptiveHoughTransformConfig.zMinBinSize = 1 * u.mm
            adaptiveHoughTransformConfig.cotThetaMinBinSize = 0.1
            adaptiveHoughTransformConfig.deduplicate = True
            seeds = acts_reco.addAdaptiveHoughTransformSeeding(
                s, adaptiveHoughTransformConfig, logLevel=logLevel
            )
        elif seedingAlgorithm == SeedingAlgorithm.Gbts:
            # output of algs changed, only one output now
            seeds = acts_reco.addGbtsSeeding(
                s,
                spacePoints,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                trackingGeometry,
                logLevel,
                layerMappingConfigFile,
                geoSelectionConfigFile,
                connectorInputConfigFile,
                lutInputConfigFile,
            )
        elif seedingAlgorithm == SeedingAlgorithm.Hashing:
            seeds, buckets = acts_reco.addHashingSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                hashingTrainingConfigArg,
                hashingAlgorithmConfigArg,
                logLevel,
            )
        elif seedingAlgorithm == SeedingAlgorithm.GridTriplet:
            seeds = acts.examples.reconstruction.addGridTripletSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
                outputSeeds=outputSeedsWithPrefix,
            )
            
        elif seedingAlgorithm == SeedingAlgorithm.OrthogonalTriplet:
            seeds = acts_reco.addOrthogonalTripletSeeding(
                s,
                spacePoints,
                seedingAlgorithmConfigArg,
                seedFinderConfigArg,
                seedFinderOptionsArg,
                seedFilterConfigArg,
                spacePointGridConfigArg,
                logLevel,
            )
        else:
            pass

        parEstimateAlg = acts.examples.TrackParamsEstimationAlgorithm(
            level=logLevel,
            inputSeeds=seeds,
            inputParticleHypotheses=perSeedParticleHypothesis,
            outputTrackParameters=outputEstimatedParameters,
            outputSeeds=outputEstimatedSeeds,
            trackingGeometry=trackingGeometry,
            magneticField=field,
            **acts.examples.defaultKWArgs(
                initialSigmas=initialSigmas,
                initialSigmaQoverPt=initialSigmaQoverPt,
                initialSigmaPtRel=initialSigmaPtRel,
                initialVarInflation=initialVarInflation,
                particleHypothesis=particleHypothesis,
            ),
        )
        s.addAlgorithm(parEstimateAlg)

        s.addAlgorithm(
            acts.examples.SeedsToProtoTracks(
                level=logLevel,
                inputSeeds=outputEstimatedSeeds,
                outputProtoTracks=outputProtoTracks,
            )
        )

        s.addAlgorithm(
            acts.examples.ProtoTracksToTracks(
                level=logLevel,
                inputProtoTracks=outputProtoTracks,
                inputTrackParameters=outputEstimatedParameters,
                inputMeasurements=inputMeasurements,
                outputTracks=outputSeedTracks,
            )
        )

        print("PF::DEBUG collection_prefix", collection_prefix)
        print("PF::DEBUG inputMeasurements", inputMeasurements)
        print("PF::DEBUG measurement_particles_map")
        print("PF::DEBUG", outputSeedParticleMatching)
        print("PF::DEBUG", outputParticleSeedMatching)

        inputMeasurementParticlesMap = "measurement_particles_map"
        inputParticleMeasurementMap = "particle_measurements_map"
        
        s.addAlgorithm(
            acts.examples.TrackTruthMatcher(
                level=logLevel,
                inputTracks=outputSeedTracks,
                inputParticles=selectedParticles,
                inputMeasurementParticlesMap=inputMeasurementParticlesMap,
                outputTrackParticleMatching=outputSeedParticleMatching,
                outputParticleTrackMatching=outputParticleSeedMatching,
                matchingRatio=1.0,
                doubleMatching=False,
            )
        )

        if outputDirRoot is not None:
            addSeedPerformanceWriters(
                s,
                outputDirRoot,
                outputSeedTracks,
                outputProtoTracks,
                selectedParticles,
                inputParticles,
                inputParticleMeasurementMap,
                outputSeedParticleMatching,
                outputParticleSeedMatching,
                parEstimateAlg.config.outputTrackParameters,
                trackFinderWriterOutName,
                trackParamsWriterOutName,
                logLevel,
            )

    return s


# This adds only the performance algorithms to check the performance on a set of seed-tracks
# - Schedule the TruthMatcher
# - Schedule the TrackFinderPerformanceWriter for plots

def addSeedPerformanceAlgorithms(
        sequence : acts.examples.Sequencer,
        outputDirRoot : Union[Path, str],
        tracks: str,
        inputParticles : str,
        inputMeasurementParticlesMap : str,
        outputTrackParticleMatching : str,
        outputParticleTrackMatching : str,
        logLevel: acts.logging.Level = None,
        ) :

    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    
    sequence.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=customLogLevel(),
            inputTracks=tracks,
            inputParticles=inputParticles,
            inputMeasurementParticlesMap=inputMeasurementParticlesMap,
            outputTrackParticleMatching=outputTrackParticleMatching,
            outputParticleTrackMatching=outputParticleTrackMatching,
            matchingRatio=1.0,
            doubleMatching=False,
        )
    )
    
    
def addSeedPerformanceWriters(
        sequence: acts.examples.Sequencer,
        outputDirRoot: Union[Path, str],
        tracks: str,
        prototracks: str,
        selectedParticles: str,
        inputParticles: str,
        inputParticleMeasurementsMap : str,
        inputTrackParticleMatching : str,
        inputParticleTrackMatching : str,
        outputTrackParameters: str,
        trackFinderWriterOutName : str,
        trackParamsWriterOutName : str,
        logLevel: acts.logging.Level = None,
):
    """Writes seeding related performance output"""
    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir(parents=True, exist_ok=True)

    
    print("PF:: Adding RootTrackFinderPerformanceWriter on tracks: ", tracks)
    
    sequence.addWriter(
        acts.examples.root.RootTrackFinderPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks=tracks,
            inputParticles=selectedParticles,
            inputTrackParticleMatching=inputTrackParticleMatching,
            inputParticleTrackMatching=inputParticleTrackMatching,
            inputParticleMeasurementsMap=inputParticleMeasurementsMap,
            effPlotToolConfig = alice3_plotting.effPlotToolConfig,
            fakePlotToolConfig = alice3_plotting.fakePlotToolConfig,
            filePath=str(outputDirRoot / trackFinderWriterOutName),
        )
    )

    sequence.addWriter(
        acts.examples.root.RootTrackParameterWriter(
            level=customLogLevel(),
            inputTrackParameters=outputTrackParameters,
            inputProtoTracks=prototracks,
            inputParticles=inputParticles,
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDirRoot / trackParamsWriterOutName),
            treeName="estimatedparams",
        )
    )



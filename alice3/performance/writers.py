from pathlib import Path
from typing import Optional, Union, List, Literal
from enum import Enum
from collections import namedtuple
import pathlib

import acts
import acts.examples
import acts.examples.reconstruction as acts_reco
from acts import UnitConstants as u


from acts.examples.reconstruction import (
    TrackSelectorConfig,
    CkfConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    )

#Load config
from job_configs import ChainConfig
cfg = ChainConfig.Config()

#Alice3 specific algorithms
from AliceActsPythonBindings import TrackTruthMatcher
from AliceActsPythonBindings import HitRemoverAlgorithm
from AliceActsPythonBindings import TrackMergerAlgorithm
#from AliceActsPythonBindings import RootTrackFitterPerformanceWriter

#Alice3 seeding
import alice3.performance.seeding as alice3_seeding

# Alice3 plotting
import alice3.performance.plotting as alice3_plotting


def get_track_selector_config(iteration: int = 0, **overrides) -> TrackSelectorConfig:

    """
    Get track selector configuration for a given iteration
    Args:
    iteration: Iteration number
    overrides: kwargs to override iteration-specific defaults
    """

    iteration_params = {
        0 : {
            'pt' : (0.05 * u.GeV,None),
            },
        
        1 : {
            'pt' : (0.05 * u.GeV, None),
            },
        
        2 : {
            'pt' : (0.05 * u.GeV, None),
            },
        
        3 : {
            'pt' : (0.05 * u.GeV, None),
            },
        }
    
    params = {
        "nMeasurementsMin" : cfg.tracking.nMeasurementsMin,
        "maxHoles"         : 2,
        "maxOutliers"      : 2,
        "maxSharedHits"    : 2,
        
        
        **iteration_params[0]
    }

    params.update(iteration_params.get(iteration,{}))

    params.update(overrides)

    return TrackSelectorConfig(**params)


def addCKFTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    trackSelectorConfig: Optional[
        Union[acts_reco.TrackSelectorConfig, List[acts_reco.TrackSelectorConfig]]
    ] = None,
    ckfConfig: acts_reco.CkfConfig = acts_reco.CkfConfig(),
    twoWay: bool = True,
    reverseSearch: bool = False,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeTrackSummary: bool = True,
    writeTrackStates: bool = False,
    writePerformance: bool = True,
    writeCovMat: bool = False,
    inputMeasurements: str = "measurement_subset",
    inputMeasurementParticlesMap: str = "measurement_particles_map",
    inputParticleMeasurementsMap: str = "particle_measurements_map",
    inputInitialTrackParameters: str = "estimatedparameters",
    inputSeeds: str = "estimatedseeds",
    outputTracks: str = "ckf_tracks",
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the seeding

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Seeding steps (returned from addSeeding)
    trackingGeometry : tracking geometry
    field : magnetic field
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for ROOT output, None triggers no output
    trackSelectorConfig : TrackSelectorConfig(loc0, loc1, time, eta, absEta, pt, phi, minMeasurements)
        TrackSelector configuration. Each range is specified as a tuple of (min,max).
        Specify as a list(TrackSelectorConfig) for eta-dependent cuts, with binning specified by absEta[1].
        Defaults of no cuts specified in Examples/Algorithms/TruthTracking/ActsExamples/TruthTracking/TrackSelector.hpp
    writeTrackSummary : bool, True
        write tracksummary_ckf.root ntuple?
    writeTrackStates : bool, False
        write trackstates_ckf.root ntuple? This can be quite large.
    writePerformance : bool, True
        write performance_fitting_ckf.root and performance_finding_ckf.root ntuples?
    writeCovMat : bool, False
        write covaraiance matrices to tracksummary_ckf.root ntuple?
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    tslist = (
        []
        if trackSelectorConfig is None
        else (
            [trackSelectorConfig]
            if isinstance(trackSelectorConfig, acts_reco.TrackSelectorConfig)
            else trackSelectorConfig
        )
    )

    if len(tslist) > 1:
        cutSets = []
        for c in tslist:
            defKW = acts_reco.trackSelectorDefaultKWArgs(c)
            defKW.pop("absEtaMax", None)
            cutSets += [acts.TrackSelector.Config(**defKW)]
    else:
        cutSets = [
            acts.TrackSelector.Config(**acts_reco.trackSelectorDefaultKWArgs(c)) for c in tslist
        ]

    if len(tslist) == 0:
        trkSelCfg = None
    elif len(tslist) == 1:
        trkSelCfg = cutSets[0]
    else:
        trkSelCfg = acts.TrackSelector.EtaBinnedConfig(
            cutSets=cutSets,
            absEtaEdges=[cutSets[0].absEtaMin] + [c.absEta[1] for c in tslist],
        )

    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=customLogLevel(),
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [
                (
                    acts.GeometryIdentifier(),
                    (
                        [],
                        [ckfConfig.chi2CutOffMeasurement],
                        [ckfConfig.chi2CutOffOutlier],
                        [ckfConfig.numMeasurementsCutOff],
                    ),
                )
            ]
        ),
        inputMeasurements=inputMeasurements,
        inputInitialTrackParameters=inputInitialTrackParameters,
        inputSeeds=(
            inputSeeds 
            if ckfConfig.seedDeduplication or ckfConfig.stayOnSeed 
            else ""
        ),
        outputTracks=outputTracks,
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field, customLogLevel()
        ),
        **acts.examples.defaultKWArgs(
            trackingGeometry=trackingGeometry,
            magneticField=field,
            trackSelectorCfg=trkSelCfg,
            maxSteps=ckfConfig.maxSteps,
            twoWay=twoWay,
            reverseSearch=reverseSearch,
            seedDeduplication=ckfConfig.seedDeduplication,
            stayOnSeed=ckfConfig.stayOnSeed,
            pixelVolumeIds=ckfConfig.pixelVolumes,
            stripVolumeIds=ckfConfig.stripVolumes,
            maxPixelHoles=ckfConfig.maxPixelHoles,
            maxStripHoles=ckfConfig.maxStripHoles,
            trimTracks=ckfConfig.trimTracks,
            constrainToVolumeIds=ckfConfig.constrainToVolumes,
            endOfWorldVolumeIds=ckfConfig.endOfWorldVolumes,
        ),
    )
    s.addAlgorithm(trackFinder)

    if outputTracks == "ckf_tracks":
        s.addWhiteboardAlias("tracks", trackFinder.config.outputTracks)

    truthMatchCfg = TrackTruthMatcher.Config()
    truthMatchCfg.inputTracks = outputTracks
    truthMatchCfg.inputParticles = "particles_selected"
    truthMatchCfg.inputMeasurementParticlesMap = inputMeasurementParticlesMap
    truthMatchCfg.outputTrackParticleMatching = f"{outputTracks}_particle_matching"
    truthMatchCfg.outputParticleTrackMatching = f"particle_{outputTracks}_matching"
    truthMatchCfg.matchingRatio = 1.0
    truthMatchCfg.doubleMatching = False
    truthMatchCfg.looperProtection = True
    truthMatchCfg.loop_absEta = 1.5
    truthMatchCfg.loop_maxPt = 0.2
    truthMatchCfg.loop_maxParticleHits = 11

    matchAlg = TrackTruthMatcher(
        config=truthMatchCfg,
        level=customLogLevel(),
    )
    s.addAlgorithm(matchAlg)

    if outputTracks == "ckf_tracks":
        s.addWhiteboardAlias(
            "track_particle_matching", truthMatchCfg.outputTrackParticleMatching
        )
        s.addWhiteboardAlias(
            "particle_track_matching", truthMatchCfg.outputParticleTrackMatching
        )

    addTrackWriters(
        s,
        name=trackFinder.config.outputTracks,
        tracks=trackFinder.config.outputTracks,
        inputTrackParticleMatching=truthMatchCfg.outputTrackParticleMatching,
        inputParticleTrackMatching=truthMatchCfg.outputParticleTrackMatching,
        inputParticleMeasurementsMap=inputParticleMeasurementsMap,
        outputDirCsv=outputDirCsv,
        outputDirRoot=outputDirRoot,
        writeSummary=writeTrackSummary,
        writeStates=writeTrackStates,
        writeFitterPerformance=False,
        writeFinderPerformance=writePerformance,
        writeCovMat=writeCovMat,
        logLevel=logLevel,
    )

    return s


def addTrackTruthMatcher(
        s: acts.examples.Sequencer,
        inputTracks: str,
        inputParticles: str,
        inputMeasurementParticlesMap: str,
        outputTrackParticleMatching: str,
        outputParticleTrackMatching: str,
        looperProtection: bool = True,
        loop_absEta: float = 1.5,
        loop_maxPt: float = 1.0,
        loop_maxParticleHits: int = 11,
        logLevel: Optional[acts.logging.Level] = None,
):
    truthMatchCfg = TrackTruthMatcher.Config()
    truthMatchCfg.inputTracks = inputTracks
    truthMatchCfg.inputParticles = inputParticles
    truthMatchCfg.inputMeasurementParticlesMap = inputMeasurementParticlesMap
    truthMatchCfg.outputTrackParticleMatching = outputTrackParticleMatching
    truthMatchCfg.outputParticleTrackMatching = outputParticleTrackMatching
    truthMatchCfg.matchingRatio = 1.0
    truthMatchCfg.doubleMatching = False
    truthMatchCfg.looperProtection = looperProtection
    truthMatchCfg.loop_absEta = loop_absEta
    truthMatchCfg.loop_maxPt = loop_maxPt
    truthMatchCfg.loop_maxParticleHits = loop_maxParticleHits

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    matchAlg = TrackTruthMatcher(
        config=truthMatchCfg,
        level=customLogLevel(),
    )

    s.addAlgorithm(matchAlg)


def addTrackWriters(
    s: acts.examples.Sequencer,
    name: str,
    tracks: str = "tracks",
    inputTrackParticleMatching: str = "track_particle_matching",
    inputParticleTrackMatching: str = "particle_track_matching",
    inputParticleMeasurementsMap: str = "particle_measurements_map",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    writeSummary: bool = True,
    writeStates: bool = False,
    writeFitterPerformance: bool = False,
    writeFinderPerformance: bool = False,
    logLevel: Optional[acts.logging.Level] = None,
    writeCovMat: bool = False,
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir(parents=True, exist_ok=True)

        if writeSummary:
            trackSummaryWriter = acts.examples.root.RootTrackSummaryWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching=inputTrackParticleMatching,
                filePath=str(outputDirRoot / f"tracksummary_{name}.root"),
                treeName="tracksummary",
                writeCovMat=writeCovMat,
            )
            s.addWriter(trackSummaryWriter)

        if writeStates:
            trackStatesWriter = acts.examples.root.RootTrackStatesWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching=inputTrackParticleMatching,
                inputSimHits="simhits",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDirRoot / f"trackstates_{name}.root"),
                treeName="trackstates",
            )
            s.addWriter(trackStatesWriter)

        if writeFitterPerformance:

            #cfg = RootTrackFitterPerformanceWriter.Config()
            #cfg.inputTracks = tracks
            #cfg.inputParticles = "particles_selected"
            #cfg.inputTrackParticleMatching="track_particle_matching"
            #cfg.resPlotToolConfig = alice3_plotting.resPlotToolConfig
            #cfg.filePath=str(outputDirRoot / f"performance_fitting_{name}.root")

            #print("PF:: ResPlotToolConfig")
            #print(cfg.resPlotToolConfig.varBinning["Eta"].bins)
            
            #trackFitterPerformanceWriter = (
            #    RootTrackFitterPerformanceWriter(
            #        cfg,
            #        level=customLogLevel())
            #)

            trackFitterPerformanceWriter = acts.examples.root.RootTrackFitterPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching=inputTrackParticleMatching,
                resPlotToolConfig=alice3_plotting.resPlotToolConfig,
                effPlotToolConfig=alice3_plotting.effPlotToolConfig,
                trackSummaryPlotToolConfig=alice3_plotting.trackSummaryPlotToolConfig,
                filePath=str(outputDirRoot / f"performance_fitting_{name}.root"),
            )
            s.addWriter(trackFitterPerformanceWriter)

        if writeFinderPerformance:
            trackFinderPerfWriter = acts.examples.root.RootTrackFinderPerformanceWriter(
                level=customLogLevel(),
                inputTracks=tracks,
                inputParticles="particles_selected",
                inputTrackParticleMatching=inputTrackParticleMatching,
                inputParticleTrackMatching=inputParticleTrackMatching,
                inputParticleMeasurementsMap=inputParticleMeasurementsMap,
                effPlotToolConfig=alice3_plotting.effPlotToolConfig,
                fakePlotToolConfig=alice3_plotting.fakePlotToolConfig,
                duplicationPlotToolConfig=alice3_plotting.duplicationPlotToolConfig,
                trackQualityPlotToolConfig=alice3_plotting.trackQualityPlotToolConfig,
                trackSummaryPlotToolConfig=alice3_plotting.trackSummaryPlotToolConfig,
                filePath=str(outputDirRoot / f"performance_finding_{name}.root"),
            )
            s.addWriter(trackFinderPerfWriter)

def addMeasurementFilterAlgorithm(
    s: acts.examples.Sequencer,
    inputTracks: str,
    inputMeasurementSubset: str,
    outputMeasurementSubset: str,
    includeOutliers: bool = False,
    logLevel: acts.logging.Level = None,
):
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    s.addAlgorithm(
        acts.examples.MeasurementFilterAlgorithm(
            level=customLogLevel(),
            inputTracks=inputTracks,
            inputMeasurementSubset=inputMeasurementSubset,
            outputMeasurementSubset=outputMeasurementSubset,
            includeOutliers=includeOutliers,
        )
    )

def addTrackPerformanceWriters(
    sequence: acts.examples.Sequencer,
    outputDirRoot: Union[Path, str],
    tracks: str,
    prototracks: str,
    selectedParticles: str,
    inputParticles: str,
    outputTrackParameters: str,
    logLevel: acts.logging.Level = None,
):
    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)
    outputDirRoot = Path(outputDirRoot)
    if not outputDirRoot.exists():
        outputDirRoot.mkdir(parents=True, exist_ok=True)

    sequence.addWriter(
        acts.examples.root.RootTrackFinderPerformanceWriter(
            level=customLogLevel(),
            inputTracks=tracks,
            inputParticles=selectedParticles,
            inputTrackParticleMatching="seed_particle_matching",
            inputParticleTrackMatching="particle_seed_matching",
            inputParticleMeasurementsMap="particle_measurements_map",
            filePath=str(outputDirRoot / "performance_seeding.root"),
        )
    )


def addTrackMerger(
        sequence: acts.examples.Sequencer,
        inputTrackCollections,
        outputTrackCollection: str,
        logLevel: acts.logging.Level = None,):

    customLogLevel = acts.examples.defaultLogging(sequence, logLevel)

    cfgMerger = TrackMergerAlgorithm.Config()
    cfgMerger.inputTrackCollections = inputTrackCollections
    cfgMerger.outputTrackCollection = outputTrackCollection

    trackMergerAlg = TrackMergerAlgorithm(
        config=cfgMerger,
        level=customLogLevel(),
    )

    sequence.addAlgorithm(trackMergerAlg)


def addIterativeTracking(
    s: acts.examples.Sequencer = None,
    geo_dir: pathlib.Path = None,
    trackingGeometry: acts.TrackingGeometry = None,
    field: Literal[acts.MagneticFieldMapRz, acts.MagneticFieldMapXyz] = None,
    iterations: int = 0,
    inputTracks: str = "ckf_tracks",
    inputMeasurementSubset: str = "measurement_subset",
    outputDir: Optional[Union[Path, str]] = None,):

    seedTrackCollectionForMerging = ["seed-tracks"]
    trackCollectionForMerging = ["ambi_tracks"]

    mergedSeedTrackCollection = "seed-tracks-merged"
    mergedTrackCollection = "ambi-tracks-merged"

    previousTracks = inputTracks
    previousMeasurementSubset = inputMeasurementSubset

    for iteration in range(1, iterations):
        prefix = f"iter_{iteration}_"

        iterOutputDir = Path(outputDir) / f"iter_{iteration}"
        iterOutputDir.mkdir(parents=True, exist_ok=True)

        currentMeasurementSubset = f"{prefix}measurement_subset"
        currentSeedTracks = f"{prefix}seed-tracks"
        currentCkfTracks = f"{prefix}ckf_tracks"
        currentAmbiTracks = f"{prefix}ambi_tracks"

        print("Maria:: Creating MeasurementSubset for iteration", iteration)
        print("Maria:: prefix =", prefix)
        print("Maria:: previousTracks =", previousTracks)
        print("Maria:: previousMeasurementSubset =", previousMeasurementSubset)
        print("Maria:: currentMeasurementSubset =", currentMeasurementSubset)
        # Each iteration of tracking uses left over hits

        addMeasurementFilterAlgorithm(
            s,
            inputTracks=previousTracks,
            inputMeasurementSubset=previousMeasurementSubset,
            outputMeasurementSubset=currentMeasurementSubset,
            includeOutliers=False,
            logLevel=acts.logging.INFO,
        )

        alice3_seeding.addSeeding(
            s,
            trackingGeometry,
            field,
            geoSelectionConfigFile=geo_dir / "../seedingConfigurations" / cfg.seeding.seedingLayers,
            seedFinderConfigArg=alice3_seeding.get_seed_finder_config(iteration),
            seedFinderOptionsArg=alice3_seeding.DefaultSeedFinderOptionsArg,
            seedFilterConfigArg=alice3_seeding.DefaultSeedFilterConfigArg,
            spacePointGridConfigArg=alice3_seeding.DefaultSpacePointGridConfigArg,
            seedingAlgorithmConfigArg=alice3_seeding.DefaultSeedingAlgorithmConfigArg,
            outputDirRoot=iterOutputDir,
            initialSigmas=[
                1 * u.mm,
                1 * u.mm,
                1 * u.degree,
                1 * u.degree,
                0.1 * u.e / u.GeV,
                1 * u.ns,
            ],
            initialSigmaPtRel=0.1,
            initialVarInflation=alice3_seeding.PavelInitialVarInflation,
            particleHypothesis=acts.ParticleHypothesis.pion,
            inputMeasurements=currentMeasurementSubset,
            outputSpacePoints="spacepoints",
            prefix=prefix,
        )

        # run CKF
        print("PF:: Running CKF! at iteration", str(iteration))
        addCKFTracks(
            s,
            trackingGeometry,
            field,
            get_track_selector_config(iteration),
            #TrackSelectorConfig(pt=(0.05 * u.MeV, None),
            #                    nMeasurementsMin=cfg.tracking.nMeasurementsMin,
            #                    maxHoles=2,
            #                    maxOutliers=2,
            #                    maxSharedHits=2),
            CkfConfig(
                seedDeduplication=True,
                stayOnSeed=True,
                chi2CutOffMeasurement=cfg.tracking.ckfChi2Measurement,
                chi2CutOffOutlier=cfg.tracking.ckfChi2Outlier,
                numMeasurementsCutOff=cfg.tracking.ckfMeasPerSurf,
            ),
            twoWay=cfg.tracking.twoWayCKF,
            outputDirRoot=iterOutputDir,
            writeTrackSummary=cfg.tracking.writeTrackSummary,
            writeTrackStates=False,
            inputMeasurements=currentMeasurementSubset,
            inputMeasurementParticlesMap="measurement_particles_map",
            inputParticleMeasurementsMap="particle_measurements_map",
            inputInitialTrackParameters=f"{prefix}estimatedparameters",
            inputSeeds=f"{prefix}estimatedseeds",
            outputTracks=currentCkfTracks,
            logLevel=acts.logging.INFO,
        )

        s.addWhiteboardAlias(f"{prefix}tracks", currentCkfTracks)

        s = addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(
                maximumSharedHits=cfg.tracking.maxSharedHits,
                nMeasurementsMin=cfg.tracking.nMeasurementsMin,
            ),
            prefix=prefix,
            outputDirRoot=iterOutputDir,
            logLevel=acts.logging.INFO,
        )

        # Add the seed tracks for merging and the measurement mapping for this iteration
        seedTrackCollectionForMerging.append(currentSeedTracks)
        trackCollectionForMerging.append(currentAmbiTracks)

        previousTracks = currentAmbiTracks
        previousMeasurementSubset = currentMeasurementSubset

        print("Maria:: DEBUG Merging iteration")
        print(seedTrackCollectionForMerging)
        print(trackCollectionForMerging)


    # Merge seeds
    addTrackMerger(
        s,
        seedTrackCollectionForMerging,
        mergedSeedTrackCollection,
        acts.logging.INFO,
    )

    # Match seeds
    addTrackTruthMatcher(
        s,
        inputTracks=mergedSeedTrackCollection,
        inputParticles="particles_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="seed_merged_particle_matching",
        outputParticleTrackMatching="particle_seed_merged_matching",
        logLevel=acts.logging.INFO,
    )

    # Merge tracks
    addTrackMerger(
        s,
        trackCollectionForMerging,
        mergedTrackCollection,
        acts.logging.INFO,
    )

    # Match tracks
    addTrackTruthMatcher(
        s,
        inputTracks=mergedTrackCollection,
        inputParticles="particles_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ambi_tracks_merged_particle_matching",
        outputParticleTrackMatching="particle_ambi_tracks_merged_matching",
        logLevel=acts.logging.INFO,
    )

    s.addWhiteboardAlias("tracks", mergedTrackCollection)
    
    #   tracksummary_ckf-tracks-merged.root
    addTrackWriters(
        s,
        name=mergedTrackCollection,
        tracks=mergedTrackCollection,
        inputTrackParticleMatching="ambi_tracks_merged_particle_matching",
        inputParticleTrackMatching="particle_ambi_tracks_merged_matching",
        inputParticleMeasurementsMap="particle_measurements_map",
        outputDirCsv=None,
        outputDirRoot=outputDir,
        writeSummary=True,
        writeStates=False,
        writeFitterPerformance=False,
        writeFinderPerformance=False,
        writeCovMat=False,
        logLevel=acts.logging.INFO,
    )

    s.addWriter(
        acts.examples.root.RootTrackFinderPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks=mergedSeedTrackCollection,
            inputParticles="particles_selected",
            inputTrackParticleMatching="seed_merged_particle_matching",
            inputParticleTrackMatching="particle_seed_merged_matching",
            inputParticleMeasurementsMap="particle_measurements_map",
            effPlotToolConfig=alice3_plotting.effPlotToolConfig,
            fakePlotToolConfig=alice3_plotting.fakePlotToolConfig,
            filePath=str(Path(outputDir) / "performance_merged_seed.root"),
        )
    )

    s.addWriter(
        acts.examples.root.RootTrackFinderPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks=mergedTrackCollection,
            inputParticles="particles_selected",
            inputTrackParticleMatching="ambi_tracks_merged_particle_matching",
            inputParticleTrackMatching="particle_ambi_tracks_merged_matching",
            inputParticleMeasurementsMap="particle_measurements_map",
            effPlotToolConfig=alice3_plotting.effPlotToolConfig,
            fakePlotToolConfig=alice3_plotting.fakePlotToolConfig,
            filePath=str(Path(outputDir) / "performance_merged_ambi_tracks.root"),
        )
    )

    return s

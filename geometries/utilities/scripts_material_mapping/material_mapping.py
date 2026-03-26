#!/usr/bin/env python3

import os
import argparse

import acts
from acts import (
    SurfaceMaterialMapper,
    VolumeMaterialMapper,
    Navigator,
    Propagator,
    StraightLineStepper,
)
from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    MaterialMapping,
    )


from acts.examples.json import (
    JsonMaterialWriter,
    JsonFormat,
)

from acts.json import MaterialMapJsonConverter
from acts.examples.root import RootMaterialTrackReader
from acts.examples.root import RootMaterialTrackWriter

# from acts.examples.odd import getOpenDataDetector
#import alice3
import pathlib
import argparse
import sys

def runMaterialMapping(
    trackingGeometry,
    decorators,
    outputDir,
    inputDir,
    mapName="material-map",
    mapFormat=JsonFormat.Json,
    mapSurface=True,
    mapVolume=True,
    readCachedSurfaceInformation=False,
    mappingStep=1,
    s=None,
):
    s = s or Sequencer(numThreads=1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    wb = WhiteBoard(acts.logging.INFO)

    context = AlgorithmContext(0, 0, wb, 0)

    for decorator in decorators:
        assert decorator.decorate(context) == ProcessCode.SUCCESS

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks="material-tracks",
            fileList=[
                os.path.join(
                    inputDir,
                    (
                        mapName + "_tracks.root"
                        if readCachedSurfaceInformation
                        else "geant4_material_tracks.root"
                    ),
                )
            ],
            readCachedSurfaceInformation=readCachedSurfaceInformation,
        )
    )

    stepper = StraightLineStepper()

    mmAlgCfg = MaterialMapping.Config(context.geoContext, context.magFieldContext)
    mmAlgCfg.trackingGeometry = trackingGeometry
    mmAlgCfg.inputMaterialTracks = "material-tracks"

    if mapSurface:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
            resolveSensitive=True,
            resolveMaterial=True,
            resolvePassive=True,
        )
        propagator = Propagator(stepper, navigator)
        mapper = SurfaceMaterialMapper(level=acts.logging.INFO, propagator=propagator)
        mmAlgCfg.materialSurfaceMapper = mapper

    if mapVolume:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
        )
        propagator = Propagator(stepper, navigator)
        mapper = VolumeMaterialMapper(
            level=acts.logging.INFO, propagator=propagator, mappingStep=mappingStep
        )
        mmAlgCfg.materialVolumeMapper = mapper

    jmConverterCfg = MaterialMapJsonConverter.Config(
        processSensitives=True,
        processApproaches=True,
        processRepresenting=True,
        processBoundaries=True,
        processVolumes=True,
        context=context.geoContext,
    )

    jmw = JsonMaterialWriter(
        level=acts.logging.INFO,
        converterCfg=jmConverterCfg,
        fileName=os.path.join(outputDir, mapName),
        writeFormat=mapFormat,
    )

    mmAlgCfg.materialWriters = [jmw]

    s.addAlgorithm(MaterialMapping(level=acts.logging.INFO, config=mmAlgCfg))

    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=mmAlgCfg.mappingMaterialCollection,
            filePath=os.path.join(
                outputDir,
                mapName + "_tracks.root",
            ),
            storeSurface=True,
            storeVolume=True,
        )
    )

    return s

parser = argparse.ArgumentParser(
    description="Geometry to json")

parser.add_argument(
    "--geo-dir",
    help="Location of the geometry",
    type=pathlib.Path,
    default=pathlib.Path.cwd()
    )

args = parser.parse_args()

            
geo_dir = args.geo_dir
print("Loading geometry from ... ",str(geo_dir))
sys.path.insert(0,str(pathlib.Path(geo_dir).resolve()))

import buildALICE3Geometry


if "__main__" == __name__:
    matDeco = acts.IMaterialDecorator.fromFile("geometry-map.json")

    p = argparse.ArgumentParser()
    
    
    detector = buildALICE3Geometry.buildALICE3Geometry(
       geo_dir, False, False, acts.logging.VERBOSE, matDeco)
        # geo_dir, True, False, acts.logging.INFO)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    mapFormat = JsonFormat.Json
    #mapName = args.outFile.split(".")[0]

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=os.getcwd(),
        inputDir=os.getcwd(),
        readCachedSurfaceInformation=False,
        #mapName=mapName,
        mapFormat=mapFormat,
    ).run()

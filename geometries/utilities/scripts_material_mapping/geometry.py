#!/usr/bin/env python3
import os
#from common import getOpenDataDetectorDirectory
#from acts.examples.odd import getOpenDataDetector
from acts.examples import (
    GenericDetector,
#    AlignedDetector,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
)

from acts.examples.json import (    
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)

import acts
from acts.json import MaterialMapJsonConverter

#import alice3
import pathlib
import argparse
import sys


def runGeometry(
    trackingGeometry,
    decorators,
    outputDir,
    events=1,
    outputObj=False,
    outputCsv=False,
    outputJson=True,
    outputRoot=False,
):

    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0
        ithread = 0

        context = AlgorithmContext(ialg, ievt, eventStore, ithread)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=os.path.join(outputDir, "csv"),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO, outputDir=os.path.join(outputDir, "obj")
            )
            writer.write(context, trackingGeometry)

        if outputJson:
            writer = JsonSurfacesWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=os.path.join(outputDir, "json"),
                writePerEvent=True,
                writeSensitive=True,
            )
            writer.write(context)

            jmConverterCfg = MaterialMapJsonConverter.Config(
                processSensitives=True,
                processApproaches=True,
                processRepresenting=True,
                processBoundaries=True,
                processVolumes=True,
                processNonMaterial=True,
                context=context.geoContext,
            )

            jmw = JsonMaterialWriter(
                level=acts.logging.VERBOSE,
                converterCfg=jmConverterCfg,
                fileName=os.path.join(outputDir, "geometry-map"),
                writeFormat=JsonFormat.Json,
            )

            jmw.write(trackingGeometry)



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
#    detector, trackingGeometry, decorators = AlignedDetector.create()

#     detector, trackingGeometry, decorators = alice3.buildALICE3Geometry(
# #    geo_dir, False, False, acts.logging.VERBOSE)
# #    geo_dir, False, False, acts.logging.DEBUG)
#     geo_dir, False, False, acts.logging.INFO)
#     # detector, trackingGeometry, decorators = GenericDetector.create()
#     # detector, trackingGeometry, decorators = getOpenDataDetector(getOpenDataDetectorDirectory())

    detector = buildALICE3Geometry.buildALICE3Geometry(
    #    geo_dir, False, False, acts.logging.VERBOSE)
        geo_dir, False, False, acts.logging.INFO)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()


    runGeometry(trackingGeometry, decorators, outputDir=os.getcwd())

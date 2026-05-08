#!/bin/bash

# Default value
GEOMETRY="invalid"
NEVENTS=4000

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --geometry)
            GEOMETRY="$2"
            shift 2
            ;;
        -n|--events)
            NEVENTS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Use the geometry variable
echo "Running with geometry: $GEOMETRY"
echo "Number of events: $NEVENTS"


# visualization of ACTS tracking geometry
python geometry.py --geo-dir $GEOMETRY
python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/GeometryVisualisationAndMaterialHandling.py --geometry geometry-map.json

# record material hits via geantino scan
OUTPUT_GEANTINO_HITS="geant4_material_tracks.root"
python material_recording.py --input ${GEOMETRY}/o2sim_geometry.gdml -n $NEVENTS \
    --output ${OUTPUT_GEANTINO_HITS}

# tune binning, map material
python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/writeMapConfig.py geometry-map.json config-map.json
python modifyConfigMap.py config-map.json 80 80
python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/configureMap.py geometry-map.json config-map-binned.json
MATERIAL_MAP_OUTPUT_BASE="material"
python material_mapping.py --geo-dir ${GEOMETRY} -n $NEVENTS \
    --input ${OUTPUT_GEANTINO_HITS} --output ${MATERIAL_MAP_OUTPUT_BASE}

# validation
python material_validation.py --geo-dir ${GEOMETRY} -n $NEVENTS \
    --map "${MATERIAL_MAP_OUTPUT_BASE}_map.json"
mkdir -p Validation
root -l -b -q "$ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/Mat_map.C(\
    \"propagation-material.root\",\
    \"${MATERIAL_MAP_OUTPUT_BASE}_mapped.root\",\
    \"Validation\")"

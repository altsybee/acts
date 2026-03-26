#!/bin/bash

# Default value
GEOMETRY="invalid"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --geometry)
            GEOMETRY="$2"
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


# visualization of ACTS tracking geometry
python geometry.py --geo-dir $GEOMETRY
python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/GeometryVisualisationAndMaterialHandling.py --geometry geometry-map.json

# record material hits via geantino scan
python material_recording.py --gdml ${GEOMETRY}/o2sim_geometry.gdml

# tune binning, map material
python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/writeMapConfig.py geometry-map.json config-map.json
python modifyConfigMap.py config-map.json 80 80
python $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/configureMap.py geometry-map.json config-map-binned.json
python material_mapping.py --geo-dir ${GEOMETRY}

# validation
python material_validation.py --geo-dir ${GEOMETRY}
mkdir Validation
root -l -b -q $ACTS_SOURCE_DIR/Examples/Scripts/MaterialMapping/Mat_map.C'("propagation-material.root","material-map_tracks.root","Validation")'

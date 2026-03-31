#!/usr/bin/env bash

# Source the ACTS setup script to set up the environment
. $ACTS_ROOT/python/setup.sh

# Absolute path to this setup.sh
export MAIN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

# Add the directory to the python path
export PYTHONPATH="${MAIN_DIR}:${PYTHONPATH}"

# Add the local python bindings to the python path
if [ -z "$ACTSO2_ROOT" ]; then
    echo "ACTSO2_ROOT is not set. Using local build directory for Python bindings."
    export PYTHONPATH="${MAIN_DIR}/build/ActsAlgorithms/python:${PYTHONPATH}"
else
    echo "Using ACTSO2_ROOT=${ACTSO2_ROOT} for Python bindings."
    export PYTHONPATH="${ACTSO2_ROOT}/python:${PYTHONPATH}"
fi

# Disable FPE errors
export ACTS_SEQUENCER_DISABLE_FPEMON=1

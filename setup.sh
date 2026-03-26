#!/usr/bin/env bash

# Absolute path to this setup.sh
export MAIN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

# Add the directory to the python path
export PYTHONPATH="${MAIN_DIR}:${PYTHONPATH}"

# Add the local python bindings to the python path
export PYTHONPATH="${MAIN_DIR}/build/ActsAlgorithms/python:${PYTHONPATH}"

# Disable FPE errors
export ACTS_SEQUENCER_DISABLE_FPEMON=1

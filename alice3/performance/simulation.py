from typing import Optional, Union, List, Literal, Tuple
from pathlib import Path

import acts
import acts.examples
import acts.examples.reconstruction as acts_reco
from acts import UnitConstants as u

from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addGenParticleSelection,
    addSimParticleSelection,
    ParticleSelectorConfig,
)


#Load config
from job_configs import ChainConfig
cfg = ChainConfig.Config()


def addDigiParticleSelection(
        s: acts.examples.Sequencer,
        config: ParticleSelectorConfig,
        measurementLayers: List[Tuple[int,int]] = [], 
        logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """
    This function steers the particle selection after digitization.

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the ParticleSelector
    config: ParticleSelectorConfig
        the particle selection configuration
    measurementLayers: list of (VolID,lyID) to where adding the measurement counter. 
    at least 1 measurement will be requested
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)
    
    # Get filtered kwargs
    kwargs = acts.examples.defaultKWArgs(**acts.examples.simulation._getParticleSelectionKWargs(config))

    # Create config object
    cfg = acts.examples.ParticleSelector.Config()

    # Apply all kwargs to config
    for key, value in kwargs.items():
        setattr(cfg, key, value)

    # Override specific fields
    cfg.inputParticles = "particles_simulated_selected"
    cfg.inputParticleMeasurementsMap = "particle_measurements_map"
    cfg.inputMeasurements = "measurements"
    cfg.outputParticles = "tmp_particles_digitized_selected"

    
    geoIDs = []
    for tup in measurementLayers:
        geoIDs.append(acts.GeometryIdentifier(volume=tup[0],layer=tup[1]))


    print("PF:: geoIDs=",geoIDs)
        
    if len(measurementLayers) > 0:
        cfg.measurementCounter.addCounter(geoIDs,1,1)
    
    
    selector = acts.examples.ParticleSelector(
        config=cfg,
        level=customLogLevel(),
    )
    
    s.addAlgorithm(selector)

    s.addWhiteboardAlias("particles_selected", selector.config.outputParticles)
    s.addWhiteboardAlias(
        "particles_digitized_selected", selector.config.outputParticles
    )


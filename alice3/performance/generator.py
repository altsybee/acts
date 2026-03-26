from typing import Optional, Union, List, Literal
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


def addAlice3ParticleGun(s : acts.examples.Sequencer = None,
                         rnd : acts.examples.RandomNumbers = None,
                         outputDir: Optional[Union[Path, str]] = None):

    addParticleGun(
        s,
        MomentumConfig(
            cfg.particleGun.gunPtRange[0] * u.GeV, cfg.particleGun.gunPtRange[1] * u.GeV, transverse=True
        ),
        EtaConfig(cfg.particleGun.gunEtaRange[0], cfg.particleGun.gunEtaRange[1], uniform=True),
        PhiConfig(cfg.particleGun.gunPhiRange[0], cfg.particleGun.gunPhiRange[1]),
        ParticleConfig(
            cfg.particleGun.gunMult,
            acts.PdgParticle(cfg.particleGun.gunPID),
            randomizeCharge=cfg.particleGun.gunRandCharge,
        ),
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(
                0.0001 * u.mm, 0.0001 * u.mm, 0.0001 * u.mm, 0.0001 * u.ns
            ),
            mean=acts.Vector4(0, 0, cfg.particleGun.meanVz * u.mm, 0),
        ),
        rnd=rnd,
        logLevel=acts.logging.INFO,
        outputDirRoot=outputDir,
)


def addAlice3Pythia8(s : acts.examples.Sequencer = None,
                     rnd : acts.examples.RandomNumbers = None,
                     outputDir: Optional[Union[Path, str]] = None):

    beam=acts.PdgParticle.eProton
    hardProcess=[]
    cmsEnergy=0.
    if cfg.pythia.system == "pp":
        beam=acts.PdgParticle.eProton,
        cmsEnergy=13.6 * acts.UnitConstants.TeV
        hardProcess=[
            "SoftQCD:inelastic = on",
            "Tune:pp = 14",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10",
        ]  # tune 14 is Monash
        pileupProcess=[
            "SoftQCD:inelastic = on",
            "Tune:pp = 14",
            "ParticleDecays:limitTau0 = on",
            "ParticleDecays:tau0Max = 10",
        ]  # tune 14 is Monash
    elif cfg.pythia.system == "PbPbMB":
        beam=acts.PdgParticle.eLead,
        cmsEnergy=5.36 * acts.UnitConstants.TeV
        hardProcess=[
                "SoftQCD:inelastic = on",
                "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
                "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
                "HeavyIon:SigFitNGen = 20",
                "ParticleDecays:limitTau0 = on",
                "ParticleDecays:tau0Max = 10",
            ]
        pileupProcess=[ # MB as pileup process
                "SoftQCD:inelastic = on",
                "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
                "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
                "HeavyIon:SigFitNGen = 20",
                "ParticleDecays:limitTau0 = on",
                "ParticleDecays:tau0Max = 10",
        ]
    elif cfg.pythia.system == "PbPbCentral":
        beam=acts.PdgParticle.eLead
        cmsEnergy=5.36 * acts.UnitConstants.TeV
        hardProcess=[
                "SoftQCD:inelastic = on",
                "HeavyIon:bWidth=0.1",
                "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
                "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
                "HeavyIon:SigFitNGen = 20",
                "ParticleDecays:limitTau0 = on",
                "ParticleDecays:tau0Max = 10",
        ]
        pileupProcess=[ # MB as pileup process
                "SoftQCD:inelastic = on",
                "HeavyIon:SigFitErr =  0.02,0.02,0.1,0.05,0.05,0.0,0.1,0.0",
                "HeavyIon:SigFitDefPar = 17.24,2.15,0.33,0.0,0.0,0.0,0.0,0.0",
                "HeavyIon:SigFitNGen = 20",
                "ParticleDecays:limitTau0 = on",
                "ParticleDecays:tau0Max = 10",
        ]
    else:
        print("Pythia process ", cfg.pythia.system, " not available!")



    s = addPythia8(
        s,
        npileup=cfg.pythia.pileup,
        beam = beam[0],
        cmsEnergy = cmsEnergy,
        hardProcess = hardProcess,
        pileupProcess = pileupProcess,
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(
                0.000125 * u.mm,
                0.000125 * u.mm,
                cfg.pythia.sigmaVz * u.mm,
                0.0001 * u.ns,  #!! makes 0 sense
                ),
            mean=acts.Vector4(0, 0, cfg.pythia.meanVz * u.mm, 0),
            ),
        rnd=rnd,
        logLevel=acts.logging.INFO,
        outputDirRoot=outputDir)
    
    return s




**In is the v3 version of the geometry, created on Feb 13, 2026,** <br>
where the barrel ML and OT layers are **cylindrical**, VD is also **purely cylindrical**.<br>
[ please unzip material-map.json.zip, ~130 Mb ]

The O2 command used: `o2-sim-serial-run5 -n 1 -g pythia8hi -m A3IP TRK TF3 --configKeyValues "TRKBase.layoutML=kCylinder;TRKBase.layoutOL=kCylinder;"`<br>
(please note: no FT3)

The full "post-processing" on top of the O2 `.root` geometry is done via `bash run_all_modifications.sh` script stored 
[here.](https://gitlab.cern.ch/alice3-trackers/wp1-simulationsandperformances/acts/geometries/-/tree/master/utilities/scripts_geometry_modification) <br>
In particular, the ML and OT disks are created and segmented using a hand-made procedure, in order to mockup a realistic detector and make sim+reco runnable in full $\eta$ acceptance with ACTS (Fatras, Geant 4). More details in: [WP1 report Jan 29, 2026](https://indico.cern.ch/event/1642644/#13-updates-on-geometry)


After this, the material giantino scan was performed, and material mapping on surfaces was done (the usual procedure).

Some QA plots are in the `plots` folder.

**(!) Note** that in `digi*.json` files in this folder **the iTOF+oTOF (barrel) layers are considered as sensitive** -- so there will be measurements on these layers in the output!
* Use `digi-smearing-config_with_TOFs_noEndcapTOFs_noTimeInTOFs.json` if no time coordinate is needed (3D tracking).  
* In `digi-smearing-config_with_TOFs_tunedSmearingForTOF_noEndcapTOFs.json`, the time resolution is taken into account for TOFs --> experiments with the 4D tracking can be done.



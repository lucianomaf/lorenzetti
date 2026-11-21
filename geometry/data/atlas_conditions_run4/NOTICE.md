# ATLAS Run 4 conditions used by Lorenzetti

The files in this folder are derived from **data of the ATLAS Collaboration**, so that Lorenzetti can use them without access to any ATLAS
service. They are **not covered by the GPL license of Lorenzetti**. Copyright CERN for the benefit of the ATLAS Collaboration. The sources
are reachable on the internet **without any account** and do not state a license; if you redistribute these files, keep this notice and cite
the sources. They describe the **ATLAS simulation and calibration**, not the detector as built for Run 4.

**Every number can be rebuilt from public sources** with the scripts of `geometry/scripts/atlas_conditions_run4/` (see its `README.md`):
the ATLAS COOL conditions database through the public Frontier server, the CernVM-FS repository `atlas-condb.cern.ch` (free client, no
account), the CREST service and the open Athena code (gitlab.cern.ch/atlas/athena, branch `release/25.0.21`, readable without login).

**One value per layer and |eta|**, the same on both sides: Lorenzetti uses no per-channel constants. EMB, EMEC and HEC: one ring per
(wheel, sampling, region, eta index) of the ATLAS cell identifier, at the nominal |eta| of the Lorenzetti cell (`abs_eta_min` =
`abs_eta_max`); FCal: bands of |eta| 0.1 wide per module (`ring` = 1000 + floor(|eta| / 0.1)). A cell takes the row of its layer and |eta|.
Constants: median over the cells of the ring, with the number of cells more than 5% away from it in `<quantity>_n_off5pct`.

## `lar_run4_constants_by_ring.csv` (1003 rings)

| Column | Meaning | Source |
|---|---|---|
| `fsampl_mc` | sampling fraction of the simulation | COOL `/LAR/ElecCalibMC/fSampl`, global tag `OFLCOND-MC21-SDR-RUN4-02`, run 350000 |
| `adc2mev_hg_mc` | MeV per ADC count, high gain, of the official Run 4 simulation with the COOL global tag (the Run 3 electronics): ramp × DAC2uA × uA2MeV / Mphys/Mcal × HVScaleCorr, as Athena 25.0.21 combines them for simulation (`LArRecUtils/python/LArADC2MeVCondAlgConfig.py`) | COOL `/LAR/ElecCalibMC/{Ramp, DAC2uA, uA2MeV, MphysOverMcal, HVScaleCorr}`, same global tag and run |
| `ramp_hg_mc` | ramp (c1) of the high gain of that simulation | COOL `/LAR/ElecCalibMC/Ramp` |
| `ramp_r4_hg`, `ramp_r4_lg` | ramp (c1; c0 is zero) of the Phase-II front-end board (FEB 2), high and low gain | CREST, global tag `COND-MC21-SDR-RUN4-06`, tag `LARElecCalibMCRamp-R4-00` |
| `noise_adc_r4_hg`, `noise_adc_r4_lg` | noise of the FEB 2 in ADC counts | CREST, tag `LARElecCalibMCNoise-R4-00` |

The FEB 2 pedestal is one value for all channels: 8192 ADC counts in both gains (CREST tag `LARElecCalibMCPedestal-R4-00`). MeV per ADC count
of the FEB 2: `adc2mev_hg_mc * ramp_r4 / ramp_hg_mc`. CREST marks `COND-MC21-SDR-RUN4-06` as scenario "test" ("Migration global tag"):
these values may change. The CREST payloads are ordered by the LAr online hash; the link to the cells is the simulation cabling
`/LAR/Identifier/OnOffIdMap`, tag `LARIdentifierOnOffIdMap-012` (the tag Athena 25.0.21 sets for simulation,
`LArCabling/python/LArCablingConfig.py`).

## `lar_run4_noise_autocorr_by_ring.csv`

The correlation of the LAr noise between samples 1 to 4 apart, per ring, for the high and the medium gain (the medium one is the one to use
for the Run 4 low gain): COOL `/LAR/ElecCalibMC/AutoCorr`, same global tag and run. The simulation stores one vector per ring; in the FCal
bands, the mean of the vectors of the band. Beyond 4 samples the correlation is not measured.

## `lar_run4_pulses_by_ring.csv`

The LAr pulse per ring and gain (high and medium), 0.5 ns apart, peak 1 at t = 0, total area zero. Source: the physics waves of the Response
Transformation Method, one per channel, COOL `/LAR/ElecCalibOfl/PhysWaves/RTM`, tag `LARElecCalibOflPhysWavesRTM-RUN2-UPD3-00`, run 520451
(data global tag `CONDBR2-BLKPA-2026-01`, cabling tag `LARIdentifierOnOffIdMap-RUN2-001`). These are the pulses of the **present** LAr
electronics: no per-channel pulse of the Phase-II electronics is public. The FCal waves in that folder are 156 ns long; after them the pulse
continues with straight lines through the 32 samples at 25 ns of the `LArPileupShape` of the channel (COOL
`/LAR/ElecCalibOfl/LArPileupShape`, tag `LARElecCalibOflLArPileupShape-RUN2-2018`), placed and scaled by least squares on the measured wave.
How the pulse of a ring is made: each wave normalised to its peak with the peak at t = 0, resampled to 0.5 ns, the mean over the cells, a
constant subtracted so that the area is zero, and the peak brought back to 1; left out are waves still above 5% of their peak at their last
point and cells whose peak time is more than 10 ns from the median of the ring.

## `tile_run4_sampling_fraction_by_ring.csv` (41 rings)

Section, sampling and tower (both sides and both PMTs together): the median over the PMTs of the sampling fraction of the ATLAS simulation
of Run 4, COOL `/TILE/OFL02/CALIB/SFR`, tag `TileOfl02CalibSfr-SIM-07` (the tag Athena 25.0.21 forces for the Run 4 geometry), run 350000,
with the drawer defaults and the channel to cell map of the open Athena code (`TileCalibUtils.cxx`, `TileCablingService.cxx`). The source keeps
the Run 1 arrangement in 32 drawers; a few PMTs have values that come from that arrangement, and the median of a ring is not affected by them.

# ATLAS Run 4 conditions used by Lorenzetti

The files in this folder are **data of the ATLAS Collaboration**, copied here so that Lorenzetti can use them without access to any
ATLAS service. They are **not covered by the GPL license of Lorenzetti**. Copyright CERN for the benefit of the ATLAS Collaboration. The
sources below are reachable without authentication and do not state a license; if you redistribute these files, keep this notice and
cite the sources.

All values were read on 20 November 2026 (CREST) and on 1 and 18 November 2026 (COOL). They describe the **ATLAS simulation** of
Run 4, not the detector as built or measured.

## `lar_run4_constants.csv`

One line per LAr cell (182 468 cells).

| Column | Meaning | Source |
|---|---|---|
| `offline_id` | ATLAS offline identifier of the cell (64 bits) | ATLAS identifier dictionary (Athena 25.0.21) |
| `online_hash` | hash of the LAr online channel; index of the CREST payloads | ATLAS cabling, `/LAR/Identifier/OnOffIdMap` (COOL, global tag `OFLCOND-MC21-SDR-RUN4-02`, run 350000), read with Athena 25.0.21 |
| `layer`, `eta` | layer (EMB_0 to FCAL_3) and eta of the cell centre | ATLAS Run 4 geometry `ATLAS-P2-RUN4-03-00-02`, Athena 25.0.21 |
| `fsampl_mc` | sampling fraction of the simulation | COOL `/LAR/ElecCalibMC/fSampl`, `OFLCOND-MC21-SDR-RUN4-02`, run 350000 |
| `adc2mev_hg_mc`, `ramp_hg_mc` | MeV per ADC count and ramp (c1) of the high gain in the official Run 4 simulation with the COOL global tag (the Run 3 electronics) | COOL `/LAR/ElecCalibMC/Ramp` and the Athena `LArADC2MeV` conditions, same tag and run |
| `ramp_r4_hg`, `ramp_r4_lg` | ramp (c1; c0 is zero) of the FEB 2, high and low gain | CREST, global tag `COND-MC21-SDR-RUN4-06`, tag `LARElecCalibMCRamp-R4-00` ("Electronic gain of FEB 2") |
| `noise_adc_r4_hg`, `noise_adc_r4_lg` | noise of the FEB 2 in ADC counts, high and low gain | CREST, tag `LARElecCalibMCNoise-R4-00` ("Noise of FEB 2 (in ADC counts)") |

The FEB 2 pedestal is one value for all channels: 8192 ADC counts in both gains (second field of the payload: 17), CREST tag
`LARElecCalibMCPedestal-R4-00` ("Pedestal of FEB 2").

MeV per ADC count of the FEB 2, as the ATLAS digitisation combines them: `adc2mev_hg_mc * ramp_r4 / ramp_hg_mc` (the other factors of the
conversion, DAC to µA, µA to MeV and Mphys/Mcal, are the same in the two global tags).

CREST marks the global tag `COND-MC21-SDR-RUN4-06` as scenario "test" and describes it as "Migration global tag": these values may change.

`geometry/scripts/atlas_conditions_run4.py check` compares the CREST columns with CREST; `rebuild` rewrites them. The COOL columns
were read with Athena on lxplus and are not rebuilt by the script.

## `tile_run4_sampling_fraction.csv`

One line per Tile PMT (9852 PMTs): the PMT identifier (ATLAS offline identifier, 64 bits), its section, side, module, tower, sampling and
PMT number, and the sampling fraction of the ATLAS simulation of Run 4: COOL `/TILE/OFL02/CALIB/SFR`, tag `TileOfl02CalibSfr-SIM-07` (the
tag Athena 25.0.21 forces for the Run 4 geometry), read at run 350000, the run of the Run 4 simulation.

Known defect of the source: the table keeps the Run 1 arrangement in 32 drawers (with the MBTS in channel 0). 32 PMTs have values (1.0
and 37.0) that, as far as we could tell, come from that arrangement and not from a physical cell.

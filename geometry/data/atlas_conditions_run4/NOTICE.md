# ATLAS Run 4 conditions used by Lorenzetti

The files in this folder are derived from **data of the ATLAS Collaboration**, copied here so that Lorenzetti can use them without access
to any ATLAS service. They are **not covered by the GPL license of Lorenzetti**. Copyright CERN for the benefit of the ATLAS Collaboration.
The sources below are reachable without authentication and do not state a license; if you redistribute these files, keep this notice and
cite the sources.

They describe the **ATLAS simulation and calibration**, not the detector as built for Run 4. Values were read on 20 November 2026 (CREST)
and on 1 and 18 November 2026 (COOL).

**One value per layer and |eta|.** Lorenzetti does not use per-channel constants: every quantity is given per layer and ring of |eta|, the
same on the A and C sides, as the median (constants) or the mean (pulses) over the cells of the ring, all phi and both sides together. A ring
is the eta ring of the ATLAS cell table for EMB, EMEC and HEC; for the FCal, whose cells are not in eta rings, a ring is a band of |eta| 0.1
wide (ring = 1000 + floor(|eta| / 0.1)). A cell takes the values of the row of its layer whose `abs_eta_min` to `abs_eta_max` contains its
|eta|.

## `lar_run4_constants_by_ring.csv`

1003 rings. For each quantity, the median over the cells of the ring and, in the column `<quantity>_n_off5pct`, how many cells of the ring
are more than 5% away from that median.

| Column | Meaning | Source |
|---|---|---|
| `fsampl_mc` | sampling fraction of the simulation | COOL `/LAR/ElecCalibMC/fSampl`, global tag `OFLCOND-MC21-SDR-RUN4-02`, run 350000 |
| `adc2mev_hg_mc`, `ramp_hg_mc` | MeV per ADC count and ramp (c1) of the high gain in the official Run 4 simulation with the COOL global tag (the Run 3 electronics) | COOL `/LAR/ElecCalibMC/Ramp` and the Athena `LArADC2MeV` conditions, same tag and run |
| `ramp_r4_hg`, `ramp_r4_lg` | ramp (c1; c0 is zero) of the Phase-II front-end board (FEB 2), high and low gain | CREST (`https://crest.cern.ch/api-v6.0`), global tag `COND-MC21-SDR-RUN4-06`, tag `LARElecCalibMCRamp-R4-00` |
| `noise_adc_r4_hg`, `noise_adc_r4_lg` | noise of the FEB 2 in ADC counts | CREST, tag `LARElecCalibMCNoise-R4-00` |

The FEB 2 pedestal is one value for all channels: 8192 ADC counts in both gains (second field of the payload: 17), CREST tag
`LARElecCalibMCPedestal-R4-00`. MeV per ADC count of the FEB 2, as the ATLAS digitisation combines them: `adc2mev_hg_mc * ramp_r4 /
ramp_hg_mc`. `adc2mev_hg_mc` is the quantity that varies most inside a ring (15 495 cells more than 5% from the median of their ring); the
other quantities are constant in almost every ring. CREST marks `COND-MC21-SDR-RUN4-06` as scenario "test" ("Migration global tag"): these
values may change.

## `lar_run4_pulses_by_ring.csv`

The LAr pulse, one per ring and gain (`high`, and `medium`, which is the one to use for the low gain of Run 4), 1 ns apart, peak 1 at t = 0,
total area zero. Columns: layer, ring, |eta| range, gain, number of cells of the ring and of cells used, mean and spread of the peak time
of the cells, the time of the first point (`t_start_ns`, from the peak), the step (`dt_ns`), the number of points and the pulse (space
separated).

Source: the physics waves of the Response Transformation Method, one per channel, COOL `/LAR/ElecCalibOfl/PhysWaves/RTM`, tag
`LARElecCalibOflPhysWavesRTM-RUN2-UPD3-00`, run 520451, about 1 ns apart (the payload files are on the public CVMFS repository
`atlas-condb.cern.ch`). These are the pulses of the **present** LAr electronics, not of the Phase-II electronics, for which no per-channel
pulse is public. The FCal waves in that folder are 156 ns long; after their last point the pulse continues with straight lines through
the 32 samples at 25 ns of the `LArPileupShape` of the channel, scaled to the wave by a least-squares fit, and is zero one sample after the
last one.

How the pulse of a ring was made: each wave of the ring normalised to its peak, the peaks aligned at t = 0, resampled to 1 ns by linear
interpolation, then the mean over the cells, a constant subtracted over the window so that the area is zero, and the peak brought back to 1.
Left out of the mean: cells without a wave, waves still above 5% of their peak at their last point, and cells whose peak time is more than
10 ns from the median of the ring.

## `tile_run4_sampling_fraction_by_ring.csv`

43 rings (section, sampling and tower; both sides and both PMTs of a cell together): the median of the sampling fraction of the ATLAS
simulation of Run 4 over the PMTs of the ring, and how many PMTs are more than 5% from it. Source: COOL `/TILE/OFL02/CALIB/SFR`, tag
`TileOfl02CalibSfr-SIM-07` (the tag Athena 25.0.21 forces for the Run 4 geometry), read at run 350000. The source keeps the Run 1 arrangement
in 32 drawers; 32 PMTs have values (1.0 and 37.0) that, as far as we could tell, come from that arrangement and not from a physical cell.
The median of the ring is not affected by them.

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

## `tile_run4_noise_by_ring.csv` (41 rings × 2 gains)

The Tile noise per sample, converted to the counts of the Phase-II front-end (FENICS), per ring (section, sampling, tower; both sides, every
phi, both PMTs) and gain. Source: the sample noise of the ATLAS Run 4 simulation, COOL `/TILE/OFL02/NOISE/SAMPLE`, tag
`TileOfl02NoiseSample-TwoGauss-24` (from the global tag `OFLCOND-MC21-SDR-RUN4-02`), run 350000. As in the ATLAS digitisation
(`TileSimAlgs/src/TileDigitsMaker.cxx`), the noise of a sample is a double Gaussian: `sigma1_counts` with probability `w_sigma1`,
`sigma2_counts` otherwise. The simulation values are those of the present electronics (10-bit ADC, high gain up to 12.5 pC, low gain up to
800 pC; ATLAS-TDR-028, pp. 95-96) and are kept in the `*_legacy_counts` columns. Conversion, keeping the charge: the quantisation of the
10-bit ADC (1/sqrt(12) count) is taken out of each sigma in quadrature, and the sigma is multiplied by the charge per count of the present
electronics over that of the FENICS (12-bit ADC, ATLAS-TDR-028, p. 96; gain ratio 40 and low gain up to 1050 pC, PoS(ICHEP2020)860, p. 4;
so the high gain goes up to 26.25 pC): × (12.5/1023)/(26.25/4095) = 1.905 in high gain, × (800/1023)/(1050/4095) = 3.050 in low gain.
Each parameter is the median over the PMTs of the ring, without the dead channels (`n_dead`).

## `hllhc_bunch_train_ip1.csv`

The bunch-crossing slots of one LHC orbit (3564) at the ATLAS interaction point, as runs of colliding (1) and empty (0) slots: 2748
colliding bunch crossings. This file is not ATLAS data: it is derived from the HL-LHC baseline 25 ns filling scheme of the standard
production mode, `25ns_2760b_2748_2492_2574_288bpi_13inj_800ns_bs200ns` (2760 bunches per beam, 2748 collisions in IP1 and IP5; HL-LHC
Technical Design Report, CERN-2020-010, Table A-3, p. 347), as published in the examples of the FillingPatterns package
(https://github.com/giadarol/FillingPatterns). That repository declares no license, so its file is not copied here; `bunch_train.py`
downloads it and checks that it gives this table.

## Numbers of the electronics and of the beam (no table)

| Quantity | Value | Source |
|---|---|---|
| Tile ADC | 12 bits, 40 MHz | ATLAS-TDR-028, p. 96 |
| Tile gains | ratio 40 (gains 0.4 and 16), low gain up to 1050 pC | PoS(ICHEP2020)860, p. 4 |
| Tile scale | MeV per count of the ATLAS Run 4 simulation × 0.525 = (26.25 pC / 12.5 pC) × (1023 / 4095) | present high gain up to 12.5 pC = 800 pC / 64 (ATLAS-TDR-028, pp. 95-96); FENICS high gain 1050 pC / 40 |
| Tile pedestal | 90 counts, both gains | twice the median pedestal of the simulation (45 counts at 10 bits, `pedestal_legacy_counts`) |
| Tile pulse | CR-RC⁴, h(t) = (t/nτ)ⁿ exp(n(1 − t/nτ)), n = 4, τ = 10 ns, peak 1 at t = nτ | the standard CR-RCⁿ shaper (e.g. Electronics 14 (2025) 493, doi:10.3390/electronics14030493); FWHM 47.2 ns, close to the 50 ns of the FENICS (PoS(ICHEP2024)885, p. 4) |
| Tile long tail (optional, off by default) | first-order AC coupling of the pulse, total area zero, τ = 50 µs by default | a choice of the model, not a measured value |
| LAr ADC | 15 bits: full scale 32767 counts, pedestal 8192 | Nuovo Cim. 48 C (2025) 122, p. 3; pedestal: CREST (above) |
| Pile-up | Poisson, mean 200 interactions per colliding bunch crossing | ATLAS-TDR-028, pp. 1 and 17 (instantaneous luminosity up to 7.5 × 10³⁴ cm⁻²s⁻¹) |

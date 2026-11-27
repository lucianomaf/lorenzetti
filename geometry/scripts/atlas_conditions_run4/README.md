# Rebuilding the ATLAS Run 4 conditions without an ATLAS account

These scripts rebuild every table of `geometry/data/atlas_conditions_run4/` from sources that anyone can reach on the internet, without a
CERN or ATLAS account. They are the proof that the tables are public; see `NOTICE.md` in the data folder for what each table holds.

## What you need

- The CernVM-FS client with the repository `atlas-condb.cern.ch` mounted (free, no account; https://cernvm.cern.ch/fs/). The scripts read
  the ATLAS conditions POOL files from it directly.
- Internet access to GitHub (the public HL-LHC filling scheme), to the ATLAS Frontier server (`atlasfrontier-ai.cern.ch:8000`, the COOL conditions database, read-only) and to CREST
  (`crest.cern.ch`, the new ATLAS conditions service). Neither asks for authentication.
- Python 3 with numpy and uproot (the Lorenzetti container has numpy; `pip install --target <dir> uproot awkward` adds uproot).
- Free-running hits files (`simu_trf.py --free-running-hits`) of the minimum-bias productions of `freerunning/production/`, the blocks
  listed in `inputs/` (one block per line, `block events seed`, the blocks file of that recipe, which makes them again event by event):
  - `inputs/a3_inclusivo_v1_cells_blocks.txt`: 41 blocks of the inclusive production `a3_inclusivo_v1` (2050 events), for the cells and
    their |eta| (`lorenzetti_cells.py`);
  - `inputs/a3_hard_v1_phase_blocks.txt`: 200 blocks of the hard production `a3_hard_v1` (10 000 events with a jet or a photon of high
    transverse momentum), for the sampling phase (`phase_by_ring.py`).

  For instance, with the environment of `freerunning/production/README.md`:
  `bash ../../../freerunning/production/run_local.sh inclusive a3_inclusivo_v1 inputs/a3_inclusivo_v1_cells_blocks.txt prod` and
  `bash ../../../freerunning/production/run_local.sh hard a3_hard_v1 inputs/a3_hard_v1_phase_blocks.txt prod`.

## Steps

```
export LZT_COND_WORK=work          # results go here
export CVMFS_ROOT=/cvmfs            # where CernVM-FS is mounted
python3 sources.py                  # resolves every COOL tag and POOL file through Frontier -> work/sources.json
python3 lorenzetti_cells.py work/lorenzetti_lar_cells.csv prod/a3_inclusivo_v1/frhits/*.FRHITS.root
python3 tile_sfr.py                 # -> work/tables/tile_run4_sampling_fraction_by_ring.csv
python3 tile_noise.py               # -> work/tables/tile_run4_noise_by_ring.csv
python3 tile_scale.py               # -> work/tables/tile_run4_scale_by_ring.csv
python3 bunch_train.py ../../data/atlas_conditions_run4/hllhc_bunch_train_ip1.csv   # -> work/tables/hllhc_bunch_train_ip1.csv, and the check
python3 make_lar_tables.py         # -> work/tables/lar_run4_{constants,noise_autocorr,pulses}_by_ring.csv
python3 phase_by_ring.py '' prod/a3_hard_v1/frhits/*.FRHITS.root   # -> work/tables/{lar,tile}_phase_by_ring.csv
```

`make_lar_tables.py` reads about 1.5 GB of RTM waves from CernVM-FS and runs in some tens of minutes.

## Files

| File | What it does |
|---|---|
| `frontier.py` | read-only SQL queries to COOL through the public Frontier |
| `sources.py` | the tag and the POOL file of every folder: simulation folders from the global tag `OFLCOND-MC21-SDR-RUN4-02` (run 350000), the RTM waves and the LArPileupShape of the data global tag `CONDBR2-BLKPA-2026-01` (run 520451), the cabling tags, CREST and the Tile folder |
| `lar.py` | readers of the LAr POOL formats without Athena, the cabling, the identifier fields, the symmetrisation of the simulation constants and the completion of the FCal waves |
| `lorenzetti_cells.py` | the LAr cells and their nominal eta, from the Lorenzetti free-running hits |
| `tile.py` | the Tile blob format, default drawers and channel to cell map of the open Athena code |
| `tile_sfr.py` | the Tile sampling fraction per ring |
| `tile_noise.py` | the Tile noise per ring and gain, converted to the FENICS counts |
| `tile_scale.py` | the Tile MeV per count and photoelectrons per GeV of the simulation, per ring |
| `bunch_train.py` | the HL-LHC bunch pattern at IP1, from the public filling scheme |
| `phase_by_ring.py` | the sampling phase per ring, from the free-running hits of Lorenzetti |
| `make_lar_tables.py` | the three LAr tables per ring |

## Optional check inside ATLAS (21 November 2026)

The tables do not depend on this check, and rebuilding them does not need it. Someone with an ATLAS account can compare them with the
values Athena 25.0.21 extracts from the same conditions; we did, with these results:

- cabling of the simulation: equal for the 182 468 LAr cells;
- sampling fraction, ramp and ADC to MeV factor of the simulation: equal for the 182 468 cells (relative difference below 1e-6);
- CREST ramp and noise linked to the cells: equal;
- noise autocorrelation: equal for the 182 468 cells and the three gains (to the 3 decimals of the reference);
- RTM waves: the 529 024 EM and HEC waves are identical, and the first 156 points of the 3524 FCal waves are identical; the FCal tail built
  here differs from the one of the reference by 0.1% of the peak (median) and 2.2% at most;
- Tile sampling fraction: the medians of the 41 rings are equal.

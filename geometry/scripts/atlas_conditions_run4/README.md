# Rebuilding the ATLAS Run 4 conditions without an ATLAS account

These scripts rebuild every table of `geometry/data/atlas_conditions_run4/` from sources that anyone can reach on the internet, without a
CERN or ATLAS account. They are the proof that the tables are public; see `NOTICE.md` in the data folder for what each table holds.

## What you need

- The CernVM-FS client with the repository `atlas-condb.cern.ch` mounted (free, no account; https://cernvm.cern.ch/fs/). The scripts read
  the ATLAS conditions POOL files from it directly.
- Internet access to GitHub (the public HL-LHC filling scheme), to the ATLAS Frontier server (`atlasfrontier-ai.cern.ch:8000`, the COOL conditions database, read-only) and to CREST
  (`crest.cern.ch`, the new ATLAS conditions service). Neither asks for authentication.
- Python 3 with numpy and uproot (the Lorenzetti container has numpy; `pip install --target <dir> uproot awkward` adds uproot).
- A few free-running hits files of a Lorenzetti minimum-bias production (`simu_trf.py --free-running-hits`), for the |eta| of the cells,
  and of a production of hard events (jet or photon of high transverse momentum), for the sampling phase.

## Steps

```
export LZT_COND_WORK=work          # results go here
export CVMFS_ROOT=/cvmfs            # where CernVM-FS is mounted
python3 sources.py                  # resolves every COOL tag and POOL file through Frontier -> work/sources.json
python3 lorenzetti_cells.py work/lorenzetti_lar_cells.csv <FRHITS files>
python3 tile_sfr.py                 # -> work/tables/tile_run4_sampling_fraction_by_ring.csv
python3 tile_noise.py               # -> work/tables/tile_run4_noise_by_ring.csv
python3 tile_scale.py               # -> work/tables/tile_run4_scale_by_ring.csv
python3 bunch_train.py ../../data/atlas_conditions_run4/hllhc_bunch_train_ip1.csv   # -> work/tables/hllhc_bunch_train_ip1.csv, and the check
python3 make_lar_tables.py         # -> work/tables/lar_run4_{constants,noise_autocorr,pulses}_by_ring.csv
python3 phase_by_ring.py '' <FRHITS files of the hard production>   # -> work/tables/{lar,tile}_phase_by_ring.csv
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

## Checks done on 21 November 2026

Against the values extracted with Athena 25.0.21 from the same conditions (inside ATLAS, for validation only):

- cabling of the simulation: equal for the 182 468 LAr cells;
- sampling fraction, ramp and ADC to MeV factor of the simulation: equal for the 182 468 cells (relative difference below 1e-6);
- CREST ramp and noise linked to the cells: equal;
- noise autocorrelation: equal for the 182 468 cells and the three gains (to the 3 decimals of the reference);
- RTM waves: the 529 024 EM and HEC waves are identical, and the first 156 points of the 3524 FCal waves are identical; the FCal tail built
  here differs from the one of the reference by 0.1% of the peak (median) and 2.2% at most;
- Tile sampling fraction: the medians of the 41 rings are equal.

# Productions of the free-running library

The minimum-bias events of the free-running library (`../library/`), as they were produced: Pythia 8 with the ATLAS A3 tune at 14 TeV,
one interaction per event, simulated with the ATLAS Run 4 calorimeters of the v2 detector and written as free-running hits (the layout of
the ATLAS HITS). Everything here runs outside CERN: the image is public, the inputs come from this repository and from CVMFS.

| File | What it does |
|---|---|
| `block.sh` | one block in one of three modes (below), inside the Lorenzetti image with this repository built; writes in the current folder |
| `run_local.sh` | the blocks of a production one after the other on one machine, in the layout the library expects |
| `condor/block_job.sh`, `condor/production.sub` | the same in HTCondor jobs, with the output on EOS (example for lxplus) |
| `minbias_inelastic_A3_14TeV.cmnd` | the Pythia card (source of each setting and the differences from ATLAS in its header) |

## Modes

| Mode | What | Output (`<name>` = `<production>_<block>`) |
|---|---|---|
| `inclusive` | generates the block (seed of the block) and simulates it | `<name>.EVT.root`, `<name>.FRHITS.root`, `<name>.ok` |
| `mark` | generates the block again with the same seed, marks the hard events and checks, particle by particle, that the events are those of the stored EVT | `<name>.hard.txt` (EventNumber and bits: 1 jet, 2 photon, 4 B hadron and lepton; 0 not hard), `<name>.TRUTH.root` (the full Pythia record) |
| `hard` | generates hard events only and simulates them | `<name>.EVT.root`, `<name>.TRUTH.root`, `<name>.FRHITS.root`, `<name>.ok` (with the efficiency of the filter) |

The hard selection is that of the ATLAS sample 800831 (`generator/filters/src/HardFilter.cxx`): an anti-kt R = 0.6 jet with pT > 35 GeV, or
a photon with pT > 8 GeV that does not come from a hadron, or a weakly decaying B hadron with pT > 5 GeV and an electron or muon with pT > 4 GeV (the |eta| ranges in the code). The library takes
the events with bit 1 or 2 out of the inclusive sample (as the veto of the ATLAS sample 900311) and mixes the hard sample in with the fraction
of hard events of the marking (`../library/README.md`).

Generation: `gen_minbias.py`, one interaction per event, all particles (`--eta-max 10`, no pT cut), vertex spread 50 mm and 175 ps. The
Geant4 seed is not an option of `simu_trf.py`: every block uses the fixed seed of the `ComponentAccumulator` (512); the blocks differ by
their events. Simulation: `simu_trf.py --detector v2` with the options listed in `block.sh`, 4 threads, a watchdog of 110 minutes.

## The productions of the library

Blocks of 50 events; the Pythia seed of a block is the base of its production plus the block number.

| Production | Mode | Blocks | Seed | Events |
|---|---|---|---|---|
| `a3_inclusivo_v1` | `inclusive`, then `mark` | `0000` to `1999` | 100000 + block | 100 000 |
| `a3_inclusivo_v2` | `inclusive`, then `mark` | `00000` to `19999` | 102000 + block | 1 000 000 |
| `a3_hard_v1` | `hard` | `0000` to `9999` | 200000 + block | 500 000 hard |

A blocks file has one block per line, `block events seed`; for instance the first production:

```bash
for k in $(seq 0 1999); do printf '%04d 50 %d\n' $k $((100000 + k)); done > v1_blocks.txt
```

## Image and build

The image is `lorenzetti/lorenzetti@sha256:800d38bf762396054e970117f13c58ed937a6f483d82ae89ce69dfd2e2707466` (created on 3 October
2026), the one the productions used; give it by its digest, not by `latest`, which moves. Build this repository in it with GeoModel
(`../../README.md`, "Build with GeoModel") and get the three external files of the v2 detector listed there (the Run 4 geometry dump, the
field map and the FCal shower libraries; the EMEC field map too).

## On one machine

```bash
docker run --rm -it -v $HOME:$HOME -w $HOME lorenzetti/lorenzetti@sha256:800d38bf762396054e970117f13c58ed937a6f483d82ae89ce69dfd2e2707466 bash
export LZT_SETUP=$HOME/lorenzetti/build/lzt_setup.sh
export LZT_RUN4_GEOMETRY=$HOME/atlas_geometry/geometry-ATLAS-P2-RUN4-01-01-00-sim.db
export LZT_FIELD_MAP=$HOME/atlas_field/bfieldmap_7730_20400_14m.root
bash lorenzetti/freerunning/production/run_local.sh inclusive a3_inclusivo_v1 v1_blocks.txt $HOME/minbias
bash lorenzetti/freerunning/production/run_local.sh mark a3_inclusivo_v1 v1_blocks.txt $HOME/minbias
```

A block of 50 inclusive events takes about an hour on 4 cores (70 to 86 s per event); a hard event about twice as long.

## In HTCondor (lxplus)

The jobs bring the build as a package: a `.tgz` of the home folder where it was made, with the repository built, GeoModel and the
geometry dump (`atlas_geometry/`), mounted at the same path in the job, because the build keeps absolute paths. `condor/make_package.sh`
makes it on the machine of the build, after checking the files the jobs need, and says how to copy it to EOS and check it there. The image goes as a `.sif`
(`apptainer build lorenzetti.sif docker://lorenzetti/lorenzetti@sha256:800d38bf…`, in a job: the interactive nodes of lxplus cannot build
it). Fill in the environment line of `condor/production.sub` and submit, from `condor/`:

```bash
condor_submit production.sub mode=inclusive production=a3_inclusivo_v1 blocks=v1_blocks.txt
```

The jobs write `<PROD_OUT>/<production>/{evt,frhits,truth,hard,ok}/`; a block already on EOS is skipped, so the same command resumes.

## Validation

The recipe was run on one machine (`run_local.sh`, the image above, a single build of this branch) with the seeds of the productions and
compared with them, event by event (22 November 2026):

| Block | Run | Result |
|---|---|---|
| `a3_inclusivo_v1_0000`, mode `inclusive` | the first 10 events | the FRHITS identical to the production in all 40 columns of `CollectionTree` |
| the same block, mode `mark` | the 10 events, checked against the EVT of the run above | the marking identical to the first 10 lines of the production |
| `a3_hard_v1_0000`, mode `hard` | the first 10 hard events (2131 generated) | the FRHITS identical to the production in all 40 columns |

The productions were simulated with the build of `cae76f1` and the hard events generated with that of `1686fd2`; the later commits change
only the event generator and the EVT format (`generator/`), so one build of this branch gives the same blocks. A run of 10 events gives the
first 10 events of the block of 50, as the comparison shows, so the short runs stand for the whole blocks.

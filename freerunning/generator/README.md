# Free-running generator

Digitised samples at 40 MHz of selected calorimeter channels with the pile-up of the HL-LHC (Run 4). The generator reads the free-running
library (`../library/`) and the public conditions in `geometry/data/atlas_conditions_run4/`.

```bash
python3 generate.py --library <library> --cell "LBA14 A8" --output seq.root          # one Tile cell (both PMTs), 9 orbits, mu 200
python3 generate.py --library <library> --cell "EMB2 0.5" --orbits 2 --output seq.root  # a LAr cell of the ring of |eta| 0.5, side and phi drawn
python3 generate.py --library <library> --cluster 0.5 --layers em --output seq.root    # a 0.075 x 0.175 window of the EM layers
python3 generate.py --library <library> --tower 5 --output seq.root                   # a projective tower of 0.1 x 2 pi / 64
python3 generate.py --library <library> --list-cells EMB2                              # the names and identifiers of a layer
```

## Smallest end-to-end run

From nothing to a sequence: one inclusive block and one hard block of 10 events, a library made of them and one sequence. It
needs this repository built with GeoModel in the Lorenzetti image and the external files of the v2 detector (`../production/README.md`,
"Image and build"). A library this small only allows a low mu: each event can be used at most 128 times, once per copy.

```bash
docker run --rm -it -e HOME=$HOME -v $HOME:$HOME -w $HOME \
  lorenzetti/lorenzetti@sha256:800d38bf762396054e970117f13c58ed937a6f483d82ae89ce69dfd2e2707466 bash
export LZT_SETUP=$HOME/lorenzetti/build/lzt_setup.sh
export LZT_RUN4_GEOMETRY=$HOME/atlas_geometry/geometry-ATLAS-P2-RUN4-01-01-00-sim.db
export LZT_FIELD_MAP=$HOME/atlas_field/bfieldmap_7730_20400_14m.root
R=$HOME/lorenzetti/freerunning
mkdir -p $HOME/freerunning_example && cd $HOME/freerunning_example
echo "0000 10 100000" > inclusive_blocks.txt    # block, events, Pythia seed
echo "0000 10 200000" > hard_blocks.txt
bash $R/production/run_local.sh inclusive example_inclusive inclusive_blocks.txt prod
bash $R/production/run_local.sh mark example_inclusive inclusive_blocks.txt prod
bash $R/production/run_local.sh hard example_hard hard_blocks.txt prod
python3 $R/library/make_list.py list.csv inclusive:prod/example_inclusive:100000 hard:prod/example_hard:200000
bash $R/library/run_local.sh list.csv library work
source $LZT_SETUP
python3 $R/generator/generate.py --library library --cell "LBA14 A8" --mu 0.05 --orbits 1 --output seq.root
```

On a desktop computer (WSL, 4 threads per block) it took about 1 h 30 min: 55 min for the inclusive block, 2.5 min for the marking,
29 min for the hard block, 2 min for the library and a few seconds for the sequence. The tests (`tests/run_tests.py`) need a larger
library, with about 100 inclusive and 100 hard events: with 10 events a sample runs out of its 128 copies.

## Options

| Option | Meaning |
|---|---|
| `--cell` | one cell: short form `"Tile A13"`, `"EMB2 0.5"` (side and phi drawn with the seed); full form `"LBA14 A8"`, `"EMB2 A eta=0.5125 phi=0.0245"`, or the identifier |
| `--cluster ETA`, `--cluster-deta`, `--cluster-dphi` | the cells of `--layers` whose centre falls in a window around \|eta\| ETA (side and phi drawn; the width in phi in steps of 2 pi / 64) |
| `--tower N` | the cells of `--layers` in tower N of 0.1 x 2 pi / 64 (LAr: centre in the tower; Tile: tower field N), side and module drawn |
| `--layers` | EMB0-3, EMEC0-3, HEC0-3, FCAL1-3, TileA, TileBC, TileD, TileE and the groups em (default of `--cluster`), had, tile, all (default of `--tower`) |
| `--mu` (200), `--orbits` (9), `--seed` (1), `--first-slot` (0) | pile-up, kept orbits, seed, slot of the filling scheme for the first sample |
| `--truth-unit adc\|visible_mev\|total_mev` | unit of the truth `a` (default: counts of the high gain) |
| `--gains hg,lg` | gains written (default both) |
| `--reuse on\|off` | when a sample of the library runs out: rotated and reflected copies (default), or stop |
| `--tile-tail on\|off`, `--tile-tail-tau` | long tail of the Tile pulse (AC coupling, zero area), off by default; time constant 50 us |
| `--hard-fraction` | hard events per interaction (default: the fraction measured by the marking, which the library writes in `samples.json`) |
| `--no-hits` | do not write the hits tree |
| `--library`, `--support`, `--output` | paths |

## Output (ROOT)

| Tree | One entry per | Content |
|---|---|---|
| `parameters` (+ `parameters_json`) | key | version, commit, command, seed, mu, orbits, filling scheme, library and samples, hard fraction, passes over the library, options |
| `channels` | channel | index, identifier, name, subdetector, layer, ring, eta, phi, phase, pedestal, counts per MeV of visible energy and noise sigma of each gain |
| `samples` | 25 ns sample | `n`, `bcid`, `orbit`, `filled`, `n_interactions`; per channel `r_hg`, `r_lg` (integer counts), `saturated_hg`, `saturated_lg`, `a` (truth) |
| `hits` | hit in the kept samples | channel, sample, visible energy, time, bunch crossing of origin, library event, copy |

## Terms

| Term | Meaning |
|---|---|
| marking | the hard-event bits of each inclusive event (1 jet, 2 photon, 4 B hadron with a lepton), written by the `mark` mode of `../production/` |
| pass | one use of every event of a sample; when the kept samples need more events than the sample has, a new pass starts, with a new order and another copy |
| copy | an event rotated by k x 2 pi / 64 in phi, reflected between the two sides, or both: 128 per event, counting the event as it is |
| ring | the cells that share their conditions: LAr, the same layer and \|eta\|; Tile, the same section, sampling and tower |
| DSID | the ATLAS dataset identifier of a simulated sample; 800831 is the high-pT minimum-bias sample whose selection `HardFilter` reproduces |

## Method

- **Bunch crossings.** HL-LHC standard filling scheme (`hllhc_bunch_train_ip1.csv`), with Poisson interactions of mean `--mu` in the
  colliding bunch crossings. Each interaction is hard with probability `--hard-fraction` and inclusive otherwise. By default this is the
  fraction of hard events in the marking of the inclusive production, read from the library.
- **Events.** In the kept orbits each event of a sample is used once per pass, in an order shuffled from the seed, with a new shuffle in
  each pass. Pass p uses copy `copy_order[p]`: the event as it is first, then the 127 copies (rotation by k x 2 pi / 64, k = 0..63, and
  reflection between the sides) in an order from the seed, the same copy in the whole detector. Rotation is exact in Tile, EMB, EMEC and
  the HEC up to |eta| 2.5; the HEC rings of 32 cells rotate by k / 2 cells; the FCal does not rotate. Reflection goes by the identifier, in
  every subdetector. The warm-up (the longest pulse; 10 tau with the Tile tail) and the bunch crossings after the kept ones draw events with
  replacement, from a stream of their own. Each takes a random copy among those the kept samples use, so no copy is read for the warm-up
  alone.
- **Signal.** Each hit adds counts x pulse(25 (n - m) + phase - t) to sample n, with the pulse peak 1 at t = 0. LAr: the pulse of the ring,
  linearly interpolated on 0.5 ns. Tile: CR-RC^4, tau = 10 ns. The phase is that of the ring, continuous. Tile photoelectron statistics as
  in Athena (TileHitVecToCntTool, type 2): Gaussian above 20 photoelectrons, Poisson with the single-photoelectron spectrum below. It acts
  per channel and bunch crossing, on the energy within +-62.5 ns, on the amplitude only; the draws for all the bunch crossings of a channel
  are made together.
- **Scale.** LAr: counts = E_visible / fsampl / (adc2mev_hg_mc x ramp_r4 / ramp_hg_mc). Tile: counts = E_visible x sampling fraction /
  (MeV per count of the Run 4 simulation x 0.525), with 40 times fewer counts in the low gain.
- **Noise and digitisation.** LAr: white Gaussian noise through an AR(4) filter that reproduces the correlation of the ring at lags 1-4,
  with the sigma of the Phase-II front end. Tile: one Gaussian, with the sigma that keeps the variance of the double Gaussian of the table.
  Rounding floor(x + 0.5), saturation below 0 and above 4095 (Tile) or 32767 (LAr); pedestals 90 and 8192.
- **Truth.** Each hit goes to the sample where the peak of its pulse falls, without photoelectron statistics.
- **Reproducibility.** The interactions and events come from one stream for the whole detector. The noise and the photoelectron
  statistics of a channel come from the seed and its identifier, so a channel has the same samples whatever the selection.
- **Known limits.** The library has 20 special Tile PMTs (C10 with one PMT in modules 39-42 and 55-58, D4 absent in EBA15 and EBC18). A
  copy that rotates a normal cell onto them loses the PMT that does not exist there, and the D5 of EBA15 and EBC18 receives only its own
  part.

## Tests

`python3 tests/run_tests.py <small library> <work folder>`, in the Lorenzetti image (no pytest). They check: the seed; a channel alone
against the same channel in a tower; one use per pass and a new order in the next one; `--reuse off`; the copies of every ring; the truth
against the hits; the noiseless samples against the pulse of the hits; the noise (pedestal, sigma, LAr correlation); the zero area of the
Tile tail; saturation.

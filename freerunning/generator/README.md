# Free-running generator

Digitised samples at 40 MHz of selected calorimeter channels with the pile-up of the HL-LHC (Run 4), from the free-running library
(`../library/`) and the public conditions of `geometry/data/atlas_conditions_run4/`. The options follow those of the ATLAS generator of
free-running sequences.

```bash
python3 generate.py --library <library> --cell "LBA14 A8" --output seq.root          # one Tile cell (both PMTs), 9 orbits, mu 200
python3 generate.py --library <library> --cell "EMB2 0.5" --orbits 2 --output seq.root  # a LAr cell of the ring of |eta| 0.5, side and phi drawn
python3 generate.py --library <library> --cluster 0.5 --layers em --output seq.root    # a 0.075 x 0.175 window of the EM layers
python3 generate.py --library <library> --tower 5 --output seq.root                   # a projective tower of 0.1 x 2 pi / 64
python3 generate.py --library <library> --list-cells EMB2                              # the names and identifiers of a layer
```

## Options

| Option | Meaning |
|---|---|
| `--cell` | one cell: short form `"Tile A13"`, `"EMB2 0.5"` (side and phi drawn with the seed); full form `"LBA14 A8"`, `"EMB2 A eta=0.5125 phi=0.0245"`, or the identifier |
| `--cluster ETA`, `--cluster-deta`, `--cluster-dphi` | the cells of `--layers` whose centre falls in a window around \|eta\| ETA (side and phi drawn; the width in phi in steps of 2 pi / 64) |
| `--tower N` | the cells of `--layers` in the tower N of 0.1 x 2 pi / 64 (LAr: centre in the tower; Tile: tower field N), side and module drawn |
| `--layers` | EMB0-3, EMEC0-3, HEC0-3, FCAL1-3, TileA, TileBC, TileD, TileE and the groups em (default of `--cluster`), had, tile, all (default of `--tower`) |
| `--mu` (200), `--orbits` (9), `--seed` (1), `--first-slot` (0) | pile-up, kept orbits, seed, slot of the filling scheme of the first sample |
| `--truth-unit adc\|visible_mev\|total_mev` | unit of the truth `a` (default: counts of the high gain) |
| `--gains hg,lg` | gains written (default both) |
| `--reuse on\|off` | when a sample of the library runs out: rotated and reflected copies (default), or stop |
| `--tile-tail on\|off`, `--tile-tail-tau` | long tail of the Tile pulse (AC coupling, zero area), off by default; time constant 50 us |
| `--no-hits` | do not write the hits tree |
| `--library`, `--support`, `--output` | paths |

## Output (ROOT)

| Tree | One entry per | Content |
|---|---|---|
| `parameters` (+ `parameters_json`) | key | version, commit, command, seed, mu, orbits, filling scheme, library and samples, hard fraction, passes over the library, options |
| `channels` | channel | index, identifier, name, subdetector, layer, ring, eta, phi, phase, pedestal, counts per MeV of visible energy and noise sigma of each gain |
| `samples` | 25 ns sample | `n`, `bcid`, `orbit`, `filled`, `n_interactions`; per channel `r_hg`, `r_lg` (integer counts), `saturated_hg`, `saturated_lg`, `a` (truth) |
| `hits` | hit of the kept samples | channel, sample, visible energy, time, bunch crossing of origin, event of the library, copy |

## Method

- **Bunch crossings.** HL-LHC standard filling scheme (`hllhc_bunch_train_ip1.csv`), Poisson interactions of mean `--mu` in the colliding
  bunch crossings; each interaction is hard with probability 0.338% (the fraction of the hard filter in the inclusive production),
  inclusive otherwise.
- **Events.** In the kept orbits each event of a sample is used once per pass, in an order shuffled from the seed, a new shuffle in each
  pass; pass p uses copy `copy_order[p]`: the event as it is first, then the 127 copies (rotation by k x 2 pi / 64, k = 0..63, and
  reflection between the sides) in an order from the seed, the same copy in the whole detector. Rotation: exact in Tile, EMB, EMEC and the
  HEC up to |eta| 2.5; k / 2 cells in the HEC rings of 32 cells; none in the FCal. Reflection: by the identifier, every subdetector. The
  warm-up (the longest pulse; 10 tau with the Tile tail) and the bunch crossings after the kept ones draw events with replacement, with a
  random copy, from a stream of their own.
- **Signal.** Each hit adds counts x pulse(25 (n - m) + phase - t) to sample n; pulse peak 1 at t = 0; LAr: the pulse of the ring, linear
  interpolation on 0.5 ns; Tile: CR-RC^4, tau = 10 ns; phase of the ring, continuous. Tile photoelectron statistics as Athena
  (TileHitVecToCntTool, type 2: Gaussian above 20 photoelectrons, Poisson with the single-photoelectron spectrum below), per channel and
  bunch crossing, energy within +-62.5 ns, on the amplitude only.
- **Scale.** LAr: counts = E_visible / fsampl / (adc2mev_hg_mc x ramp_r4 / ramp_hg_mc). Tile: counts = E_visible x sampling fraction / (MeV
  per count of the Run 4 simulation x 0.525), low gain 40 times fewer.
- **Noise and digitisation.** LAr: white Gaussian through an AR(4) filter that reproduces the correlation of the ring at lags 1-4, sigma of
  the Phase-II front end. Tile: one Gaussian, the sigma that keeps the variance of the double Gaussian of the table. Rounding floor(x + 0.5),
  saturation below 0 and above 4095 (Tile) or 32767 (LAr); pedestals 90 and 8192.
- **Truth.** Each hit goes to the sample where the peak of its pulse falls; no photoelectron statistics.
- **Reproducibility.** The interactions and events come from one stream for the whole detector; the noise and the photoelectron statistics
  of a channel from the seed and its identifier: a channel has the same samples whatever the selection.
- **Known limits.** The 20 special Tile PMTs of the library (C10 with one PMT in modules 39-42 and 55-58, D4 absent in EBA15 and EBC18):
  a copy that rotates a normal cell onto them loses the PMT that does not exist there, and the D5 of EBA15 and EBC18 receives only its own
  part.

## Tests

`python3 tests/run_tests.py <small library> <work folder>` (Lorenzetti image; no pytest): seed, a channel alone against the same channel in
a tower, one use per pass and a new order in the next one, `--reuse off`, the copies of every ring, the truth against the hits, the
noiseless samples against the pulse of the hits, the noise (pedestal, sigma, LAr correlation), the zero area of the Tile tail, saturation.

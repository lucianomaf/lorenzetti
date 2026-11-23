# Free-running library

The calorimeter hits of a minimum-bias production of Lorenzetti, reorganised by channel for a free-running sequence generator. The
generator reads one channel at a time: every hit of that channel, in every event of the library, without reading the other channels.

## Input

FRHITS files (`simu_trf.py --free-running-hits --free-running-atlas-format`), one per block of 50 events, of two samples:

- **inclusive**: Pythia 8 minimum bias (tune A3), every event simulated. Each block has a hard-event marking file (`<event> <bits>` per
  line; bit 1 jet, bit 2 photon, bit 4 B hadron with a lepton). Events with a jet or a photon leave the inclusive sample, since the hard
  sample stands for them; events with only the B-with-lepton bit stay.
- **hard**: the same Pythia 8 minimum bias, with only the events that pass the hard filter sent to Geant4. Every event enters.

The input list has one block per line, `sample,production,block,seed,FRHITS,marking` (marking `-` in the hard sample), with the
inclusive sample first. The order of the lines is the order of the global event index. `make_list.py` writes it from the folders of the productions
(`../production/README.md`), which hold `frhits/` and, for an inclusive production, `hard/` with the markings:

```bash
python3 make_list.py list.csv inclusive:<folder>/a3_inclusivo_v1:100000 inclusive:<folder>/a3_inclusivo_v2:102000 hard:<folder>/a3_hard_v1:200000
```

## Layout

| Path | Content |
|---|---|
| `per_channel/<subdetector>_<NNN>.root` | pieces of `tile`, `emb`, `emec`, `hec`, `fcal`, each a contiguous range of whole channels: tree `channels` (`channel`, `start`, `n`) and tree `hits` (`event`, `energy`, `time`), sorted by channel, event and time |
| `support/per_channel_index.json` | per subdetector, for each piece: file, first and last channel, number of channels and of hits |
| `support/samples.json` | name, first global index and number of events of each sample; in the hard sample also `hard_fraction`, the events with any bit of the marking over the marked inclusive events, which the generator uses as the fraction of hard interactions |
| `support/productions.json`, `support/manifest.root` | per global index: production, block, event number in the block, Pythia seed of the block, hard-event bits (-1 in the hard sample) |
| `support/tile_cells.csv` | per PMT: identifier, section, side, module, tower, sampling, PMT, eta, phi |
| `support/lar_cells.csv` | per cell: identifier, subdetector, layer, ring, eta, phi; `layer` and `ring` are the keys of the LAr tables of `geometry/data/atlas_conditions_run4/` |

- **channel**: the identifier written by Lorenzetti. Tile: the 64-bit `cellID` of the hit, which is the PMT. LAr: the upper 32 bits of
  the `cellID`, the ATLAS cell identifier.
- **event**: global index in the library (int32).
- **energy**: visible MeV, copied from the FRHITS (float32).
- **time**: ns, copied from the FRHITS (float32). The time of flight from the interaction point is already removed.
- **eta, phi**: the nominal centre of the cell. The D5 PMTs of the two special extended-barrel modules (EBA15, EBC18) also read the D4
  region, as in ATLAS, so their hits carry two cell centres; the table keeps the D5 one.
- **ring**: Tile, section, sampling and tower. LAr, from `geometry/scripts/atlas_conditions_run4/lar.py`, with the FCal in bands of |eta|
  0.1 wide of the cell.

## Steps

| Step | Script | Output |
|---|---|---|
| events | `make_support.py events <input list> <library> <blocks.csv>` | `samples.json`, `productions.json`, `manifest.root` and the list of block jobs |
| block (one job per block) | `block_to_channels.py <FRHITS> <marking or -> <first index> <events> <out.root>` | the hits of the block by channel, one tree per subdetector, and the channel trees |
| plan | `merge_pieces.py plan <block files> <library> <pieces.csv> [hits per piece]` | `per_channel_index.json` and the list of piece jobs |
| piece (one job per piece) | `merge_pieces.py piece <block files> <library> <subdetector> <n>` | `per_channel/<subdetector>_<NNN>.root` |
| cells | `make_support.py cells <block files> <library>` | `tile_cells.csv`, `lar_cells.csv`; checks that every ring is in the conditions tables |
| check | `check_library.py <input list> <library>` | per channel, the number of hits and the energy sum against the FRHITS, read again without the code of the other steps |

All steps use only ROOT (PyROOT, `RDataFrame`) and numpy, as in the Lorenzetti image. To run them all in sequence on one machine:

```bash
docker run --rm -v $PWD:/lzt:ro -v <data>:/w lorenzetti/lorenzetti:latest \
  bash /lzt/freerunning/library/run_local.sh /w/list.csv /w/library /w/work
```

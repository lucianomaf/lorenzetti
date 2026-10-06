#!/usr/bin/env python3

import multiprocessing
import argparse
import sys

from pathlib                import Path
from GaugiKernel.constants  import MINUTES
from GaugiKernel            import LoggingLevel, get_argparser_formatter
from G4Kernel               import ComponentAccumulator, EventReader
from RootStreamBuilder      import recordable
from CaloHitBuilder         import CaloHitBuilder
from RootStreamBuilder      import RootStreamHITMaker

from geometry import DetectorConstruction_v1
from reco import update_args_from_file,merge_args_from_file


"""
Script: simu_trf.py
Purpose: Runs the full Geant4 simulation step in the Lorenzetti framework.
         Reads generated event files (EVT), simulates particle interactions with the detector,
         and produces Hit files (HIT).
Usage:
    simu_trf.py -i input.EVT.root -o output.HIT.root -nt <threads>
"""

def parse_args():
    """
    Parses command-line arguments for the simulation job.

    Returns:
        argparse.Namespace: Arguments specifying inputs/outputs, thread count,
                            magnetic field options, and execution hooks.
    """
    parser = argparse.ArgumentParser(
        description='',
        add_help=False,
        formatter_class=get_argparser_formatter())

    parser.add_argument('-i', '--input-file', action='store',
                            dest='input_file', required=True,
                            help="The input file or folder to run the job")
    parser.add_argument('-o', '--output-file', action='store',
                        dest='output_file', required=True,
                        help="The output file.")
    parser.add_argument('--nov', '--number-of-events', action='store',
                        dest='number_of_events', required=False,
                        type=int, default=-1,
                        help="The total number of events to run.")
    parser.add_argument('-nt', '--number-of-threads', action='store',
                        dest='number_of_threads', required=False,
                        type=int, default=multiprocessing.cpu_count(),
                        help="The number of threads")
    parser.add_argument('--enable-magnetic-field', action='store_true',
                        dest='enable_magnetic_field', required=False,
                        help="Enable the magnetic field.")
    parser.add_argument('--active-energy-only', action='store_true',
                        dest='active_energy_only', required=False,
                        help="Keep only the energy deposited in the active medium of the calorimeters "
                             "(liquid argon, plastic scintillator), as in the ATLAS hits. By default the "
                             "whole energy deposited in the cell volume, absorber included, is kept.")
    parser.add_argument('--birks-law', action='store_true',
                        dest='birks_law', required=False,
                        help="Apply Birks' law to the energy deposited in the plastic scintillator and in the "
                             "liquid argon, with the constants of the ATLAS simulation. Use it together with "
                             "--active-energy-only to follow the ATLAS hits.")
    parser.add_argument('--enable-solenoid-field', action='store_true',
                        dest='enable_solenoid_field', required=False,
                        help="Enable a uniform 2 T axial field only inside the ATLAS central solenoid "
                             "(R < 1.23 m, |z| < 2.9 m), with no field in the calorimeters. "
                             "Cannot be combined with --enable-magnetic-field. Generate the events with "
                             "--pt-min-charged 0 --pt-min-neutral 0 to follow ATLAS.")
    parser.add_argument('--tile-atlas-geometry', action='store_true',
                        dest='tile_atlas_geometry', required=False,
                        help="Build the tile calorimeter as in ATLAS: scintillating tiles normal to the beam line, "
                             "stacked along z in periods of 18 mm (14 mm steel, 3 mm tile, 1 mm clearance), "
                             "with the ATLAS layer radii and z extent (long barrel |z| < 2808 mm, extended barrel "
                             "3554 < |z| < 6110 mm, the ITC aluminium block in between). Without --tile-atlas-cells "
                             "the cells keep the default eta x phi segmentation. Off by default. Use the same option "
                             "in digit_trf.py.")
    parser.add_argument('--tile-atlas-cells', action='store_true',
                        dest='tile_atlas_cells', required=False,
                        help="Use the ATLAS tile cells (A1-A16, BC1-BC8, B9, D0-D6; JINST 3 (2008) S08003, fig. 5.12) "
                             "instead of the eta x phi grid, so that every point of the tile volume belongs to a "
                             "cell. Needs --tile-atlas-geometry. Hits only: the digitization does not support these "
                             "cells yet. Off by default.")
    parser.add_argument('--tile-dual-readout', action='store_true',
                        dest='tile_dual_readout', required=False,
                        help="Read each tile cell with two photomultipliers, as in ATLAS: the energy of each step is "
                             "shared between the PMT 0 and the PMT 1 (the one on the side of larger phi) according "
                             "to its azimuthal position across the module, and both are stored with the cell "
                             "(edep_pmt0, edep_pmt1). Needs --tile-atlas-cells. Off by default.")
    parser.add_argument('--tile-pmt-split', action='store', choices=['ushape', 'linear'], default='ushape',
                        dest='tile_pmt_split', required=False,
                        help="Sharing between the two PMTs with --tile-dual-readout: 'ushape' (default), the "
                             "measured U-shape of the ATLAS simulation (Athena, TileGeoG4SDCalc::"
                             "Tile_1D_profileRescaled); 'linear', 0.5 +- 0.2 from the centre to the edges of the "
                             "module, the two PMTs adding up to the cell energy.")
    parser.add_argument('--atlas-material-in-front', action='store_true',
                        dest='atlas_material_in_front', required=False,
                        help="Material in front of the barrel electromagnetic calorimeter as in ATLAS, in radiation "
                             "lengths at normal incidence: the solenoid (0.66 X0), the cryostat wall in front of the "
                             "presampler (0.76 X0 of aluminium instead of 0.45) and the material between the presampler "
                             "and the accordion (0.60 X0); JINST 3 (2008) S08003, fig. 5.1 and sec. 2.1.1. The inner "
                             "detector is not included. Simulation only (the cells do not change). Off by default.")
    parser.add_argument('--tile-atlas-itc', action='store_true',
                        dest='tile_atlas_itc', required=False,
                        help="Build the gap between the tile barrel and the extended barrel as in ATLAS instead of the "
                             "aluminium block: the plug of the ITC (D4 and C10, steel and scintillator, passive) and the "
                             "cables and services (aluminium with 5%% of the volume); JINST 3 (2008) S08003, sec. 5.5 and "
                             "fig. 5.12. Needs --tile-atlas-geometry. Simulation only. Off by default.")
    parser.add_argument('--atlas-emb', action='store_true',
                        dest='atlas_emb', required=False,
                        help="Barrel EM calorimeter with the ATLAS absorber composition: lead of 1.53 mm below "
                             "|eta| = 0.8 and 1.13 mm above, two 0.2 mm steel sheets, glue and electrode, liquid argon "
                             "of 2 x 2.1 mm, 470 mm of depth (JINST 3 (2008) S08003, sec. 5.2 and fig. 5.1). Radial "
                             "shells, not the accordion. Use the same option in digit_trf.py. Off by default.")
    parser.add_argument('--atlas-hec', action='store_true',
                        dest='atlas_hec', required=False,
                        help="HEC as in ATLAS, with the nominal values of the ATLAS geometry (without the shift of the installed detector): "
                             "four samplings HEC0-HEC3 in two wheels (z = 4277-5093.5 and 5134-6095 mm), copper plates "
                             "of 25 and 50 mm with first plates of 12.5 and 25 mm, 8.5 mm gaps, and the cells of the "
                             "ATLAS identifier dictionary found from the radius of the readout pads of each block, as "
                             "in the ATLAS simulation (JINST 3 (2008) S08003, sec. 5.3; Athena HECGeometry). "
                             "Use the same option in digit_trf.py. Off by default.")
    parser.add_argument('--atlas-emec', action='store_true',
                        dest='atlas_emec', required=False,
                        help="EM end-cap as in ATLAS, with the nominal values of the ATLAS geometry (without the shift of the "
                             "installed detector): two wheels (1.375-2.5 and 2.5-3.2 in eta, 3 mm apart) from z = 3702 to "
                             "4216 mm, the cones approximated by steps; radial bands of 100 mm made of layers of an "
                             "equivalent absorber (lead of 1.7 or 2.2 mm, steel, prepreg and electrode) and liquid "
                             "argon with the gap of ATLAS at each radius; EMEC3 only in the outer wheel; the cells of "
                             "the ATLAS identifier dictionary found as in the ATLAS simulation (JINST 3 (2008) S08003, "
                             "sec. 5.2; Athena LArG4EC EnergyCalculator). The end-cap presampler is moved to just in "
                             "front of the EMEC (provisional). Use the same option in digit_trf.py. Off by default.")
    parser.add_argument('--atlas-endcap-cryostat', action='store_true',
                        dest='atlas_endcap_cryostat', required=False,
                        help="Outer cylinders of the end-cap cryostats as in ATLAS: warm vessel of 20 mm of aluminium "
                             "(r = 2260-2280 mm), cold vessel of 35 mm (r = 2140-2175 mm) and the liquid argon between "
                             "the EMEC and the cold vessel (dead material), from z = 3717.5 mm to the back of each vessel; "
                             "LAr calorimeter TDR, CERN/LHCC 96-41, fig. 3-3 and 5-i, and JINST 3 (2008) S08003, sec. "
                             "5.4. Simulation only (the cells do not change). Off by default.")
    parser.add_argument('--atlas-barrel-cryostat', action='store_true',
                        dest='atlas_barrel_cryostat', required=False,
                        help="Outer part of the barrel cryostat as in ATLAS instead of the two 100 mm aluminium shells: "
                             "cold and warm outer cylinders of 30 mm (r = 2140-2170 and 2220-2250 mm), end walls, the "
                             "steps for the feed-throughs and the dead liquid argon between the EM barrel and the cold "
                             "vessel; LAr calorimeter TDR, CERN/LHCC 96-41, fig. 3-2. Needs --tile-atlas-itc. "
                             "Simulation only (the cells do not change). Off by default.")
    parser.add_argument('-t', '--timeout', action='store',
                        dest='timeout', required=False, type=int, default=240,
                        help="Event timeout in minutes")
    parser.add_argument('-l', '--output-level', action='store',
                        dest='output_level', required=False,
                        type=str, default='INFO',
                        help="The output level messenger.")
    parser.add_argument('--pre-init', action='store',
                        dest='pre_init', required=False, default="''",
                        help="The preinit command")
    parser.add_argument('--pre-exec', action='store',
                        dest='pre_exec', required=False, default="''",
                        help="The preexec command")
    parser.add_argument('--post-exec', action='store',
                        dest='post_exec', required=False, default="''",
                        help="The postexec command")
    parser.add_argument('--save-all-hits', action='store_true',
                        dest='save_all_hits', required=False,
                        help="Save all hits into the output file.")
    parser.add_argument('--free-running-hits', action='store_true',
                        dest='free_running_hits', required=False,
                        help="Write the tile hits in the layout of the ATLAS HITS ntuple (tree CollectionTree: "
                             "EventNumber, TileCalHit_cellID/energy/time/eta/phi/sampling/side/module/tower, plus the "
                             "extra TileCalHit_tglobal; the LAr branches are written empty) instead of the standard HIT "
                             "stream. As in the ATLAS simulation: one entry per PMT, with the ATLAS identifier (pmt_id, "
                             "64-bit), the visible energy (active medium and Birks' law, switched on by this option), the "
                             "time of TileGeoG4SDCalc (time of flight and fibre corrections, PMT delays, hits after "
                             "350.5 ns kept at 99995 ns) and its time bins (0.5 ns within +-75.25 ns, 5 ns outside). "
                             "Needs --tile-atlas-geometry and --tile-atlas-cells; switches on --tile-dual-readout (U-shape). "
                             "Off by default.")
    parser.add_argument('--neutron-time-cut', action='store', type=float, default=None,
                        dest='neutron_time_cut', required=False,
                        help="Kill every neutron this time (ns) after the start of the event, as the ATLAS simulation does "
                             "with 150 ns (Sim.NeutronTimeCut; arXiv:1005.4568, sec. 5.5). Without it, the default of "
                             "Geant4 (10 us in FTFP_BERT). --free-running-hits sets 150 ns unless given here.")
    parser.add_argument('--free-running-no-binning', action='store_true',
                        dest='free_running_no_binning', required=False,
                        help="With --free-running-hits: one entry per PMT and Geant4 step, without the time bins.")
    parser.add_argument('--free-running-keep-hits', action='store_true',
                        dest='free_running_keep_hits', required=False,
                        help="With --free-running-hits, for checks: keep also the standard HIT stream (tree "
                             "CollectionTree) and write the free-running hits in the tree FreeRunningTree.")
    parser.add_argument('--dry-run', action='store_true',
                        dest='dry_run', required=False,
                        help="Run the script without executing the main logic.")
    return merge_args_from_file(parser)


def main(logging_level: str,
         input_file: str | Path,
         output_file: str | Path,
         pre_init: str,
         pre_exec: str,
         post_exec: str,
         enable_magnetic_field: bool,
         enable_solenoid_field: bool,
         tile_atlas_geometry: bool,
         tile_atlas_cells: bool,
         atlas_material_in_front: bool,
         tile_atlas_itc: bool,
         atlas_emb: bool,
         atlas_hec: bool,
         atlas_emec: bool,
         atlas_endcap_cryostat: bool,
         atlas_barrel_cryostat: bool,
         tile_dual_readout: int,
         active_energy_only: bool,
         birks_law: bool,
         save_all_hits : bool,
         free_running_hits : bool,
         free_running_binning : bool,
         free_running_keep_hits : bool,
         neutron_time_cut : float,
         timeout: int,
         number_of_events: int,
         number_of_threads: int,
         dry_run: bool,
         ):
    """
    Main function to drive the Geant4 simulation.

    Constructs the simulation environment using `ComponentAccumulator`:
    1. Initializes Detector Construction (DetectorConstruction_v1).
    2. Sets up the Event Reader to stream events from input.
    3. Configures the CaloHitBuilder to collect energy deposits.
    4. Merges components and executes the run loop.

    Args:
        logging_level (str): Verbosity level.
        input_file (str | Path): Path to the input event file.
        output_file (str | Path): Path to the output hit file.
        pre_init (str): Python code to execute before initialization.
        pre_exec (str): Python code to execute before the run loop.
        post_exec (str): Python code to execute after the run loop.
        enable_magnetic_field (bool): Toggle for the detector magnetic field.
        enable_solenoid_field (bool): Toggle for the field confined to the solenoid volume.
        tile_atlas_geometry (bool): Build the ATLAS-like tile calorimeter.
        tile_atlas_cells (bool): Use the ATLAS tile cells (needs tile_atlas_geometry).
        atlas_material_in_front (bool): Material in front of the barrel EM calorimeter as in ATLAS.
        tile_atlas_itc (bool): The gap between the tile barrel and the extended barrel as in ATLAS.
        atlas_emb (bool): Barrel EM calorimeter with the ATLAS absorber composition.
        atlas_hec (bool): HEC as in ATLAS (four samplings, nominal geometry and cells of ATLAS).
        atlas_emec (bool): EM end-cap as in ATLAS (two wheels, nominal geometry, composition and cells of ATLAS).
        atlas_endcap_cryostat (bool): Outer cylinders of the end-cap cryostats as in ATLAS.
        atlas_barrel_cryostat (bool): Outer part of the barrel cryostat as in ATLAS.
        tile_dual_readout (int): Two PMTs per tile cell: 0 = off, 1 = ATLAS U-shape, 2 = linear (needs tile_atlas_cells).
        active_energy_only (bool): Keep only the energy deposited in the active medium.
        birks_law (bool): Apply Birks' law in the scintillator and in the liquid argon.
        save_all_hits (bool): If True, saves all hits regardless of Region of Interest (RoI).
        free_running_hits (bool): Write the tile hits in the layout of the ATLAS HITS ntuple instead of the HIT stream.
        free_running_binning (bool): Time bins of the ATLAS simulation for the free-running hits.
        free_running_keep_hits (bool): Keep also the standard HIT stream (free-running hits in the tree FreeRunningTree).
        neutron_time_cut (float): Time (ns) after which neutrons are killed; 0 keeps the Geant4 default.
        timeout (int): Timeout in minutes.
        number_of_events (int): Number of events to process.
        number_of_threads (int): Number of Geant4 threads.
        dry_run (bool): If True, sets up but does not execute the run.
    """

    if isinstance(input_file, Path):
        input_file = str(input_file)
    if isinstance(output_file, Path):
        output_file = str(output_file)

    
    outputLevel = LoggingLevel.toC(logging_level)
    exec(pre_init)

    acc = ComponentAccumulator("ComponentAccumulator", 
                               DetectorConstruction_v1( "ATLAS", UseMagneticField=enable_magnetic_field,
                                                        UseSolenoidField=enable_solenoid_field,
                                                        TileAtlasGeometry=tile_atlas_geometry,
                                                        TileAtlasCells=tile_atlas_cells,
                                                        AtlasMaterialInFront=atlas_material_in_front,
                                                        TileAtlasItc=tile_atlas_itc,
                                                        AtlasEmb=atlas_emb,
                                                        AtlasHec=atlas_hec,
                                                        AtlasEmec=atlas_emec,
                                                        AtlasEndcapCryostat=atlas_endcap_cryostat,
                                                        AtlasBarrelCryostat=atlas_barrel_cryostat),
                               NumberOfThreads=number_of_threads,
                               OutputFile=output_file,
                               Timeout=timeout * MINUTES,
                               NeutronTimeCut=neutron_time_cut)

    gun = EventReader("EventReader", input_file,
                      # outputs
                      OutputEventKey=recordable("Events"),
                      OutputTruthKey=recordable("Particles"),
                      OutputSeedKey=recordable("Seeds"),
                      )
    calorimeter = CaloHitBuilder("CaloHitBuilder",
                                 HistogramPath="Expert/Hits",
                                 OutputLevel=outputLevel,
                                 OutputHitsKey=recordable("Hits"),
                                 ActiveEnergyOnly=active_energy_only,
                                 BirksLaw=birks_law,
                                 TileDualReadout=tile_dual_readout,
                                 FreeRunningHits=free_running_hits,
                                 FreeRunningBinning=free_running_binning,
                                 FreeRunningNtupleName=("FreeRunningTree" if free_running_keep_hits else "CollectionTree"),
                                 InputEventKey=recordable("Events")
                                 )
    
    gun.merge(acc)
    calorimeter.merge(acc)
    if free_running_hits and not free_running_keep_hits:
        # only the free-running hits (CaloFreeRunningHitWriter, added by the CaloHitBuilder)
        HIT = None
    else:
        HIT = RootStreamHITMaker("RootStreamHITMaker",
                                 OutputLevel=outputLevel,
                                 OnlyRoI= not save_all_hits,
                                 InputHitsKey=recordable("Hits"),
                                 InputEventKey=recordable("Events"),
                                 InputTruthKey=recordable("Particles"),
                                 InputSeedsKey=recordable("Seeds"),
                                 )
    if HIT is not None:
        acc += HIT
    
    exec(pre_exec)
    if not dry_run:
        acc.run(number_of_events)
        exec(post_exec)



if __name__ == "__main__":
    parser=parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if Path(args.input_file).is_dir():
        raise IsADirectoryError(f"Input file '{args.input_file}' was expected to be a file, "
                                 "but it is a directory. Please provide a list of files instead.")
    
    args = update_args_from_file(args)
    print(f"input file: {args.input_file}")
    if args.enable_magnetic_field and args.enable_solenoid_field:
        parser.error("--enable-magnetic-field and --enable-solenoid-field cannot be used together")
    if args.tile_atlas_cells and not args.tile_atlas_geometry:
        parser.error("--tile-atlas-cells needs --tile-atlas-geometry")
    if args.tile_dual_readout and not args.tile_atlas_cells:
        parser.error("--tile-dual-readout needs --tile-atlas-cells")
    if args.tile_atlas_itc and not args.tile_atlas_geometry:
        parser.error("--tile-atlas-itc needs --tile-atlas-geometry")
    if args.free_running_hits:
        # free-running hits as in the ATLAS simulation: ATLAS tile cells, active medium, Birks' law, U-shape dual readout
        if not (args.tile_atlas_geometry and args.tile_atlas_cells):
            parser.error("--free-running-hits needs --tile-atlas-geometry and --tile-atlas-cells")
        if args.tile_dual_readout and args.tile_pmt_split != 'ushape':
            parser.error("--free-running-hits uses the U-shape sharing of the ATLAS simulation (--tile-pmt-split ushape)")
        args.active_energy_only = True
        args.birks_law = True
        args.tile_dual_readout = True
        args.tile_pmt_split = 'ushape'
        if args.neutron_time_cut is None:
            args.neutron_time_cut = 150.
    elif args.free_running_no_binning or args.free_running_keep_hits:
        parser.error("--free-running-no-binning and --free-running-keep-hits need --free-running-hits")
    print(f"output file: {args.output_file}")
    print(f"number of threads: {args.number_of_threads}")

    main(
             logging_level         = args.output_level,
             input_file            = args.input_file,
             output_file           = args.output_file,
             pre_init              = args.pre_init,
             pre_exec              = args.pre_exec,
             post_exec             = args.post_exec,
             enable_magnetic_field = args.enable_magnetic_field,
             enable_solenoid_field = args.enable_solenoid_field,
             tile_atlas_geometry   = args.tile_atlas_geometry,
             tile_atlas_cells      = args.tile_atlas_cells,
             atlas_material_in_front = args.atlas_material_in_front,
             tile_atlas_itc        = args.tile_atlas_itc,
             atlas_emb             = args.atlas_emb,
             atlas_hec             = args.atlas_hec,
             atlas_emec            = args.atlas_emec,
             atlas_endcap_cryostat = args.atlas_endcap_cryostat,
             atlas_barrel_cryostat = args.atlas_barrel_cryostat,
             tile_dual_readout     = (0 if not args.tile_dual_readout else
                                      (1 if args.tile_pmt_split == 'ushape' else 2)),
             active_energy_only    = args.active_energy_only,
             birks_law             = args.birks_law,
             save_all_hits         = args.save_all_hits,
             free_running_hits     = args.free_running_hits,
             free_running_binning  = not args.free_running_no_binning,
             free_running_keep_hits = args.free_running_keep_hits,
             neutron_time_cut      = args.neutron_time_cut or 0,
             timeout               = args.timeout,
             number_of_events      = args.number_of_events,
             number_of_threads     = args.number_of_threads,
             dry_run               = args.dry_run,
        )




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
                             "with the ATLAS layer radii. The cells keep the default eta x phi segmentation. "
                             "Off by default. Use the same option in digit_trf.py.")
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
         active_energy_only: bool,
         birks_law: bool,
         save_all_hits : bool,
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
        active_energy_only (bool): Keep only the energy deposited in the active medium.
        birks_law (bool): Apply Birks' law in the scintillator and in the liquid argon.
        save_all_hits (bool): If True, saves all hits regardless of Region of Interest (RoI).
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
                                                        TileAtlasGeometry=tile_atlas_geometry),
                               NumberOfThreads=number_of_threads,
                               OutputFile=output_file,
                               Timeout=timeout * MINUTES)

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
                                 BirksLaw=birks_law
                                 )
    
    gun.merge(acc)
    calorimeter.merge(acc)
    HIT = RootStreamHITMaker("RootStreamHITMaker",
                             OutputLevel=outputLevel,
                             OnlyRoI= not save_all_hits,
                             InputHitsKey=recordable("Hits"),
                             InputEventKey=recordable("Events"),
                             InputTruthKey=recordable("Particles"),
                             InputSeedsKey=recordable("Seeds"),
                             )
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
             active_energy_only    = args.active_energy_only,
             birks_law             = args.birks_law,
             save_all_hits         = args.save_all_hits,
             timeout               = args.timeout,
             number_of_events      = args.number_of_events,
             number_of_threads     = args.number_of_threads,
             dry_run               = args.dry_run,
        )




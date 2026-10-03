#!/usr/bin/env python3
import argparse
import sys
import os

from pathlib            import Path
from typing             import List
from CaloCellBuilder    import CaloCellBuilder
from GaugiKernel        import LoggingLevel, get_argparser_formatter
from GaugiKernel        import ComponentAccumulator
from RootStreamBuilder  import RootStreamHITReader, recordable
from RootStreamBuilder  import RootStreamESDMaker

from reco.reco_job import merge_args, update_args, create_parallel_job
from geometry import DetectorConstruction_v1


"""
Script: digit_trf.py
Purpose: Performs the digitization step in the simulation chain.
         Converts Geant4 energy hits into digital signals (cells) by simulating
         electronic logical pulses, noise, and cross-talk.
Usage:
    digit_trf.py -i input.HIT.root -o output.ESD.root
"""

def parse_args():
    """
    Parses command-line arguments for the digitization job.

    Returns:
        argparse.Namespace: Configuration arguments including logging level and execution hooks.
    """
    # create the top-level parser
    parser = argparse.ArgumentParser(
        description='',
        formatter_class=get_argparser_formatter(),
        add_help=False)

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

    parser.add_argument('--tile-atlas-geometry', action='store_true',
                        dest='tile_atlas_geometry', required=False,
                        help="Use the ATLAS-like tile calorimeter geometry (see simu_trf.py). It must match "
                             "the option used in the simulation, since the tile cells depend on it.")
    parser.add_argument('--atlas-emb', action='store_true',
                        dest='atlas_emb', required=False,
                        help="Barrel EM calorimeter with the ATLAS absorber composition (see simu_trf.py). It must "
                             "match the option used in the simulation, since the barrel cells depend on it.")
    parser.add_argument('--tile-atlas-cells', action='store_true',
                        dest='tile_atlas_cells', required=False,
                        help="Hits simulated with --tile-atlas-cells (ATLAS tile cells). Not supported yet: the "
                             "digitization only knows the eta x phi grid, so the job stops with an error.")

   
    parser = merge_args(parser)

    return parser


def main(events : List[int],
         logging_level: str,
         input_file: str | Path,
         output_file: str | Path,
         pre_init: str,
         pre_exec: str,
         post_exec: str,
         tile_atlas_geometry: bool = False,
         atlas_emb: bool = False,
        ):
    """
    Main function for the digitization process.

    Reads Hits from the input file, simulates the calorimeter readout electronics
    (CaloCellBuilder), and produces an Event Summary Data (ESD) file containing
    calorimeter cells.

    Args:
        events (List[int]): List of event indices to process.
        logging_level (str): Logging verbosity.
        input_file (str | Path): Path to input HIT file.
        output_file (str | Path): Path to output ESD file.
        pre_init (str): Hook for pre-initialization code.
        pre_exec (str): Hook for pre-execution code.
        post_exec (str): Hook for post-execution code.
        tile_atlas_geometry (bool): Use the ATLAS-like tile calorimeter geometry.
        atlas_emb (bool): Barrel EM calorimeter with the ATLAS absorber composition.
    """

    if isinstance(input_file, Path):
        input_file = str(input_file)
    if isinstance(output_file, Path):
        output_file = str(output_file)

    outputLevel = LoggingLevel.toC(logging_level)

    exec(pre_init)

    acc = ComponentAccumulator("ComponentAccumulator", output_file)

    # the reader must be first in sequence
    reader = RootStreamHITReader("HITReader",
                                 InputFile=input_file,
                                 OutputHitsKey=recordable("Hits"),
                                 OutputEventKey=recordable("Events"),
                                 OutputTruthKey=recordable("Particles"),
                                 OutputSeedsKey=recordable("Seeds"),
                                 OutputLevel=outputLevel,
                                 )

    reader.merge(acc)

    # digitalization!    
    calorimeter = CaloCellBuilder("CaloCellBuilder", 
                                  DetectorConstruction_v1("ATLAS", TileAtlasGeometry=tile_atlas_geometry, AtlasEmb=atlas_emb),
                                  HistogramPath="Expert/Cells",
                                  OutputLevel=outputLevel,
                                  InputHitsKey=recordable("Hits"),
                                  OutputCellsKey=recordable("Cells"),
                                  OutputTruthCellsKey=recordable("TruthCells"),
                                  InputEventKey=recordable("Events"),
    )
    calorimeter.merge(acc)

    ESD = RootStreamESDMaker("RootStreamESDMaker",
                             InputCellsKey=recordable("Cells"),
                             InputCellsTruthKey=recordable("TruthCells"),
                             InputEventKey=recordable("Events"),
                             InputTruthKey=recordable("Particles"),
                             InputSeedsKey=recordable("Seeds"),
                             OutputLevel=outputLevel)
    acc += ESD
    
    exec(pre_exec)
    acc.run(events)
    exec(post_exec)


    


if __name__ == "__main__":
    parser=parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    if args.tile_atlas_cells:
        parser.error("--tile-atlas-cells: the digitization of the ATLAS tile cells is not supported yet "
                     "(the cells exist only in the hit file). Do not digitize these hits with the default cells.")
    args = update_args(args)
    pool  = create_parallel_job(args)
    pool( main, 
         logging_level    = args.output_level,
         pre_init         = args.pre_init,
         pre_exec         = args.pre_exec,
         post_exec        = args.post_exec,
         tile_atlas_geometry = args.tile_atlas_geometry,
         atlas_emb        = args.atlas_emb,
         )

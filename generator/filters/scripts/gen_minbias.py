#!/usr/bin/env python3

import argparse
import random
import sys
import os

from math           import ceil
from typing         import List
from joblib         import Parallel, delayed
from evtgen         import Pythia8
from filters        import FixedRegion
from GenKernel      import EventTape
from GaugiKernel    import get_argparser_formatter
from GaugiKernel    import LoggingLevel

from reco.gen_job import merge_args, update_args, create_parallel_job


datapath    = os.environ["LORENZETTI_EVTGEN_DATA_DIR"]
PILEUP_FILE = f'{datapath}/minbias_config.cmnd'

def parse_args():

    parser = argparse.ArgumentParser(
        description='',
        formatter_class=get_argparser_formatter(),
        add_help=False)

   
    parser.add_argument('--run-number', action='store',
                        dest='run_number', required=False,
                        type=int, default=0,
                        help="The run number.")
    parser.add_argument('--output-level', action='store',
                        dest='output_level', required=False,
                        type=str, default="INFO",
                        help="The output level messenger.")
    parser.add_argument('--eta-max', action='store', 
                        dest='eta_max', required = False, 
                        type=float, default=3.2,
                        help = "Maximum Eta.")
    parser.add_argument('--pileup-per-bunch-crossing', action='store',
                        dest='pileup_per_bunch_crossing', required=False,
                        type=float, default=-1,
                        help="The fixed number of pileup for each bunch crossing.")
    parser.add_argument('--bc-id-start', action='store',
                        dest='bc_id_start', required=False,
                        type=int, default=-21,
                        help="The bunch crossing id start.")
    parser.add_argument('--bc-id-end', action='store',
                        dest='bc_id_end', required=False,
                        type=int, default=4,
                        help="The bunch crossing id end.")
    parser.add_argument('--bc-duration', action='store',
                        dest='bc_duration', required=False,
                        type=int, default=25,
                        help="The bunch crossing duration (in nanoseconds).")
    parser.add_argument('--pileup-file', action='store',
                        dest='pileup_file', required=False,
                        type=str, default=PILEUP_FILE,
                        help="The pythia pileup file configuration.")
    parser.add_argument('--pt-min-charged', action='store',
                        dest='pt_min_charged', required=False,
                        type=float, default=0.7,
                        help="Minimum pT (GeV) of the charged particles passed to the simulation. "
                             "The cut stands in for the solenoid field when the field is not simulated; "
                             "use 0 (as ATLAS, which applies no particle cut) when the field is simulated.")
    parser.add_argument('--pt-min-neutral', action='store',
                        dest='pt_min_neutral', required=False,
                        type=float, default=0.05,
                        help="Minimum pT (GeV) of the neutral particles passed to the simulation. "
                             "Use 0 to keep every particle, as ATLAS does.")
    parser.add_argument('--sigma-z', action='store',
                        dest='sigma_z', required=False,
                        type=float, default=56,
                        help="Gaussian spread (RMS, in mm) of the collision vertex along the beam axis.")
    parser.add_argument('--sigma-t', action='store',
                        dest='sigma_t', required=False,
                        type=float, default=200,
                        help="Gaussian spread (RMS, in ps) of the collision time.")
    parser.add_argument('--hard-filter', action='store',
                        dest='hard_filter', required=False,
                        type=str, default='none', choices=['none', 'mark', 'keep'],
                        help="Selection of the ATLAS high-pT minimum bias sample (DSID 800831: truth jet "
                             "R = 0.6 with pT > 35 GeV, or prompt photon pT > 8 GeV, or B hadron pT > 5 GeV with "
                             "a lepton pT > 4 GeV). 'mark' flags each event (branch p_hard) without changing the "
                             "generated sequence; 'keep' writes only the events that pass (branch n_tried counts "
                             "the events generated).")
    parser.add_argument('--truth', action='store_true',
                        dest='truth', required=False,
                        help="Also write the full generator record of each event (all particles, status, "
                             "parents and vertices) to a HepMC3 ROOT file next to the output "
                             "(X.EVT.root -> X.TRUTH.root).")
    
    return merge_args(parser)



def main(events: List[int],
         logging_level: str,
         output_file: str,
         run_number: int,
         seed: int,
         eta_max: float,
         pileup_per_bunch_crossing : float,
         mb_file: str,
         bc_id_start: int,
         bc_id_end: int,
         pt_min_charged: float,
         pt_min_neutral: float,
         sigma_z: float,
         sigma_t: float,
         hard_filter: str,
         truth: bool):

    outputLevel = LoggingLevel.toC(logging_level)

    # one truth file per job (X.EVT.<job>.root -> X.EVT.<job>.TRUTH.root), joined after the jobs
    truth_file = output_file[:-len('.root')] + '.TRUTH.root' if truth else ""

    tape = EventTape("EventTape", OutputFile=output_file,
                     RunNumber=run_number)

    # dummy seed to store all pileup particles
    tape+=FixedRegion("Seed", Eta=0, Phi=0 )
  
    from filters import Pileup
    pileup = Pileup("Pileup",
                    Pythia8("MBGenerator", 
                            File=mb_file,
                            Seed=seed),
                    EtaMax=eta_max,
                    Select=2,
                    PileupAvg=0,
                    PileupSigma=0,
                    PileupPerBunch=pileup_per_bunch_crossing,
                    BunchIdStart=bc_id_start,
                    BunchIdEnd=bc_id_end,
                    OutputLevel=outputLevel,
                    DeltaEta=999,
                    DeltaPhi=999,
                    PtMinCharged=pt_min_charged,
                    PtMinNeutral=pt_min_neutral,
                    VertexSigmaZ=sigma_z,
                    VertexSigmaT=sigma_t,
                    HardFilter={'none': 0, 'mark': 1, 'keep': 2}[hard_filter],
                    TruthFile=truth_file,
    )
    tape += pileup
    tape.run(events)



if __name__ == "__main__":
    parser=parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    args = update_args(args)    
    job  = create_parallel_job(args)
    job(main, 
        logging_level=args.output_level,
        run_number=args.run_number,
        eta_max=args.eta_max,
        pileup_per_bunch_crossing=args.pileup_per_bunch_crossing,
        mb_file=args.pileup_file,
        bc_id_start=args.bc_id_start,
        bc_id_end=args.bc_id_end,
        pt_min_charged=args.pt_min_charged,
        pt_min_neutral=args.pt_min_neutral,
        sigma_z=args.sigma_z,
        sigma_t=args.sigma_t,
        hard_filter=args.hard_filter,
        truth=args.truth
    )

    if args.truth and not args.dry_run:
        # join the truth files of the jobs into X.TRUTH.root (X.EVT.root -> X.TRUTH.root)
        import glob, shutil, subprocess
        base  = args.output_file[:-len('.root')]
        final = (base[:-len('.EVT')] if base.endswith('.EVT') else base) + '.TRUTH.root'
        parts = sorted( glob.glob(base + '.*.TRUTH.root'), key=lambda f: int(f.split('.')[-3]) )
        if len(parts) == 1:
            shutil.move(parts[0], final)
        elif len(parts) > 1:
            subprocess.run(['hadd', '-f', final] + parts, check=True)
            for f in parts: os.remove(f)

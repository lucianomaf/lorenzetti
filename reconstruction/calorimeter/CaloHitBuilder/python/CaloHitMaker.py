

__all__ = ["CaloHitMaker"]


from GaugiKernel import Cpp, LoggingLevel
from GaugiKernel.macros import *
import ROOT


class CaloHitMaker( Cpp ):


  def __init__( self, name, sampling,
                OutputCollectionKey  : str    = "Hits",
                OutputLevel          : int    = LoggingLevel.toC('INFO'),
                DetailedHistograms   : bool   = False,
                HistogramPath        : str    = "/CaloHitMaker",
                SamplingNoiseStd     : float  = 0,
                ActiveEnergyOnly     : bool   = False,
                BirksLaw             : bool   = False,
                TileDualReadout      : int    = 0,
                FreeRunningHits      : bool   = False,
                FreeRunningListKey   : str    = "FreeRunningHits",
                FreeRunningBinning   : bool   = True,
              ):
                    
    Cpp.__init__(self, ROOT.CaloHitMaker(name) )
    self.Tools = []
    self.setProperty( "OutputCollectionKey"     , OutputCollectionKey         )
    self.setProperty( "EtaBins"                 , sampling.sensitive().EtaBins)
    self.setProperty( "PhiBins"                 , sampling.sensitive().PhiBins)
    self.setProperty( "RMin"                    , sampling.volume().RMin      )
    self.setProperty( "RMax"                    , sampling.volume().RMax      )
    self.setProperty( "ZMin"                    , sampling.volume().ZMin      )
    self.setProperty( "ZMax"                    , sampling.volume().ZMax      )
    self.setProperty( "Sampling"                , sampling.Sampling           )
    self.setProperty( "Segment"                 , sampling.sensitive().Segment)
    self.setProperty( "Detector"                , sampling.Detector           )
    self.setProperty( "BunchIdStart"            , sampling.BunchIdStart       )
    self.setProperty( "BunchIdEnd"              , sampling.BunchIdEnd         )
    self.setProperty( "BunchDuration"           , 25                          )
    self.setProperty( "SamplingNoiseStd"        , SamplingNoiseStd            )
    # Keep only the energy deposited in the active medium (liquid argon, scintillator)
    self.setProperty( "ActiveEnergyOnly"        , ActiveEnergyOnly            )
    # Apply Birks' law in the scintillator and in the liquid argon (ATLAS constants)
    self.setProperty( "BirksLaw"                , BirksLaw                    )
    # Dual readout of the tile cells (two PMTs per cell): 0 = off, 1 = ATLAS U-shape, 2 = linear sharing
    self.setProperty( "TileDualReadout"         , TileDualReadout             )
    # Free-running hits (simu_trf.py --free-running-hits): every deposit of the ATLAS tile cells also goes, per PMT, to a
    # list written in the layout of the ATLAS HITS ntuple (CaloFreeRunningHitWriter); with FreeRunningBinning the deposits
    # of a PMT are summed in the time bins of the ATLAS simulation (TileSimHit). Only the tile samplings with the ATLAS cells.
    self.setProperty( "FreeRunningHits"         , FreeRunningHits             )
    self.setProperty( "FreeRunningListKey"      , FreeRunningListKey          )
    self.setProperty( "FreeRunningBinning"      , FreeRunningBinning          )
    # Cells given as (r, z) boxes instead of the eta x phi grid (ATLAS tile cells, geometry/python/v1/TILE.py).
    # Without them (the default) the properties stay empty and the eta x phi grid is used. Only set when present:
    # setProperty cannot convert an empty list.
    cells = getattr( sampling.sensitive(), "Cells", None )
    if cells:
      self.setProperty( "CellEta"               , [float(x) for x in cells["Eta"]]      )
      self.setProperty( "CellDeltaEta"          , [float(x) for x in cells["DeltaEta"]] )
      if "BoxRMin" in cells:
        self.setProperty( "CellBoxRMin"         , [float(x) for x in cells["BoxRMin"]]  )
        self.setProperty( "CellBoxRMax"         , [float(x) for x in cells["BoxRMax"]]  )
        self.setProperty( "CellBoxZMin"         , [float(x) for x in cells["BoxZMin"]]  )
        self.setProperty( "CellBoxZMax"         , [float(x) for x in cells["BoxZMax"]]  )
        self.setProperty( "CellBoxIndex"        , [int(x) for x in cells["BoxCell"]]    )
      # ATLAS-like EM end-cap (geometry/python/v1/EMEC.py): the compartment of this maker, its eta index, the tables
      # of the compartment boundaries and the volumes of its wheel (see CaloHitMaker::emecFindCell)
      if "EmecCompartment" in cells:
        self.setProperty( "EmecCompartment"     , int(cells["EmecCompartment"])           )
        self.setProperty( "EmecEtaScale"        , float(cells["EmecEtaScale"])            )
        self.setProperty( "EmecEtaOffset"       , float(cells["EmecEtaOffset"])           )
        self.setProperty( "EmecMaxEta"          , int(cells["EmecMaxEta"])                )
        self.setProperty( "EmecFocalShift"      , float(cells["EmecFocalShift"])          )
        for key in ("EmecZSep12", "EmecZSep23", "EmecZInner"):   # integers, 1e-4 cm
          self.setProperty( key                 , [int(x) for x in cells[key]]             )
        for key in ("EmecWheelRMin", "EmecWheelRMax", "EmecWheelZMin", "EmecWheelZMax"):
          self.setProperty( key                 , [float(x) for x in cells[key]]           )
      # ATLAS-like barrel EM (geometry/python/v1/ECAL.py): sampling and region of the dictionary read by this maker
      # (see CaloHitMaker::embFindCell)
      if "EmbSampling" in cells:
        self.setProperty( "EmbMode"             , int(cells["EmbMode"])                   )
        self.setProperty( "EmbSampling"         , int(cells["EmbSampling"])               )
        self.setProperty( "EmbRegion"           , int(cells["EmbRegion"])                 )
      # ATLAS-like barrel presampler, cell from the gap of the electrode (EmbMode = 3): radius of the middle of the active layer
      if "PsbR0" in cells:
        self.setProperty( "PsbR0"               , float(cells["PsbR0"])                   )
      # radius on the axis of the module (ATLAS-like HEC, geometry/python/v1/HEC.py); absent: r
      if "RadiusModules" in cells:
        self.setProperty( "CellRadiusModules"   , int(cells["RadiusModules"])           )
      if "AtlasSection" in cells:
        self.setProperty( "CellAtlasSection"    , [int(x) for x in cells["AtlasSection"]]  )
        self.setProperty( "CellAtlasTower"      , [int(x) for x in cells["AtlasTower"]]    )
        self.setProperty( "CellAtlasSampling"   , [int(x) for x in cells["AtlasSampling"]] )
        self.setProperty( "CellRCentre"         , [float(x) for x in cells["RCentre"]]     )
        self.setProperty( "CellZCentre"         , [float(x) for x in cells["ZCentre"]]     )
    self.setProperty( "DetailedHistograms"      , DetailedHistograms          )
    self.setProperty( "HistogramPath"           , HistogramPath               )
    self.setProperty( "OutputLevel"             , OutputLevel                 )

 

  def core(self):
    # Attach all tools before return the core
    for tool in self.Tools:
      self._core.push_back(tool.core())
    return self._core


  def __add__( self, tool ):
    self.Tools += tool
    return self
  


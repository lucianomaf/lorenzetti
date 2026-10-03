__all__ = [
           "DetectorConstruction_v1", 
           ]


import ROOT

from typing import List
from prettytable import PrettyTable
from pprint import pprint
from tqdm import tqdm

from GaugiKernel.constants import *
from GaugiKernel import Cpp
from GaugiKernel.macros import *

from geometry.v1.PhysicalVolume   import Plates
from geometry.v1.ECAL             import getLArBarrelCfg
from geometry.v1.TILE             import getTileBarrelCfg, getTileExtendedCfg
from geometry.v1.EMEC             import getLArEMECCfg
from geometry.v1.HEC              import getHECCfg
from geometry.v1.DeadMaterials    import getDMVolumesCfg, getCrackVolumesCfg, getAtlasEndcapCryostatCfg
#from geometry.detectors.Tracking      import *


def flatten(nested_list : List) -> List:
  """Flatten a nested list."""
  result = []
  for item in nested_list:
    if isinstance(item, list):
      result.extend(flatten(item))
    else:
      result.append(item)
  return result


class DetectorConstruction_v1( Cpp ):

  def __init__( self, 
                name              : str, 
                UseMagneticField  : bool=False,
                UseSolenoidField  : bool=False,
                CutOnPhi          : bool=False,
                TileAtlasGeometry : bool=False,
                TileAtlasCells    : bool=False,
                AtlasMaterialInFront : bool=False,
                TileAtlasItc      : bool=False,
                AtlasEmb          : bool=False,
                AtlasEndcapCryostat : bool=False,
              ):

    Cpp.__init__(self, ROOT.DetectorConstruction_v1(name) )

    self.setProperty( "UseMagneticField", UseMagneticField  )
    # 2 T field only inside the ATLAS central solenoid volume (R < 1.23 m, |z| < 2.9 m);
    # cannot be combined with UseMagneticField (2 T in the whole world).
    self.setProperty( "UseSolenoidField", UseSolenoidField  )
    self.setProperty( "CutOnPhi"        , CutOnPhi          )

    self.samplings = []
    self.volumes = []
    # Center
    
    #volumes.extend( getPixelBarrelCfg()   )
    # Barrel EM calorimeter with the ATLAS absorber composition (lead 1.53/1.13 mm, steel, glue and electrode) and depth
    # (see geometry/python/v1/ECAL.py). The same value must be used in simulation and digitization.
    self.samplings.extend( getLArBarrelCfg(atlas_emb=AtlasEmb)   )
    # ATLAS-like tile calorimeter (tiles normal to the beam line, 18 mm period, ATLAS layer radii and z extent;
    # see geometry/python/v1/TILE.py). The same value must be used in simulation and digitization. With it the
    # dead material next to the tile calorimeter (inner aluminium shell, ITC block) follows the ATLAS z extent.
    # ATLAS tile cells (A1-A16, BC1-BC8, B9, D0-D6) instead of the eta x phi grid; needs TileAtlasGeometry.
    if TileAtlasCells and not TileAtlasGeometry:
      raise ValueError("TileAtlasCells needs TileAtlasGeometry.")
    self.samplings.extend( getTileBarrelCfg(atlas_geometry=TileAtlasGeometry, atlas_cells=TileAtlasCells)  )
    # Material in front of the barrel electromagnetic calorimeter as in ATLAS (solenoid, cryostat wall and the material
    # between the presampler and the accordion; see geometry/python/v1/DeadMaterials.py). Simulation only.
    self.volumes.extend( getDMVolumesCfg(tile_atlas_geometry=TileAtlasGeometry, atlas_material_in_front=AtlasMaterialInFront) )
    # Right side (A)
    self.samplings.extend( getTileExtendedCfg(atlas_geometry=TileAtlasGeometry, atlas_cells=TileAtlasCells)    )
    self.samplings.extend( getLArEMECCfg()         ) 
    self.samplings.extend( getHECCfg()             )
    # The gap between the tile barrel and the extended barrel as in ATLAS (plug of the ITC and services) instead of the
    # aluminium block; needs TileAtlasGeometry (see geometry/python/v1/DeadMaterials.py). Simulation only.
    self.volumes.extend( getCrackVolumesCfg(tile_atlas_geometry=TileAtlasGeometry, tile_atlas_itc=TileAtlasItc) )
    # Left side (B)
    self.samplings.extend( getTileExtendedCfg(left_side = True, atlas_geometry=TileAtlasGeometry, atlas_cells=TileAtlasCells) )
    self.samplings.extend( getLArEMECCfg(left_side=True)        ) 
    self.samplings.extend( getHECCfg(left_side=True)            )    
    self.volumes.extend( getCrackVolumesCfg(left_side=True, tile_atlas_geometry=TileAtlasGeometry, tile_atlas_itc=TileAtlasItc) )
    # Outer cylinders of the end-cap cryostats as in ATLAS (warm and cold vessels and the liquid argon between the EMEC
    # and the cold vessel; see geometry/python/v1/DeadMaterials.py). Simulation only.
    if AtlasEndcapCryostat:
      self.volumes.extend( getAtlasEndcapCryostatCfg() )
      self.volumes.extend( getAtlasEndcapCryostatCfg(left_side=True) )
    self.samplings = flatten(self.samplings)
    
  
  def compile(self):
    # Create all volumes inside of the detector
    
    volumes = [pv for pv in self.volumes]
    volumes.extend( [samp.volume() for samp in self.samplings] )
          
    for pv in tqdm( volumes, desc="Compiling  volumes...", ncols=70):
      self._core.AddVolume( pv.Name, pv.Plates, pv.AbsorberMaterial, pv.GapMaterial, 
                             # layer
                             pv.NofLayers, 
                             pv.AbsorberThickness, 
                             pv.GapThickness,
                             pv.LayerClearance,
                             # dimensions
                             pv.RMin, pv.RMax, pv.ZSize, 
                             pv.X, pv.Y, pv.Z,
                             # production cuts
                             pv.Cuts.ElectronCut, 
                             pv.Cuts.PositronCut, 
                             pv.Cuts.GammaCut, 
                             pv.Cuts.PhotonCut
                             )
      if pv.AbsorberSplitZ > 0:
        self._core.SetAbsorberSplit( pv.Name, pv.AbsorberSplitZ, pv.AbsorberMaterial2, pv.AbsorberThickness2 )

  def summary(self):

      print('Display all calorimeter samplings...')
      t = PrettyTable(["Name", "Plates", "z",'Zmin','Zmax', "Rmin", 
                       "Rmax", "abso","gap", "deta", "dphi", "EtaMin", 
                       "EtaMax", "N_bins", "Container"])
      # Add all volumes that came from a sampling detector and has a sensitive parameter
      for samp in self.samplings:
        pv = samp.volume(); sv = samp.sensitive()
        t.add_row( [pv.Name,
                    Plates.tostring(pv.Plates),pv.ZSize,pv.ZMin,pv.ZMax,pv.RMin,pv.RMax,
                    pv.AbsorberMaterial,pv.GapMaterial,
                    round(sv.DeltaEta,4) ,
                    round(sv.DeltaPhi,4) ,
                    sv.EtaMin,sv.EtaMax, 
                    len(sv.EtaBins)*len(sv.PhiBins) if sv.DeltaEta else len(sv.ZBins)*len(sv.PhiBins)  ,
                    samp.CollectionKey
                  ])
      print(t)
      print('Display all non-sensitive volumes...')
      t = PrettyTable(["Name", "Plates", "z",'Zmin','Zmax', "Rmin", 
                       "Rmax", "abso","gap"])
      # Add ither volumes that not came from a sampling detector (extra volumes only)
      for pv in self.volumes:
        t.add_row([pv.Name, Plates.tostring(pv.Plates),pv.ZSize, pv.ZMin, pv.ZMax, 
                   pv.RMin, pv.RMax, pv.AbsorberMaterial, pv.GapMaterial]) 
      print(t)

  
  
  def get_ui_commands(self) -> List[str]:

    commands = [
     "/vis/open OGL 600x600-0+0"
    ,"/vis/viewer/set/autoRefresh false"
    ,"/vis/verbose errors"
    ,"/vis/drawVolume"
    ,"/vis/viewer/set/viewpointVector 1 0 0"
    ,"/vis/viewer/set/lightsVector 1 0 0"
    ,"/vis/viewer/set/style wireframe"
    ,"/vis/viewer/set/auxiliaryEdge true"
    ,"/vis/viewer/set/lineSegmentsPerCircle 100"
    ,"/vis/scene/add/trajectories smooth"
    ,"/vis/modeling/trajectories/create/drawByCharge"
    ,"/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true"
    ,"/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2"
    ,"/vis/scene/endOfEventAction accumulate"
    ,"/vis/geometry/set/visibility World 0 false"
    ,"#/vis/geometry/set/visibility World 0 true"
    ,"/vis/ogl/set/displayListLimit 10000000"
    ]
  
    def _add_volume_vis_commands(name, color, visualization) -> List[str]:
      vis_command = [
         f"/vis/geometry/set/colour {name} 0 {color}"
        ,f"/vis/geometry/set/colour {name}_Layer 0 {color}"
        ,f"/vis/geometry/set/colour {name}_Abso 0 {color}"
        ,f"/vis/geometry/set/colour {name}_Gap 0 {color}"
        ,f"/vis/geometry/set/visibility {name} 0 {visualization}"
        ,f"/vis/geometry/set/visibility {name}_Layer 0 {visualization}"
        ,f"/vis/geometry/set/visibility {name}_Abso 0 {visualization}"
        ,f"/vis/geometry/set/visibility {name}_Gap 0 {visualization}"
      ]
      return vis_command

    for samp in self.samplings:
      vol = samp.volume()
      commands.extend( _add_volume_vis_commands( vol.name(), vol.Color, 'true' if vol.Visualization else 'false') )
    for vol in self.volumes:
      commands.extend( _add_volume_vis_commands( vol.name(), vol.Color, 'true' if vol.Visualization else 'false') )

    commands.extend( [
             "/vis/viewer/set/autoRefresh true"
            ,"/vis/verbose warnings"
            ])

    return commands




if __name__ == "__main__":
    atlas = DetectorConstruction_v1("ATLAS")
    atlas.summary()
    #atlas.compile()
    #pprint(atlas.create_visualization_commands())






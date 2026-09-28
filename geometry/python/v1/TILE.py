
__all__ = ["getTileBarrelCfg","getTileExtendedCfg"]

from GaugiKernel.constants import m,cm,mm,MeV,pi
from CaloCell.CaloDefs import Detector, CaloSampling

from .Calorimeter import Calorimeter
from .PhysicalVolume import PhysicalVolume, Plates
from .SensitiveDetector import SensitiveDetector

import os


# ATLAS-like tile calorimeter (option TileAtlasGeometry, off by default).
# In ATLAS the scintillating tiles are placed radially and normal to the beam line, stacked along z
# in periods of 18 mm: 5 mm master plate + 4 mm spacer plate + 5 mm master plate (14 mm steel), a
# 3 mm tile and 1 mm left by the 3 mm tile in its 4 mm pocket, steel:scintillator 4.67:1
# (TILECAL, hep-ex/9904032, sec. 2; JINST 3 (2008) S08003, sec. 5.3.1.2). Here each layer is built
# with vertical plates (Fe 14 mm + scintillator 3 mm + 1 mm empty clearance). The staggering of the
# tiles between radial rows is not modelled, nor the support girder at the outer radius.
# Layer radii (A, BC, D in the barrel; A, B, D in the extended barrel) from the ATLAS cell geometry;
# with this period they give 1.45/4.07/1.84 interaction lengths in the barrel and 1.45/2.62/3.30 in
# the extended barrel at normal incidence, against about 1.5/4.1/1.8 and 1.5/2.6/3.3 quoted in
# JINST 3 (2008) S08003 (sec. 5.3.1.4 and chap. 1). The cells keep the eta x phi segmentation below.
TILE_ATLAS_ABSORBER  = 14*mm
TILE_ATLAS_TILE      = 3*mm
TILE_ATLAS_CLEARANCE = 1*mm
TILE_ATLAS_BARREL_RADII   = [2300*mm, 2600*mm, 3440*mm, 3820*mm]
TILE_ATLAS_EXTENDED_RADII = [2300*mm, 2600*mm, 3140*mm, 3820*mm]


def _getAtlasTileLayers(radii, zsize):
    """
    Layout of the three ATLAS-like tile layers: vertical plates with the ATLAS period, as many whole
    periods as fit in zsize (the remainder is left empty, half at each end).
    """
    period = TILE_ATLAS_ABSORBER + TILE_ATLAS_TILE + TILE_ATLAS_CLEARANCE
    nperiods = int(zsize // period)
    return [ dict( Plates            = Plates.Vertical,
                   NofLayers         = nperiods,
                   AbsorberThickness = TILE_ATLAS_ABSORBER,
                   GapThickness      = TILE_ATLAS_TILE,
                   LayerClearance    = TILE_ATLAS_CLEARANCE,
                   RMin              = radii[i],
                   RMax              = radii[i+1] ) for i in range(3) ]


def getTileBarrelCfg(atlas_geometry=False):
    """
    Defines the geometry and readout configuration for the Tile Calorimeter (TileCal) Barrel.

    Constructs the physical volumes (Iron absorber + Scintillator gap) for the
    three longitudinal layers of the Tile Barrel and assigns readout parameters.

    Args:
        atlas_geometry (bool): If True, builds the ATLAS-like layout (tiles normal to the beam line,
                               18 mm period, ATLAS layer radii; see TILE_ATLAS_* above).

    Returns:
        List[Calorimeter]: A list of configured Calorimeter detector objects for the Tile Barrel.
    """

    basepath = os.environ['LORENZETTI_GEOMETRY_DATA_DIR']

    # TileCal 
    endcap_start = 3704.*mm
    gap_between_extended_barrel = 68.*cm
    tile_barrel_start = 0*m
    tile_barrel_end = (endcap_start - gap_between_extended_barrel)
    tile_barrel_z = (tile_barrel_start + tile_barrel_end) * 2

    if atlas_geometry:
        layers = _getAtlasTileLayers( TILE_ATLAS_BARREL_RADII, tile_barrel_z )
    else:
        r1 = 228.3*cm + 4*(6.0*cm + 4.0*cm)
        r2 = r1 + 11*(6.2*cm + 3.8*cm)
        r3 = r2 + 5*(6.2*cm + 3.8*cm)
        layers = [ dict( Plates=Plates.Horizontal, NofLayers=4 , AbsorberThickness=6.0*cm, GapThickness=4.0*cm, RMin=228.3*cm, RMax=r1 ),
                   dict( Plates=Plates.Horizontal, NofLayers=11, AbsorberThickness=6.2*cm, GapThickness=3.8*cm, RMin=r1      , RMax=r2 ),
                   dict( Plates=Plates.Horizontal, NofLayers=5 , AbsorberThickness=6.2*cm, GapThickness=3.8*cm, RMin=r2      , RMax=r3 ) ]


    tilecal1_pv =  PhysicalVolume( Name               = "TILE::TileCal1",
                                   AbsorberMaterial   = "G4_Fe", # absorber
                                   GapMaterial        = "PLASTIC SCINTILLATOR", # gap
                                   **layers[0], # plates, layers, thicknesses and radii
                                   ZSize              = tile_barrel_z,
                                   X=0,Y=0,Z=0, # x,y,z (center in 0,0,0)
                                   Visualization = True,
                                   Color         = 'salmon'
                                   )

    tilecal2_pv =  PhysicalVolume( Name               = "TILE::TileCal2",
                                   AbsorberMaterial   = "G4_Fe", # absorber
                                   GapMaterial        = "PLASTIC SCINTILLATOR", # gap
                                   **layers[1], # plates, layers, thicknesses and radii
                                   ZSize              = tile_barrel_z,
                                   X=0,Y=0,Z=0, # x,y,z (center in 0,0,0)
                                   Visualization = True,
                                   Color         = 'violetred'
                                   )

    tilecal3_pv =  PhysicalVolume( Name               = "TILE::TileCal3",
                                   AbsorberMaterial   = "G4_Fe", # absorber
                                   GapMaterial        = "PLASTIC SCINTILLATOR", # gap
                                   **layers[2], # plates, layers, thicknesses and radii
                                   ZSize              = tile_barrel_z,
                                   X=0,Y=0,Z=0, # x,y,z (center in 0,0,0)
                                   Visualization = True,
                                   Color         = 'orangered'
                                   )


    #tilecal1_pv.Cuts = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
    #tilecal2_pv.Cuts = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
    #tilecal3_pv.Cuts = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)


    tilecal1_sv = SensitiveDetector( tilecal1_pv, DeltaEta = 0.1  , DeltaPhi = pi/32 )
    tilecal2_sv = SensitiveDetector( tilecal2_pv, DeltaEta = 0.1  , DeltaPhi = pi/32 )
    tilecal3_sv = SensitiveDetector( tilecal3_pv, DeltaEta = 0.2  , DeltaPhi = pi/32 )




    # Configure the electronic frontend and the detector parameters
    tilecal1_det = Calorimeter( tilecal1_sv, -6, 4, -3, # sensitive volume, bunch start, bunch end, sampling start,
                                CollectionKey   = "Collection_TileCal1", # collection key
                                Detector        = Detector.TILE, # detector type
                                Sampling        = CaloSampling.TileCal1, # sampling type
                                Shaper          = basepath + "/pulseTile.dat", # pulse shaper
                                Noise           = 20*MeV, # electronic noise
                                Samples         = 7, # how many samples
                                OFWeightsEnergy = [-0.3781,   -0.3572,    0.1808,    0.8125,   0.2767,    -0.2056,    -0.3292], # optimal filter parameters
                                OFWeightsTime   = [0.197277, -2.336666, -18.474749, -0.247639, 14.832119, 4.627806, 1.401852], # values get from COOL, Field 1. 
                              )

    # Configure the electronic frontend and the detector parameters
    tilecal2_det = Calorimeter( tilecal2_sv, -6, 4, -3, # sensitive volume, bunch start, bunch end, sampling start,
                                CollectionKey = "Collection_TileCal2", # collection key
                                Detector      = Detector.TILE, # detector type
                                Sampling      = CaloSampling.TileCal2, # sampling type
                                Shaper        = basepath + "/pulseTile.dat", # pulse shaper
                                Noise         = 20*MeV, # electronic noise
                                Samples       = 7, # how many samples
                                OFWeightsEnergy = [-0.3781,   -0.3572,    0.1808,    0.8125,   0.2767,    -0.2056,    -0.3292], # optimal filter parameters
                                OFWeightsTime   = [0.197277, -2.336666, -18.474749, -0.247639, 14.832119, 4.627806, 1.401852], # values get from COOL, Field 1. 
                              )

    # Configure the electronic frontend and the detector parameters
    tilecal3_det = Calorimeter( tilecal3_sv, -6, 4, -3, # sensitive volume, bunch start, bunch end, sampling start,
                                CollectionKey = "Collection_TileCal3", # collection key
                                Detector      = Detector.TILE, # detector type
                                Sampling      = CaloSampling.TileCal3, # sampling type
                                Shaper        = basepath + "/pulseTile.dat", # pulse shaper
                                Noise         = 20*MeV, # electronic noise
                                Samples       = 7, # how many samples
                                OFWeightsEnergy = [-0.3781,   -0.3572,    0.1808,    0.8125,   0.2767,    -0.2056,    -0.3292], # optimal filter parameters
                                OFWeightsTime   = [0.197277, -2.336666, -18.474749, -0.247639, 14.832119, 4.627806, 1.401852], # values get from COOL, Field 1. 
                              )



    return [tilecal1_det, tilecal2_det, tilecal3_det]






def getTileExtendedCfg(left_side=False, atlas_geometry=False):
    """
    Defines the geometry and readout configuration for the Tile Calorimeter Extended Barrel.

    Constructs the physical volumes and assigns readout parameters for the Extended Barrel.
    Since the detector is symmetric but has distinct physical volumes for A-side and C-side,
    this function allows configuring either side.

    Args:
        left_side (bool): If True, configures the C-side (negative z). 
                          If False, configures the A-side (positive z).
        atlas_geometry (bool): If True, builds the ATLAS-like layout (tiles normal to the beam line,
                               18 mm period, ATLAS layer radii; see TILE_ATLAS_* above).

    Returns:
        List[Calorimeter]: A list of configured Calorimeter detector objects for the Tile Extended Barrel.
    """


    # TileExt Barrel
    endcap_start = 3.704*m
    extended_barrel_start = endcap_start
    extended_barrel_zsize = 2.83*m
    basepath = os.environ['LORENZETTI_GEOMETRY_DATA_DIR']
    sign = -1 if left_side else 1
    side_name = 'B' if left_side else 'A'

    if atlas_geometry:
        layers = _getAtlasTileLayers( TILE_ATLAS_EXTENDED_RADII, extended_barrel_zsize )
    else:
        r1 = 228.3*cm + 4*(6.0*cm + 4.0*cm)
        r2 = r1 + 11*(6.2*cm + 3.8*cm)
        r3 = r2 + 5*(6.2*cm + 3.8*cm)
        layers = [ dict( Plates=Plates.Horizontal, NofLayers=4 , AbsorberThickness=6.0*cm, GapThickness=4.0*cm, RMin=228.3*cm, RMax=r1 ),
                   dict( Plates=Plates.Horizontal, NofLayers=11, AbsorberThickness=6.2*cm, GapThickness=3.8*cm, RMin=r1      , RMax=r2 ),
                   dict( Plates=Plates.Horizontal, NofLayers=5 , AbsorberThickness=6.2*cm, GapThickness=3.8*cm, RMin=r2      , RMax=r3 ) ]


    tilecalExt1_pv =  PhysicalVolume( Name               = "TILE::TileCalExt1::"+side_name,
                                      AbsorberMaterial   = "G4_Fe", # absorber
                                      GapMaterial        = "PLASTIC SCINTILLATOR", # gap
                                      **layers[0], # plates, layers, thicknesses and radii
                                      ZSize              = extended_barrel_zsize  ,
                                      X=0,Y=0,Z=sign*(extended_barrel_start + 0.5*extended_barrel_zsize), # x,y,z 
                                      Visualization = True,
                                      Color         = 'salmon',
                                    )
    tilecalExt2_pv =  PhysicalVolume( Name               = "TILE::TileCalExt2::"+side_name,
                                      AbsorberMaterial   = "G4_Fe", # absorber
                                      GapMaterial        = "PLASTIC SCINTILLATOR", # gap
                                      **layers[1], # plates, layers, thicknesses and radii
                                      ZSize              = extended_barrel_zsize ,
                                      X=0,Y=0,Z=sign*(extended_barrel_start + 0.5*extended_barrel_zsize), # x,y,z 
                                      Visualization = True,
                                      Color         = 'violetred'
                                    )

    tilecalExt3_pv =  PhysicalVolume( Name               = "TILE::TileCalExt3::"+side_name,
                                      AbsorberMaterial   = "G4_Fe", # absorber
                                      GapMaterial        = "PLASTIC SCINTILLATOR", # gap
                                      **layers[2], # plates, layers, thicknesses and radii
                                      ZSize              = extended_barrel_zsize ,
                                      X=0,Y=0,Z=sign*(extended_barrel_start + 0.5*extended_barrel_zsize), # x,y,z 
                                      Visualization = True,
                                      Color         = 'orangered'
                                    )


    tilecalExt1_sv = SensitiveDetector( tilecalExt1_pv, DeltaEta = 0.1  , DeltaPhi = pi/32  )
    tilecalExt2_sv = SensitiveDetector( tilecalExt2_pv, DeltaEta = 0.1  , DeltaPhi = pi/32  )
    tilecalExt3_sv = SensitiveDetector( tilecalExt3_pv, DeltaEta = 0.2  , DeltaPhi = pi/32  )


    # Configure the electronic frontend and the detector parameters
    tilecalExt1_det = Calorimeter( tilecalExt1_sv, -6, 4, -3, # sensitive volume, bunch start, bunch end, sampling start,
                                   CollectionKey   = "Collection_TileExt1_"+side_name, # collection key
                                   Detector        = Detector.TILE, # detector type
                                   Sampling        = CaloSampling.TileExt1, # sampling type
                                   Shaper          = basepath + "/pulseTile.dat", # pulse shaper
                                   Noise           = 20*MeV, # electronic noise
                                   Samples         = 7, # how many samples
                                   OFWeightsEnergy = [-0.3781,   -0.3572,    0.1808,    0.8125,   0.2767,    -0.2056,    -0.3292], # optimal filter parameters
                                   OFWeightsTime   = [0.197277, -2.336666, -18.474749, -0.247639, 14.832119, 4.627806, 1.401852], # values get from COOL, Field 1. 
                                 )

    # Configure the electronic frontend and the detector parameters
    tilecalExt2_det = Calorimeter( tilecalExt2_sv, -6, 4, -3, # sensitive volume, bunch start, bunch end, sampling start,
                                   CollectionKey   = "Collection_TileExt2_"+side_name, # collection key
                                   Detector        = Detector.TILE, # detector type
                                   Sampling        = CaloSampling.TileExt2, # sampling type
                                   Shaper          = basepath + "/pulseTile.dat", # pulse shaper
                                   Noise           = 20*MeV, # electronic noise
                                   Samples         = 7, # how many samples
                                   OFWeightsEnergy = [-0.3781,   -0.3572,    0.1808,    0.8125,   0.2767,    -0.2056,    -0.3292], # optimal filter parameters
                                   OFWeightsTime   = [0.197277, -2.336666, -18.474749, -0.247639, 14.832119, 4.627806, 1.401852], # values get from COOL, Field 1. 
                                 )

    # Configure the electronic frontend and the detector parameters
    tilecalExt3_det = Calorimeter( tilecalExt3_sv, -6, 4, -3, # sensitive volume, bunch start, bunch end, sampling start,
                                   CollectionKey   = "Collection_TileExt3_"+side_name, # collection key
                                   Detector        = Detector.TILE, # detector type
                                   Sampling        = CaloSampling.TileExt3, # sampling type
                                   Shaper          = basepath + "/pulseTile.dat", # pulse shaper
                                   Noise           = 20*MeV, # electronic noise
                                   Samples         = 7, # how many samples
                                   OFWeightsEnergy = [-0.3781,   -0.3572,    0.1808,    0.8125,   0.2767,    -0.2056,    -0.3292], # optimal filter parameters
                                   OFWeightsTime   = [0.197277, -2.336666, -18.474749, -0.247639, 14.832119, 4.627806, 1.401852], # values get from COOL, Field 1. 
                                 )
    
    return [tilecalExt1_det, tilecalExt2_det, tilecalExt3_det]


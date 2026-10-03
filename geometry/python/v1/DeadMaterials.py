
__all__ = ["getCrackVolumesCfg", "getDMVolumesCfg"]

from GaugiKernel.constants import m,cm,mm
from .PhysicalVolume import PhysicalVolume, Plates, ProductionCuts
from .TILE import TILE_ATLAS_BARREL_HALF_Z, TILE_ATLAS_EXTENDED_Z_START


def getCrackVolumesCfg(left_side=False, tile_atlas_geometry=False, tile_atlas_itc=False):
    """
    Dead material in the crack between the barrel and the endcaps.

    Args:
        left_side (bool): If True, the C-side (negative z).
        tile_atlas_geometry (bool): If True, the aluminium block that stands for the ITC fills the gap of the
                                    ATLAS-like tile calorimeter (TILE_ATLAS_BARREL_HALF_Z to
                                    TILE_ATLAS_EXTENDED_Z_START) instead of the default gap.
        tile_atlas_itc (bool): If True (needs tile_atlas_geometry), the gap between the tile barrel and the extended
                               barrel is built as in ATLAS instead of the aluminium block: the plug of the ITC (D4 and
                               C10, passive) and the cables and services (see getAtlasItcCfg).
    """

    sign = -1 if left_side else 1
    side_name = 'B' if left_side else 'A'

    endcap_start = 3704.*mm
    gap_between_extended_barrel = 68.*cm
    crack_material_size = gap_between_extended_barrel
    crack_start = endcap_start - gap_between_extended_barrel
    if tile_atlas_geometry:
        crack_start = TILE_ATLAS_BARREL_HALF_Z
        crack_material_size = TILE_ATLAS_EXTENDED_Z_START - TILE_ATLAS_BARREL_HALF_Z
    crack_em_start = 3400*mm + 35*mm/2



    nlayers=2; absorber=4.3*cm; gap=10.7*cm; zsize=nlayers*(absorber+gap)
    crack_em_pv =  PhysicalVolume( Name               = "DM::Crack::EM::"+side_name, 
                                    Plates             = Plates.Vertical, # Logical type
                                    AbsorberMaterial   = "G4_Al", # absorber
                                    GapMaterial        = "liquidArgon", # gap
                                    NofLayers          = nlayers, # layers
                                    AbsorberThickness  = absorber, # abso
                                    GapThickness       = gap, # gap
                                    RMin               = 900*mm, # radio min,
                                    RMax               = 2232*mm, # radio max 
                                    ZSize              = zsize,
                                    X=0,Y=0,Z=sign * (crack_em_start + zsize/2),
                                    Visualization = True,
                                    Color         = 'gray'
                                    )


    nlayers=1; absorber=200*cm; gap=1*mm; rsize=nlayers*(absorber+gap)
    crack_tile_pv =  PhysicalVolume( Name               = "DM::Crack::TILE::"+side_name, 
                                     Plates             = Plates.Horizontal, # Logical type
                                     AbsorberMaterial   = "G4_Al", # absorber
                                     GapMaterial        = "Vacuum", # gap
                                     NofLayers          = nlayers, # layers
                                     AbsorberThickness  = absorber, # abso
                                     GapThickness       = gap, # gap
                                     RMin               = 228.3*cm, #crack_em_pv.RMax, #228.3*cm, # radio min,
                                     RMax               = 228.3*cm + rsize, # radio max 
                                     ZSize              = crack_material_size,
                                     X=0,Y=0,Z=sign * (crack_start + crack_material_size/2),
                                     Visualization = True,
                                     Color         = 'gray'
                                     )


    crack_em_pv.Cuts     = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
    crack_tile_pv.Cuts   = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)

    if tile_atlas_itc:
        if not tile_atlas_geometry:
            raise ValueError("tile_atlas_itc needs tile_atlas_geometry.")
        return [crack_em_pv] + getAtlasItcCfg(left_side=left_side)
    return [crack_em_pv, crack_tile_pv]


# The gap between the tile barrel (|z| < 2808 mm) and the extended barrel (|z| > 3554 mm) as in ATLAS. The gap "is filled
# with cables and services for the inner detector as well as power supplies and services for the barrel liquid-argon
# calorimeter"; "at the outer radius of the detector, a reduced section of a standard tile-calorimeter sub-module, the
# plug, provides additional coverage" (JINST 3 (2008) S08003, sec. 5.5). Boxes from fig. 5.12 of the same paper:
# D4 z 3225-3535 mm, r 3461-3850 mm; C10 z 3444-3531 mm, r 2988-3457 mm.
# - Plug: tiles normal to the beam line with the period of the ATLAS-like tile calorimeter (14 mm steel, 3 mm
#   scintillator, 1 mm clearance), passive (no cells): D4 with 17 periods (z 3227-3533 mm), C10 with 4 (3451.5-3523.5 mm).
# - Cables and services: the mass is not given by the paper; aluminium with 5% of the volume, from a fit of the
#   interaction lengths at the end of the tile calorimeter to fig. 5.2 of the paper for |eta| = 0.85 to 1.05 (the fit
#   gives 0 to 5%), in three boxes around the plug, in shells of 50 mm or less.
# Not modelled: the gap and cryostat scintillators E1-E4.
ATLAS_ITC_SERVICES_FRACTION = 0.05

def getAtlasItcCfg(left_side=False):
    sign = -1 if left_side else 1
    side_name = 'B' if left_side else 'A'
    period = 18*mm; steel = 14*mm; tile = 3*mm; clearance = 1*mm
    volumes = []
    for name, nper, zc, rmin, rmax in (("D4", 17, 3380*mm, 3461*mm, 3850*mm), ("C10", 4, 3487.5*mm, 2988*mm, 3457*mm)):
        pv = PhysicalVolume( Name               = "DM::ITC::"+name+"::"+side_name,
                             Plates             = Plates.Vertical,
                             AbsorberMaterial   = "G4_Fe",
                             GapMaterial        = "PLASTIC SCINTILLATOR",
                             NofLayers          = nper,
                             AbsorberThickness  = steel,
                             GapThickness       = tile,
                             LayerClearance     = clearance,
                             RMin               = rmin,
                             RMax               = rmax,
                             ZSize              = nper*period,
                             X=0,Y=0,Z=sign*zc,
                             Visualization = True,
                             Color         = 'gray' )
        pv.Cuts = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
        volumes.append(pv)
    for i, (z1, z2, r1, r2) in enumerate(((2808*mm, 3225*mm, 2283*mm, 3850*mm), (3225*mm, 3444*mm, 2283*mm, 3461*mm),
                                          (3444*mm, 3554*mm, 2283*mm, 2988*mm))):
        nlayers = max(1, int(round((r2 - r1) / (50*mm)))); layer = (r2 - r1) / nlayers
        pv = PhysicalVolume( Name               = "DM::ITC::Services%d::%s" % (i+1, side_name),
                             Plates             = Plates.Horizontal,
                             AbsorberMaterial   = "G4_Al",
                             GapMaterial        = "Vacuum",
                             NofLayers          = nlayers,
                             AbsorberThickness  = ATLAS_ITC_SERVICES_FRACTION*layer,
                             GapThickness       = (1 - ATLAS_ITC_SERVICES_FRACTION)*layer,
                             RMin               = r1,
                             RMax               = r2,
                             ZSize              = z2 - z1,
                             X=0,Y=0,Z=sign*(z1 + z2)/2,
                             Visualization = True,
                             Color         = 'gray' )
        pv.Cuts = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
        volumes.append(pv)
    return volumes


def getDMVolumesCfg(tile_atlas_geometry=False, atlas_material_in_front=False):
    """
    Dead material around the barrel calorimeters.

    Args:
        tile_atlas_geometry (bool): If True, the aluminium shell inside the tile barrel follows the length of the
                                    ATLAS-like tile barrel (|z| < TILE_ATLAS_BARREL_HALF_Z).
        atlas_material_in_front (bool): If True, the material in front of the barrel electromagnetic calorimeter
                                    follows ATLAS (see getAtlasMaterialInFrontCfg): the solenoid, the cryostat wall
                                    in front of the presampler and the material between the presampler and the
                                    accordion. Without it (the default) only 40 mm of aluminium stand in front of the
                                    presampler.
    """

    ecal_barrel_start = 0*m
    ecal_barrel_end   = 3.4*m
    ecal_barrel_z     = (ecal_barrel_start + ecal_barrel_end) * 2

    endcap_start = 3704.*mm
    gap_between_extended_barrel = 68.*cm

    tile_barrel_start = 0*m
    tile_barrel_end = (endcap_start - gap_between_extended_barrel)
    tile_barrel_z = (tile_barrel_start + tile_barrel_end) * 2
    if tile_atlas_geometry:
        tile_barrel_z = 2*TILE_ATLAS_BARREL_HALF_Z

    psb_rmin = 1460*mm

    nlayers=2; absorber=2*cm; gap=2*cm; rsize=nlayers*(absorber+gap)
    if atlas_material_in_front:
        # Same 80 mm envelope, with 2 x 33.8 = 67.6 mm of aluminium (0.76 X0) instead of 40 mm (see getAtlasMaterialInFrontCfg)
        absorber=ATLAS_CRYOSTAT_AL_THICKNESS/nlayers; gap=rsize/nlayers - absorber
    dm_pv =  PhysicalVolume( Name               = "DM::PS::Boundary", 
                             Plates             = Plates.Horizontal, # Logical type
                             AbsorberMaterial   = "G4_Al", # absorber
                             GapMaterial        = "Vacuum", # gap
                             NofLayers          = nlayers, # layers
                             AbsorberThickness  = absorber, # abso
                             GapThickness       = gap, # gap
                             RMin               = psb_rmin - rsize, #crack_em_pv.RMax, #228.3*cm, # radio min,
                             RMax               = psb_rmin, # radio max 
                             ZSize              = ecal_barrel_z,
                             X=0,Y=0,Z=0,
                             Visualization = True,
                             Color         = 'gray'
                             )


    nlayers=1; absorber=10*cm; gap=3*mm; rsize=nlayers*(absorber+gap)
    ecal_boundary_pv =   PhysicalVolume( Name               = "DM::LAr::Boundary", 
                                         Plates             = Plates.Horizontal, # Logical type
                                         AbsorberMaterial   = "G4_Al", # absorber
                                         GapMaterial        = "Vacuum", # gap
                                         NofLayers          = nlayers, # layers
                                         AbsorberThickness  = absorber, # abso
                                         GapThickness       = gap, # gap
                                         RMin               = 198*cm, #crack_em_pv.RMax, #228.3*cm, # radio min,
                                         RMax               = 198*cm + rsize, # radio max 
                                         ZSize              = ecal_barrel_z,
                                         X=0,Y=0,Z=0,
                                         Visualization = True,
                                         Color         = 'gray'
                                         )

    nlayers=1; absorber=10*cm; gap=3*mm; rsize=nlayers*(absorber+gap)
    tilecal_boundary_pv =   PhysicalVolume( Name               = "DM::TILE::Boundary", 
                                            Plates             = Plates.Horizontal, # Logical type
                                            AbsorberMaterial   = "G4_Al", # absorber
                                            GapMaterial        = "Vacuum", # gap
                                            NofLayers          = nlayers, # layers
                                            AbsorberThickness  = absorber, # abso
                                            GapThickness       = gap, # gap
                                            RMin               = 218*cm, #crack_em_pv.RMax, #228.3*cm, # radio min,
                                            RMax               = 218*cm + rsize, # radio max 
                                            ZSize              = tile_barrel_z,
                                            X=0,Y=0,Z=0,
                                            Visualization = True,
                                            Color         = 'gray'
                                            )


    dm_pv.Cuts                  = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
    ecal_boundary_pv.Cuts       = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
    tilecal_boundary_pv.Cuts    = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)

    volumes = [dm_pv, ecal_boundary_pv, tilecal_boundary_pv]
    if atlas_material_in_front:
        volumes.extend( getAtlasMaterialInFrontCfg() )
    return volumes


# Material in front of the barrel electromagnetic calorimeter as in ATLAS, in radiation lengths at normal incidence
# (X0 from PDG 2024: Al 88.97 mm, Cu 14.36 mm). Where the room is too short for aluminium alone, a copper plate is added
# so that the shell has the ATLAS X0 in the available space.
# - Cryostat wall in front of the presampler: 0.76 X0 = 67.6 mm of Al. ATLAS has 1.86 X0 in front of the presampler at
#   eta = 0 (JINST 3 (2008) S08003, fig. 5.1, top left), of which about 0.44 X0 are the inner detector (taken from the
#   Run 2 inner detector, ATLAS-TDR-030, fig. 2.6) and 0.66 X0 the solenoid (JINST 3 S08003, sec. 2.1.1). The cryostat
#   is made of two aluminium vessels (JINST 3 S08003, sec. 5.4); here its aluminium stays in the 80 mm envelope in front
#   of the presampler (r = 1380 to 1460 mm), outside the volume of the solenoid field.
# - Solenoid: 0.66 X0 (JINST 3 S08003, sec. 2.1.1) in the 50 mm of the coil (inner and outer diameter 2.46 and 2.56 m,
#   length 5.8 m; JINST 3 S08003, table 2.1): 1.68 mm of Cu and 48.32 mm of Al.
# - Between the presampler and the accordion: 0.60 X0, the 0.68 X0 of ATLAS between the presampler and the accordion
#   (fig. 5.1, top left, at eta = 0) less the 11 mm of liquid argon of the presampler; in the 28.9 mm between the
#   presampler (r = 1471.01 mm) and the first layer (r = 1500 mm): 4.71 mm of Cu and 24.19 mm of Al.
# Not modelled: the tapered walls of the cryostat, the end-cap cryostat walls and the services; the inner detector
# (a separate step).
ATLAS_CRYOSTAT_AL_THICKNESS = 67.6*mm

def getAtlasMaterialInFrontCfg():
    """
    Solenoid and material between the presampler and the accordion of the barrel, as in ATLAS (see above).
    The cryostat wall is set in getDMVolumesCfg (DM::PS::Boundary).
    """
    solenoid_pv = PhysicalVolume( Name               = "DM::Solenoid",
                                  Plates             = Plates.Horizontal,
                                  AbsorberMaterial   = "G4_Cu",
                                  GapMaterial        = "G4_Al",
                                  NofLayers          = 1,
                                  AbsorberThickness  = 1.68*mm,
                                  GapThickness       = 48.32*mm,
                                  RMin               = 1230.05*mm,   # just outside the solenoid field volume (r < 1230 mm)
                                  RMax               = 1280.05*mm,
                                  ZSize              = 2*2900*mm,
                                  X=0,Y=0,Z=0,
                                  Visualization = True,
                                  Color         = 'gray'
                                  )
    ps_accordion_pv = PhysicalVolume( Name               = "DM::PS::Accordion",
                                      Plates             = Plates.Horizontal,
                                      AbsorberMaterial   = "G4_Cu",
                                      GapMaterial        = "G4_Al",
                                      NofLayers          = 1,
                                      AbsorberThickness  = 4.71*mm,
                                      GapThickness       = 24.19*mm,
                                      RMin               = 1471.05*mm,   # the presampler ends at 1471.01 mm
                                      RMax               = 1499.95*mm,   # the first layer starts at 1500 mm
                                      ZSize              = 2*3.4*m,
                                      X=0,Y=0,Z=0,
                                      Visualization = True,
                                      Color         = 'gray'
                                      )
    solenoid_pv.Cuts     = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
    ps_accordion_pv.Cuts = ProductionCuts(ElectronCut = 1, PositronCut = 1, GammaCut = 1)
    return [solenoid_pv, ps_accordion_pv]

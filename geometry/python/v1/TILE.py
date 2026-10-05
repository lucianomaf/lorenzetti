
__all__ = ["getTileBarrelCfg","getTileExtendedCfg","getAtlasTileCells"]

from GaugiKernel.constants import m,cm,mm,MeV,pi
from CaloCell.CaloDefs import Detector, CaloSampling

from .Calorimeter import Calorimeter
from .PhysicalVolume import PhysicalVolume, Plates
from .SensitiveDetector import SensitiveDetector

import os
import math


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
# Longitudinal extent of the ATLAS-like tile calorimeter, from the ATLAS cell layout (JINST 3 (2008) S08003,
# fig. 5.12: long barrel up to |z| = 2816 mm, extended barrel from 3554 to 6113 mm), rounded to whole 18 mm
# periods so that no tile is cut: long barrel |z| < 2808 mm (312 periods), extended barrel 3554 < |z| < 6110 mm
# (142 periods). The aluminium block that stands for the ITC fills the space in between (see DeadMaterials).
TILE_ATLAS_BARREL_HALF_Z     = 2808*mm
TILE_ATLAS_EXTENDED_Z_START  = 3554*mm
TILE_ATLAS_EXTENDED_Z_END    = 6110*mm


# ATLAS tile cells (option TileAtlasCells, off by default; needs TileAtlasGeometry). In ATLAS a cell is the group of
# tiles read out by the same photomultiplier, with fixed z boundaries in each radial row (JINST 3 (2008) S08003,
# sec. 5.3.1.4, "the fibre grouping is used to define a three-dimensional cell structure", and fig. 5.12). The
# boundaries below were read from fig. 5.12 (vector drawing, z scale from its ruler; about 2 mm precision), in mm, on
# the positive side; the negative side is the mirror image. Each boundary is moved to the nearest boundary between two
# 18 mm periods, so that every tile belongs to a single cell, as in ATLAS.
# Long barrel: A1-A10, B1-B9 and C1-C8 read together as BC1-BC8 (B9 only in the B row), D0 (across z = 0) to D3.
# Extended barrel: A12-A16, B11-B15, D5 and D6. The gap/crack cells (C10, D4, E1-E4) are not built.
TILE_ATLAS_CELLS_BARREL_A = [271.9, 510.4, 767.8, 1029.0, 1305.3, 1585.5, 1865.6, 2195.0, 2531.9]  # A1/A2 ... A9/A10
TILE_ATLAS_CELLS_BARREL_B = [306.0, 593.7, 870.0, 1169.1, 1468.1, 1797.5, 2138.2, 2505.4]          # B1/B2 ... B8/B9
TILE_ATLAS_CELLS_BARREL_C = [340.0, 669.4, 998.7, 1339.4, 1687.7, 2058.7, 2456.2]                  # C1/C2 ... C7/C8
TILE_ATLAS_CELLS_BARREL_D = [389.2, 1127.4, 1903.5]                                                # D0 half length, D1/D2, D2/D3
TILE_ATLAS_CELLS_BC_RADIUS = 2990*mm   # B/C row boundary inside the BC layer
TILE_ATLAS_CELLS_EXT_A = [3720.6, 4174.8, 4682.1, 5234.8]   # A12/A13 ... A15/A16
TILE_ATLAS_CELLS_EXT_B = [3803.8, 4320.2, 4886.5, 5473.3]   # B11/B12 ... B14/B15
TILE_ATLAS_CELLS_EXT_D = [4742.7]                           # D5/D6


def _roundToPeriod(z, z0):
    """Nearest boundary between two 18 mm periods, for periods starting at z0."""
    period = TILE_ATLAS_ABSORBER + TILE_ATLAS_TILE + TILE_ATLAS_CLEARANCE
    return z0 + round((z - z0) / period) * period


def _cellEta(boxes):
    """Eta of the geometric centre of the cell (area centroid of its boxes in the r-z plane) and delta eta taken at
    the mean radius of the cell, over the full z extent of the cell."""
    area = sum((b[1]-b[0])*(b[3]-b[2]) for b in boxes)
    rc = sum((b[1]-b[0])*(b[3]-b[2])*0.5*(b[0]+b[1]) for b in boxes) / area
    zc = sum((b[1]-b[0])*(b[3]-b[2])*0.5*(b[2]+b[3]) for b in boxes) / area
    rmean = 0.5*(min(b[0] for b in boxes) + max(b[1] for b in boxes))
    zmin = min(b[2] for b in boxes); zmax = max(b[3] for b in boxes)
    return math.asinh(zc/rc), math.asinh(zmax/rmean) - math.asinh(zmin/rmean)


def _cellsFromRows(rows):
    """
    rows: list of (rmin, rmax, [(name, zmin, zmax), ...]); cells with the same name in different rows are one cell.
    Returns the table used by CaloHitMaker: one entry per box (r and z limits, cell index) and one per cell (name,
    eta, delta eta), cells in increasing z of their centre.
    """
    boxes = {}
    for rmin, rmax, cells in rows:
        for name, zmin, zmax in cells:
            boxes.setdefault(name, []).append((rmin, rmax, zmin, zmax))
    def zcentre(name):
        b = boxes[name]; area = sum((x[1]-x[0])*(x[3]-x[2]) for x in b)
        return sum((x[1]-x[0])*(x[3]-x[2])*0.5*(x[2]+x[3]) for x in b) / area
    names = sorted(boxes, key=zcentre)
    table = dict(Names=[], Eta=[], DeltaEta=[], BoxRMin=[], BoxRMax=[], BoxZMin=[], BoxZMax=[], BoxCell=[],
                 AtlasSection=[], AtlasTower=[], AtlasSampling=[], RCentre=[], ZCentre=[])
    for i, name in enumerate(names):
        eta, deta = _cellEta(boxes[name])
        table['Names'].append(name); table['Eta'].append(round(eta, 5)); table['DeltaEta'].append(round(deta, 5))
        # ATLAS identifier fields and area centroid in (r, z) of the cell (free-running hits)
        section, tower, sampling = atlasTileIdFields(name)
        b = boxes[name]; area = sum((x[1]-x[0])*(x[3]-x[2]) for x in b)
        table['AtlasSection'].append(section); table['AtlasTower'].append(tower); table['AtlasSampling'].append(sampling)
        table['RCentre'].append(sum((x[1]-x[0])*(x[3]-x[2])*0.5*(x[0]+x[1]) for x in b) / area)
        table['ZCentre'].append(zcentre(name))
        for rmin, rmax, zmin, zmax in boxes[name]:
            table['BoxRMin'].append(rmin); table['BoxRMax'].append(rmax)
            table['BoxZMin'].append(zmin); table['BoxZMax'].append(zmax); table['BoxCell'].append(i)
    return table


def _barrelRow(prefix, edges, zend, first=1, crosses_zero=False):
    """Cells of one barrel row, both sides. The side is in the name: '+' for z > 0, '-' for z < 0."""
    e = [_roundToPeriod(x, -TILE_ATLAS_BARREL_HALF_Z) for x in edges]
    cells = []
    if crosses_zero:   # D0: one cell from -e[0] to e[0], no side
        cells.append((f"{prefix}{first}", -e[0], e[0]))
        first += 1
        bounds = e + [zend]
    else:
        bounds = [0.0] + e + [zend]
    for k in range(len(bounds)-1):
        n = first + k
        cells.append((f"{prefix}{n}+", bounds[k], bounds[k+1]))
        cells.append((f"{prefix}{n}-", -bounds[k+1], -bounds[k]))
    return cells


def getAtlasTileCells(sampling_name, side=1):
    """
    ATLAS tile cells for one sampling of the ATLAS-like tile calorimeter, as the (r, z) boxes read by CaloHitMaker.
    sampling_name: 'TileCal1', 'TileCal2', 'TileCal3' (long barrel) or 'TileExt1', 'TileExt2', 'TileExt3';
    side: +1 or -1 for the extended barrel.
    """
    zb = TILE_ATLAS_BARREL_HALF_Z
    rb = TILE_ATLAS_BARREL_RADII
    re = TILE_ATLAS_EXTENDED_RADII
    if sampling_name == 'TileCal1':
        return _cellsFromRows([(rb[0], rb[1], _barrelRow('A', TILE_ATLAS_CELLS_BARREL_A, zb))])
    if sampling_name == 'TileCal2':
        brow = [(n.replace('B', 'BC', 1) if not n.startswith('B9') else n, z0, z1)
                for n, z0, z1 in _barrelRow('B', TILE_ATLAS_CELLS_BARREL_B, zb)]
        crow = [(n.replace('C', 'BC', 1), z0, z1) for n, z0, z1 in _barrelRow('C', TILE_ATLAS_CELLS_BARREL_C, zb)]
        return _cellsFromRows([(rb[1], TILE_ATLAS_CELLS_BC_RADIUS, brow), (TILE_ATLAS_CELLS_BC_RADIUS, rb[2], crow)])
    if sampling_name == 'TileCal3':
        return _cellsFromRows([(rb[2], rb[3], _barrelRow('D', TILE_ATLAS_CELLS_BARREL_D, zb, first=0, crosses_zero=True))])
    ext = {'TileExt1': ('A', TILE_ATLAS_CELLS_EXT_A, 12, re[0], re[1]),
           'TileExt2': ('B', TILE_ATLAS_CELLS_EXT_B, 11, re[1], re[2]),
           'TileExt3': ('D', TILE_ATLAS_CELLS_EXT_D, 5,  re[2], re[3])}
    prefix, edges, first, rmin, rmax = ext[sampling_name]
    z0, z1 = TILE_ATLAS_EXTENDED_Z_START, TILE_ATLAS_EXTENDED_Z_END
    bounds = [z0] + [_roundToPeriod(x, z0) for x in edges] + [z1]
    s = '+' if side > 0 else '-'
    cells = []
    for k in range(len(bounds)-1):
        lo, hi = (bounds[k], bounds[k+1]) if side > 0 else (-bounds[k+1], -bounds[k])
        cells.append((f"{prefix}{first+k}{s}", lo, hi))
    return _cellsFromRows([(rmin, rmax, cells)])


# ATLAS identifier fields (section, tower, sampling) of each ATLAS tile cell, used by the free-running hits of simu_trf.py.
# Fields and allowed combinations from the public ATLAS identifier dictionary, Athena (Apache 2.0),
# DetectorDescription/IdDictParser/data/IdDictTileCalorimeter.xml: section 1 = long barrel, 2 = extended barrel;
# sampling 0 = A, 1 = BC (B in the extended barrel), 2 = D; tower = the 0.1 eta tower of the cell. Allowed (section, tower):
# samplings, from the same dictionary: barrel tower 0: 0-2 (D0 only on side +1), 1, 3, 5, 7, 8: 0-1, 2, 4, 6: 0-2, 9: 0;
# extended barrel tower 10: 1-2, 11, 13, 14: 0-1, 12: 0-2, 15: 0.
_ATLAS_TILE_ALLOWED = {1: {0: (0, 1, 2), 1: (0, 1), 2: (0, 1, 2), 3: (0, 1), 4: (0, 1, 2), 5: (0, 1), 6: (0, 1, 2), 7: (0, 1),
                           8: (0, 1), 9: (0,)},
                       2: {10: (1, 2), 11: (0, 1), 12: (0, 1, 2), 13: (0, 1), 14: (0, 1), 15: (0,)}}

def atlasTileIdFields(name):
    """(section, tower, sampling) of the ATLAS identifier of an ATLAS tile cell, from its name (with or without the side
    sign): A1-A10, BC1-BC8, B9, D0-D3 in the long barrel; A12-A16, B11-B15, D5, D6 in the extended barrel."""
    n = name.rstrip('+-')
    letters = n.rstrip('0123456789'); number = int(n[len(letters):])
    if letters == 'A':
        fields = (1, number - 1, 0) if number <= 10 else (2, number - 1, 0)
    elif letters in ('BC', 'B'):
        fields = (1, number - 1, 1) if number <= 9 else (2, number - 1, 1)
    elif letters == 'D':
        fields = (1, 2*number, 2) if number <= 3 else (2, 10 + 2*(number - 5), 2)
    else:
        raise ValueError(f"not an ATLAS tile cell: {name}")
    section, tower, sampling = fields
    if sampling not in _ATLAS_TILE_ALLOWED.get(section, {}).get(tower, ()):
        raise ValueError(f"{name}: ({section}, {tower}, {sampling}) is not in the ATLAS identifier dictionary")
    return fields


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


def getTileBarrelCfg(atlas_geometry=False, atlas_cells=False):
    """
    Defines the geometry and readout configuration for the Tile Calorimeter (TileCal) Barrel.

    Constructs the physical volumes (Iron absorber + Scintillator gap) for the
    three longitudinal layers of the Tile Barrel and assigns readout parameters.

    Args:
        atlas_geometry (bool): If True, builds the ATLAS-like layout (tiles normal to the beam line,
                               18 mm period, ATLAS layer radii and z extent; see TILE_ATLAS_* above).
        atlas_cells (bool): If True (needs atlas_geometry), the cells are the ATLAS cells A1-A10, BC1-BC8, B9 and
                            D0-D3 (see TILE_ATLAS_CELLS_* above) instead of the eta x phi grid.

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
        tile_barrel_z = 2*TILE_ATLAS_BARREL_HALF_Z
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
    if atlas_cells:
        if not atlas_geometry:
            raise ValueError("The ATLAS tile cells (atlas_cells) need the ATLAS-like tile geometry (atlas_geometry).")
        tilecal1_sv.Cells = getAtlasTileCells("TileCal1")
        tilecal2_sv.Cells = getAtlasTileCells("TileCal2")
        tilecal3_sv.Cells = getAtlasTileCells("TileCal3")




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






def getTileExtendedCfg(left_side=False, atlas_geometry=False, atlas_cells=False):
    """
    Defines the geometry and readout configuration for the Tile Calorimeter Extended Barrel.

    Constructs the physical volumes and assigns readout parameters for the Extended Barrel.
    Since the detector is symmetric but has distinct physical volumes for A-side and C-side,
    this function allows configuring either side.

    Args:
        left_side (bool): If True, configures the C-side (negative z). 
                          If False, configures the A-side (positive z).
        atlas_geometry (bool): If True, builds the ATLAS-like layout (tiles normal to the beam line,
                               18 mm period, ATLAS layer radii and z extent; see TILE_ATLAS_* above).
        atlas_cells (bool): If True (needs atlas_geometry), the cells are the ATLAS cells A12-A16, B11-B15, D5 and
                            D6 (see TILE_ATLAS_CELLS_* above) instead of the eta x phi grid.

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
        extended_barrel_start = TILE_ATLAS_EXTENDED_Z_START
        extended_barrel_zsize = TILE_ATLAS_EXTENDED_Z_END - TILE_ATLAS_EXTENDED_Z_START
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
    if atlas_cells:
        if not atlas_geometry:
            raise ValueError("The ATLAS tile cells (atlas_cells) need the ATLAS-like tile geometry (atlas_geometry).")
        tilecalExt1_sv.Cells = getAtlasTileCells("TileExt1", sign)
        tilecalExt2_sv.Cells = getAtlasTileCells("TileExt2", sign)
        tilecalExt3_sv.Cells = getAtlasTileCells("TileExt3", sign)


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


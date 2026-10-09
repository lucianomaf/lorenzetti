
__all__ = ["getLArBarrelCfg", "getAtlasEmbVolumesCfg", "atlasEmbCells", "atlasEmbLayers"]


import os
import numpy as np
from CaloCell.CaloDefs import Detector, CaloSampling
from GaugiKernel.constants import m,cm,mm,MeV,pi

from .Calorimeter import Calorimeter
from .PhysicalVolume import PhysicalVolume, Plates
from .SensitiveDetector import SensitiveDetector


def getLArBarrelCfg(atlas_emb=False, atlas_emb_cells=False):
  """
  Defines the geometry and readout configuration for the Liquid Argon (LAr) Barrel Calorimeter.

  Args:
      atlas_emb (bool): If True, the three layers have the absorber composition and the depth of the ATLAS barrel
                        (see the comment before the ATLAS layers below); the presampler and the readout do not change.
      atlas_emb_cells (bool): Needs atlas_emb. The three layers end at |z| = 3165 mm as in ATLAS, with steps that follow
                        |eta| = 1.475 near the inner radius, and have the cells of the ATLAS identifier dictionary (see
                        _getAtlasEmbCellsCfg and the notes before it). Their volumes come from getAtlasEmbVolumesCfg.
                        The presampler does not change.

  Constructs the physical volumes (PreSampler, Back, Middle, Strips) and assigns
  readout parameters such as pulse shapes, noise levels, and Optimal Filter weights.

  Returns:
      List[Calorimeter]: A list of configured Calorimeter detector objects for the LAr Barrel.
  """


  if atlas_emb_cells and not atlas_emb:
    raise ValueError("atlas_emb_cells needs atlas_emb.")

  basepath = os.environ['LORENZETTI_GEOMETRY_DATA_DIR']

  # ECal
  ecal_barrel_start = 0*m
  ecal_barrel_end   = 3.4*m
  ecal_barrel_z     = (ecal_barrel_start + ecal_barrel_end) * 2


  psb_pv =  PhysicalVolume( Name               = "LAr::PSB",
                            Plates             = Plates.Horizontal, # Logical type
                            AbsorberMaterial   = "Vacuum", # absorber
                            GapMaterial        = "liquidArgon", # gap
                            NofLayers          = 1, # layers
                            AbsorberThickness  =  0.01*mm, # abso
                            GapThickness       = 1.1*cm, # gap
                            RMin               = 146*cm, # radio min,
                            RMax               = 146*cm + 1*(0.01*mm + 1.1*cm), # radio max 
                            ZSize              = ecal_barrel_z ,# z (3.4 left and 3.4 right)
                            X=0,Y=0,Z=0, # x,y,z (center in 0,0,0)
                            # Visualization
                            Visualization = True,
                            Color         = 'orange'
                          )

  emb1_pv = PhysicalVolume( Name               = "LAr::EMB1", 
                            Plates             = Plates.Horizontal, # Logical type
                            AbsorberMaterial   = "G4_Pb", # absorber
                            GapMaterial        = "liquidArgon", # gap
                            NofLayers          = 16, # layers
                            AbsorberThickness  = 1.51*mm, # abso
                            GapThickness       = 4.49*mm, # gap
                            RMin               = 150*cm, # radio min,
                            RMax               = 150*cm + 16*(1.51*mm + 4.49*mm), # radio max 
                            ZSize              = ecal_barrel_z ,# z (3.4 left and 3.4 right)
                            X=0,Y=0,Z=0, # x,y,z (center in 0,0,0)
                            Visualization = True,
                            Color         = 'aquamarine'
                          )

  emb2_pv = PhysicalVolume( Name               = "LAr::EMB2", 
                            Plates             = Plates.Horizontal, # Logical type
                            AbsorberMaterial   = "G4_Pb", # absorber
                            GapMaterial        = "liquidArgon", # gap
                            NofLayers          = 55, # layers
                            AbsorberThickness  = 1.7*mm, # abso
                            GapThickness       = 4.3*mm, # gap
                            RMin               = emb1_pv.RMax, # radio min,
                            RMax               = emb1_pv.RMax + 55*(1.7*mm + 4.3*mm), # radio max 
                            ZSize              = ecal_barrel_z ,# z (3.4 left and 3.4 right)
                            X=0,Y=0,Z=0, # x,y,z (center in 0,0,0)
                            Visualization = True,
                            Color         = 'cornflowerblue'
                          )

  emb3_pv = PhysicalVolume( Name               = "LAr::EMB3", 
                            Plates             = Plates.Horizontal, # Logical type
                            AbsorberMaterial   = "G4_Pb", # absorber
                            GapMaterial        = "liquidArgon", # gap
                            NofLayers          = 9, # layers
                            AbsorberThickness  = 1.7*mm, # abso
                            GapThickness       = 4.3*mm, # gap
                            RMin               = emb2_pv.RMax, # radio min,
                            RMax               = emb2_pv.RMax + 9*(1.7*mm + 4.3*mm), # radio max 
                            ZSize              = ecal_barrel_z ,# z (3.4 left and 3.4 right)
                            X=0,Y=0,Z=0, # x,y,z (center in 0,0,0)
                            Visualization = True,
                            Color         = 'cyan'
                          )



  # ATLAS-like barrel layers (option atlas_emb). Absorber of ATLAS (JINST 3 (2008) S08003, sec. 5.2.1): lead of 1.53 mm
  # for |eta| < 0.8 and 1.13 mm for |eta| > 0.8, two 0.2 mm stainless-steel sheets (as iron) and the glue and the readout
  # electrode as 0.832 mm of polyimide, fitted to the X0 of the accordion at eta = 0 in fig. 5.1 (22.4 X0); liquid argon of
  # 2 x 2.1 mm (sec. 5.2.2). Radial shells with a period of 6.962 mm (the accordion is not modelled); where the lead is
  # thinner the argon takes the room (the number of absorbers and the folds are fixed). The change of lead at |eta| = 0.8
  # is a cut in z at the mean radius of each layer, z = r * sinh(0.8) (about +-0.04 in eta at the inner and outer
  # radius). Layers as in the default geometry (14, 47 and 7 periods), ending at 1973.4 mm (470 mm of depth in fig. 5.4).
  if atlas_emb:
    period = 6.962*mm; abso153 = (1.53 + 0.4 + 0.832)*mm; abso113 = (1.13 + 0.4 + 0.832)*mm
    rmin = 150*cm; layers = []
    for name, nlayers, color in (("LAr::EMB1", 14, 'aquamarine'), ("LAr::EMB2", 47, 'cornflowerblue'), ("LAr::EMB3", 7, 'cyan')):
      rmax = rmin + nlayers*period
      layers.append( PhysicalVolume( Name               = name,
                                     Plates             = Plates.Horizontal,
                                     AbsorberMaterial   = "ATLAS_EMB_ABSORBER_153",
                                     GapMaterial        = "liquidArgon",
                                     NofLayers          = nlayers,
                                     AbsorberThickness  = abso153,
                                     GapThickness       = period - abso153,
                                     AbsorberSplitZ     = 0.5*(rmin + rmax)*np.sinh(0.8),
                                     AbsorberMaterial2  = "ATLAS_EMB_ABSORBER_113",
                                     AbsorberThickness2 = abso113,
                                     RMin               = rmin,
                                     RMax               = rmax,
                                     ZSize              = ecal_barrel_z,
                                     X=0,Y=0,Z=0,
                                     Visualization = True,
                                     Color         = color ) )
      rmin = rmax
    emb1_pv, emb2_pv, emb3_pv = layers

  psb_sv  = SensitiveDetector( psb_pv , DeltaEta = 0.025  , DeltaPhi = pi/32  )
  emb1_sv = SensitiveDetector( emb1_pv, DeltaEta = 0.00325, DeltaPhi = pi/32  )
  emb2_sv = SensitiveDetector( emb2_pv, DeltaEta = 0.025  , DeltaPhi = pi/128 )
  emb3_sv = SensitiveDetector( emb3_pv, DeltaEta = 0.050  , DeltaPhi = pi/128 )



  # Configure the electronic frontend and the detector parameters
  psb_det  = Calorimeter( psb_sv, -21, 3, -2, # sensitive volume, bunch start, bunch end, sampling start,
                          CollectionKey   = "Collection_PSB", # collection key
                          Detector        = Detector.TTEM, # detector type
                          Sampling        = CaloSampling.PSB, # sampling type
                          Shaper          = basepath + "/pulseLar.dat", # pulse shaper
                          Noise           = 90*MeV, # electronic noise
                          Samples         = 5, # how many samples
                          OFWeightsEnergy = [-0.0000853580,    0.265132,    0.594162,     0.389505,     0.124353], # optimal filter parameters for energy estimation
                          OFWeightsTime   = [-0.0000853580,   -12.870312690734863, -27.39136505126953, 8.075883865356445, 13.768877029418945] # EMB1 ATLAS sample to test (fake first number to add to 5) 
                        )
  # Configure the electronic frontend and the detector parameters
  emb1_det = Calorimeter( emb1_sv, -21, 3, -2, # sensitive volume, bunch start, bunch end, sampling start,
                          CollectionKey   = "Collection_EMB1", # collection key
                          Detector        = Detector.TTEM, # detector type
                          Sampling        = CaloSampling.EMB1, # sampling type
                          Shaper          = basepath + "/pulseLar.dat", # pulse shaper
                          Noise           = 26*MeV, # electronic noise
                          Samples         = 5, # how many samples
                          OFWeightsEnergy = [-0.0000853580,    0.265132,    0.594162,     0.389505,     0.124353], # optimal filter parameters for energy estimation
                          OFWeightsTime   = [-0.0000853580,   -12.870312690734863, -27.39136505126953, 8.075883865356445, 13.768877029418945] # EMB1 ATLAS sample to test (fake first number to add to 5) 
                        )
  # Configure the electronic frontend and the detector parameters
  emb2_det = Calorimeter( emb2_sv, -21, 3, -2, # sensitive volume, bunch start, bunch end, sampling start,
                          CollectionKey   = "Collection_EMB2", # collection key
                          Detector        = Detector.TTEM, # detector type
                          Sampling        = CaloSampling.EMB2, # sampling type
                          Shaper          = basepath + "/pulseLar.dat", # pulse shaper
                          Noise           = 60*MeV, # electronic noise
                          Samples         = 5, # how many samples
                          OFWeightsEnergy = [-0.0000853580,    0.265132,    0.594162,     0.389505,     0.124353], # optimal filter parameters for energy estimation
                          OFWeightsTime   = [-0.0000853580,   -12.870312690734863, -27.39136505126953, 8.075883865356445, 13.768877029418945] # EMB1 ATLAS sample to test (fake first number to add to 5) 
                        )
  # Configure the electronic frontend and the detector parameters
  emb3_det = Calorimeter( emb3_sv, -21, 3, -2, # sensitive volume, bunch start, bunch end, sampling start,
                          CollectionKey   = "Collection_EMB3", # collection key
                          Detector        = Detector.TTEM, # detector type
                          Sampling        = CaloSampling.EMB3, # sampling type
                          Shaper          = basepath + "/pulseLar.dat", # pulse shaper
                          Noise           = 40*MeV, # electronic noise
                          Samples         = 5, # how many samples
                          OFWeightsEnergy = [-0.0000853580,    0.265132,    0.594162,     0.389505,     0.124353], # optimal filter parameters for energy estimation
                          OFWeightsTime   = [-0.0000853580,   -12.870312690734863, -27.39136505126953, 8.075883865356445, 13.768877029418945] # EMB1 ATLAS sample to test (fake first number to add to 5) 
                        )


  if atlas_emb_cells:
    # ATLAS-like length and cells of the three layers (see _getAtlasEmbCellsCfg); the presampler is unchanged
    return [psb_det] + _getAtlasEmbCellsCfg()
  return [psb_det, emb1_det, emb2_det, emb3_det]



# ---------------------------------------------------------------------------------------------------------------------
# ATLAS-like length and cells of the barrel EM calorimeter (option atlas_emb_cells, needs atlas_emb; simu_trf.py and
# digit_trf.py --atlas-emb-cells; off by default). F41 conserto 17 / F45.
#
# - Length: the active EM barrel of ATLAS ends at |z| = 3165 mm (ATLAS geometry database, BarrelGeometry-03, ZMAX =
#   316.5 cm, read by the ATLAS geometry as LArEMBMotherZmax, LArGeoBarrel/src/BarrelConstruction.cxx; the ATLAS
#   simulation keeps 3150 mm as fiducial limit and notes that it "should be changed in database to become 3164",
#   LArG4Barrel/src/LArBarrelGeometry.cxx). Near the inner radius the edge follows |eta| = 1.475 (BarrelGeometry-03
#   ETACUT; RMINHIGHZ = 154.8 cm), approximated by steps of one layer (one period of 6.962 mm): where r sinh(1.475) is
#   below 3165 mm at the middle radius r of a layer, the layer ends at |z| = r sinh(1.475) (the first four layers of the
#   first sampling). The three layers keep the radii, the absorbers and the change of lead of atlas_emb (same cut in z).
# - Cells: the regions of the ATLAS identifier dictionary (IdDictLArCalorimeter_DC3-05-Comm-01.xml) in each half barrel:
#   EMB1 region 0, strips of 0.003125 with indices 1-447 (the strip 0 does not exist), 64 phi slices, and region 1, three
#   cells of 0.025 in 1.4-1.475, 256 phi slices; EMB2 region 0, 56 cells of 0.025 up to 1.4, and region 1, one cell of
#   0.075 in 1.4-1.475, 256 phi slices; EMB3, 27 cells of 0.05 up to 1.35, 256 phi slices. One hit maker per sampling,
#   region and side (Segment = region; side A z > 0, side B z < 0). The cell of a step comes from the pseudorapidity of
#   its point seen from the origin, with the index formulas of the ATLAS simulation (LArBarrelGeometry.cxx,
#   CalculateIdentifier), kept only within the dictionary (CaloHitMaker::embFindCell); the eta written in the hit is the
#   nominal centre of the cell, eta0 + (i - i_min + 0.5) delta eta. The radial boundaries between the samplings that follow
#   eta and the rule of ATLAS that sends the back sampling above |eta| = 1.325 to the middle one are not modelled: the
#   sampling is the layer of atlas_emb.
# - phi: slices of 2 pi / N on the -pi..pi grid of the Lorenzetti barrel (no half-absorber shift).
# Not changed by this option: the barrel presampler, the end-cap presampler and the material behind the end of the barrel.
#
EMB_ATLAS_Z_MAX   = 3165*mm      # BarrelGeometry-03 ZMAX 316.5 cm
EMB_ATLAS_ETA_CUT = 1.475        # BarrelGeometry-03 ETACUT
EMB_ATLAS_Z_EPS   = 1e-6*mm      # the hit makers of side A take z > EMB_ATLAS_Z_EPS, those of side B z <= -EMB_ATLAS_Z_EPS
EMB_ATLAS_RMIN    = 150*cm       # as atlas_emb
EMB_ATLAS_PERIOD  = 6.962*mm     # as atlas_emb
EMB_ATLAS_LAYERS  = ( ("LAr::EMB1", 14, 'aquamarine'), ("LAr::EMB2", 47, 'cornflowerblue'), ("LAr::EMB3", 7, 'cyan') )
#                     sampling region eta0      delta eta first last  phi slices   (identifier dictionary, per side)
EMB_ATLAS_REGIONS = [ (1,       0,     0.003125, 0.003125, 1,    447,  64 ),
                      (1,       1,     1.4,      0.025,    0,    2,    256),
                      (2,       0,     0.0,      0.025,    0,    55,   256),
                      (2,       1,     1.4,      0.075,    0,    0,    256),
                      (3,       0,     0.0,      0.05,     0,    26,   256) ]


def atlasEmbLayers():
    """(name, periods, colour, rmin, rmax, z of the change of lead) of the three layers of atlas_emb."""
    out = []; rmin = EMB_ATLAS_RMIN
    for name, n, color in EMB_ATLAS_LAYERS:
        rmax = rmin + n*EMB_ATLAS_PERIOD
        out.append( (name, n, color, rmin, rmax, 0.5*(rmin + rmax)*np.sinh(0.8)) )
        rmin = rmax
    return out


def atlasEmbLayerZ(r):
    """End in |z| of a period whose middle radius is r: 3165 mm, or r sinh(1.475) where this is smaller."""
    return min(EMB_ATLAS_Z_MAX, r*np.sinh(EMB_ATLAS_ETA_CUT))


def getAtlasEmbVolumesCfg():
    """
    Volumes of the three layers of the ATLAS-like barrel (see above): |z| up to 3165 mm, and the first periods of EMB1
    as steps that follow |eta| = 1.475. Built in place of the layers of atlas_emb; their cells come from the samplings.
    """
    abso153 = (1.53 + 0.4 + 0.832)*mm; abso113 = (1.13 + 0.4 + 0.832)*mm
    vols = []
    for name, n, color, rmin, rmax, zsplit in atlasEmbLayers():
        edges = [rmin + k*EMB_ATLAS_PERIOD for k in range(n)] + [rmax]
        groups = []   # [first period, |z| end, number of periods]
        for k in range(n):
            z = atlasEmbLayerZ(0.5*(edges[k] + edges[k+1]))
            if groups and groups[-1][1] == z: groups[-1][2] += 1
            else: groups.append([k, z, 1])
        for i, (k0, z, nk) in enumerate(groups):
            pvname = name if z == EMB_ATLAS_Z_MAX else f"{name}::Step{i+1}"
            vols.append( PhysicalVolume( Name               = pvname,
                                         Plates             = Plates.Horizontal,
                                         AbsorberMaterial   = "ATLAS_EMB_ABSORBER_153",
                                         GapMaterial        = "liquidArgon",
                                         NofLayers          = nk,
                                         AbsorberThickness  = abso153,
                                         GapThickness       = EMB_ATLAS_PERIOD - abso153,
                                         AbsorberSplitZ     = zsplit,
                                         AbsorberMaterial2  = "ATLAS_EMB_ABSORBER_113",
                                         AbsorberThickness2 = abso113,
                                         RMin               = edges[k0],
                                         RMax               = edges[k0 + nk],
                                         ZSize              = 2*z,
                                         X=0,Y=0,Z=0,
                                         Visualization = True,
                                         Color         = color ) )
    return vols


def atlasEmbCells(sampling, region, side=1):
    """
    Cell table of one region of one sampling of the ATLAS-like barrel for CaloHitMaker and CaloCellMaker: the nominal eta
    and delta eta of each cell (first to last index of the dictionary), the eta edges of the region and the sampling and
    region read by CaloHitMaker::embFindCell. side: +1 (A) or -1 (B). Returns (table, number of phi slices).
    """
    for s, reg, eta0, deta, first, last, nphi in EMB_ATLAS_REGIONS:
        if (s, reg) == (sampling, region): break
    else:
        raise RuntimeError(f"EMB{sampling} region {region}: not in the dictionary")
    table = dict(Names=[], Eta=[], DeltaEta=[], EmbSampling=sampling, EmbRegion=region)
    for i in range(first, last + 1):
        table['Names'].append(f"EMB{sampling}/{region}/{i}{'+' if side > 0 else '-'}")
        table['Eta'].append(side * (eta0 + deta * (i - first + 0.5))); table['DeltaEta'].append(deta)
    edges = [eta0 + deta * k for k in range(last - first + 2)]
    table['EtaEdges'] = edges if side > 0 else sorted(-e for e in edges)
    return table, nphi


def _getAtlasEmbCellsCfg():
    """Hit makers of the ATLAS-like barrel: [emb1, emb2, emb3], lists of Calorimeter, two sides each (see above)."""
    basepath = os.environ['LORENZETTI_GEOMETRY_DATA_DIR']
    noise = {1: 26*MeV, 2: 60*MeV, 3: 40*MeV}           # those of the default EMB1-3
    sampling_enum = {1: CaloSampling.EMB1, 2: CaloSampling.EMB2, 3: CaloSampling.EMB3}
    abso153 = (1.53 + 0.4 + 0.832)*mm
    out = []
    for s, (name, n, color, rmin, rmax, zsplit) in zip((1, 2, 3), atlasEmbLayers()):
        dets = []
        for sign, side_name in ((1, 'A'), (-1, 'B')):
            z0, z1 = (EMB_ATLAS_Z_EPS, EMB_ATLAS_Z_MAX) if sign > 0 else (-EMB_ATLAS_Z_MAX, -EMB_ATLAS_Z_EPS)
            # Limits of the hit makers of this sampling on this side (not built: the volumes come from
            # getAtlasEmbVolumesCfg); within them the cell comes from eta (CaloHitMaker::embFindCell).
            pv = PhysicalVolume( Name               = f"{name}::{side_name}",
                                 Plates             = Plates.Horizontal,
                                 AbsorberMaterial   = "ATLAS_EMB_ABSORBER_153",
                                 GapMaterial        = "liquidArgon",
                                 NofLayers          = n,
                                 AbsorberThickness  = abso153,
                                 GapThickness       = EMB_ATLAS_PERIOD - abso153,
                                 RMin               = rmin,
                                 RMax               = rmax,
                                 ZSize              = z1 - z0,
                                 X=0, Y=0, Z=0.5 * (z0 + z1),
                                 Visualization = False,
                                 Color         = color )
            pv.Envelope = True
            for s_, region, eta0, deta, first, last, nphi in EMB_ATLAS_REGIONS:
                if s_ != s: continue
                cells, nphi = atlasEmbCells(s, region, side=sign)
                edges = cells.pop('EtaEdges')
                lo, hi = sorted(abs(e) for e in (edges[0], edges[-1]))
                sv = SensitiveDetector( pv, EtaMin = sign*lo, EtaMax = sign*hi, Segment = region, DeltaEta = deta,
                                        DeltaPhi = 2*pi/nphi )
                if len(sv.PhiBins) - 1 != nphi:
                    raise RuntimeError(f"EMB{s} region {region}: {len(sv.PhiBins) - 1} phi slices, expected {nphi}")
                sv.EtaBins = edges; sv.EtaMin = edges[0]; sv.EtaMax = edges[-1]
                sv.Cells = cells
                dets.append( Calorimeter( sv, -21, 3, -2,
                                  CollectionKey   = f"Collection_EMB{s}_{region}" + side_name,
                                  Detector        = Detector.TTEM,
                                  Sampling        = sampling_enum[s],
                                  Shaper          = basepath + "/pulseLar.dat",
                                  Noise           = noise[s],
                                  Samples         = 5,
                                  OFWeightsEnergy = [-0.0000853580,    0.265132,    0.594162,     0.389505,     0.124353],
                                  OFWeightsTime   = [-0.0000853580,   -12.870312690734863, -27.39136505126953, 8.075883865356445, 13.768877029418945]
                                ) )
        out.append(dets)
    return out



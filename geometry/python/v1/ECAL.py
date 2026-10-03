
__all__ = ["getLArBarrelCfg"]


import os
import numpy as np
from CaloCell.CaloDefs import Detector, CaloSampling
from GaugiKernel.constants import m,cm,mm,MeV,pi

from .Calorimeter import Calorimeter
from .PhysicalVolume import PhysicalVolume, Plates
from .SensitiveDetector import SensitiveDetector


def getLArBarrelCfg(atlas_emb=False):
  """
  Defines the geometry and readout configuration for the Liquid Argon (LAr) Barrel Calorimeter.

  Args:
      atlas_emb (bool): If True, the three layers have the absorber composition and the depth of the ATLAS barrel
                        (see the comment before the ATLAS layers below); the presampler and the readout do not change.

  Constructs the physical volumes (PreSampler, Back, Middle, Strips) and assigns
  readout parameters such as pulse shapes, noise levels, and Optimal Filter weights.

  Returns:
      List[Calorimeter]: A list of configured Calorimeter detector objects for the LAr Barrel.
  """


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


  return [psb_det, emb1_det, emb2_det, emb3_det]



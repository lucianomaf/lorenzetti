"""
Material in front of the accordion of the barrel as in ATLAS, piece by piece (F41 conserto 17c, generated on 09/10/2026
by Simulacao_MinBias/41_consertos_lorenzetti/c17_emb_ps_atlas/c17c_ps_mae/gerar_frente_barril.py; do not edit by hand).

Sources: the public ATLAS geometry of the GeoModel project, geometry-ATLAS-R3S-2021-03-00-00.db
(https://geomodel.web.cern.ch/atlas-geometry-data/old/), for the cryostat walls, the solenoid, the materials and the
front electronics of the EM barrel; the open Athena code (LArGeoBarrel/src/BarrelPresamplerConstruction.cxx and
BarrelConstruction.cxx) for the presampler mother. Pieces that do not cover the whole phi range are layers with the
phi-averaged mixture (the rest is liquid argon, as in the ATLAS cryostat).
"""

# name: (density in g/cm3, {element symbol: mass fraction})
ATLAS_BARREL_FRONT_MIXTURES = {
    'ATLAS_STD_G10': (2.000003, {'C': 0.235440, 'H': 0.019758, 'F': 0.744802}),
    'ATLAS_PS_Shell': (1.900003, {'H': 0.081008, 'C': 0.551609, 'O': 0.367383}),
    'ATLAS_PS_Prepreg1': (1.900003, {'H': 0.081008, 'C': 0.551609, 'O': 0.367383}),
    'ATLAS_PS_Prepreg2': (1.898215, {'H': 0.080797, 'C': 0.550170, 'O': 0.366424, 'Ar': 0.002609}),
    'ATLAS_PS_Board': (1.989374, {'H': 0.057483, 'C': 0.391417, 'O': 0.260692, 'Cu': 0.130289, 'Ar': 0.160119}),
    'ATLAS_PS_Connectics_1': (1.524777, {'Cu': 0.023512, 'H': 0.018550, 'C': 0.126315, 'O': 0.084128, 'Ar': 0.747495}),
    'ATLAS_PS_Connectics_2': (1.617491, {'Cu': 0.070016, 'H': 0.020255, 'C': 0.137921, 'O': 0.091858, 'Ar': 0.679950}),
    'ATLAS_PS_Connectics_3': (1.743535, {'Cu': 0.125305, 'H': 0.022281, 'C': 0.151720, 'O': 0.101049, 'Ar': 0.599645}),
    'ATLAS_PS_Connectics_4': (1.902909, {'Cu': 0.184728, 'H': 0.024459, 'C': 0.166551, 'O': 0.110927, 'Ar': 0.513334}),
    'ATLAS_PS_Connectics_5': (2.095614, {'Cu': 0.244508, 'H': 0.026651, 'C': 0.181471, 'O': 0.120863, 'Ar': 0.426507}),
    'ATLAS_PS_Connectics_6': (2.321650, {'Cu': 0.301980, 'H': 0.028757, 'C': 0.195815, 'O': 0.130417, 'Ar': 0.343031}),
    'ATLAS_PS_Connectics_7': (2.581015, {'Cu': 0.355524, 'H': 0.030720, 'C': 0.209179, 'O': 0.139317, 'Ar': 0.265260}),
    'ATLAS_PS_Connectics_8': (2.873711, {'Cu': 0.404341, 'H': 0.032509, 'C': 0.221363, 'O': 0.147432, 'Ar': 0.194355}),
    'ATLAS_PS_Plate': (1.887689, {'H': 0.079545, 'C': 0.541642, 'O': 0.360745, 'Ar': 0.018069}),
    'ATLAS_PS_Rails': (1.481626, {'H': 0.017649, 'C': 0.120175, 'O': 0.080039, 'Ar': 0.782138}),
    'ATLAS_EMBF_CablesIn_4': (1.396017, {'Cu': 0.000012, 'C': 0.000005, 'H': 0.000000, 'O': 0.000002, 'N': 0.000001, 'Ar': 0.999980}),
    'ATLAS_EMBF_CablesIn_5': (1.498833, {'Cu': 0.078763, 'C': 0.033364, 'H': 0.001273, 'O': 0.010100, 'N': 0.003537, 'Ar': 0.872963}),
    'ATLAS_EMBF_CablesIn_6': (1.698294, {'Cu': 0.204345, 'C': 0.086560, 'H': 0.003302, 'O': 0.026205, 'N': 0.009176, 'Ar': 0.670412}),
    'ATLAS_EMBF_CablesIn_7': (1.897754, {'Cu': 0.303528, 'C': 0.128574, 'H': 0.004904, 'O': 0.038924, 'N': 0.013630, 'Ar': 0.510438}),
    'ATLAS_EMBF_CablesIn_8': (2.097215, {'Cu': 0.383846, 'C': 0.162597, 'H': 0.006202, 'O': 0.049224, 'N': 0.017237, 'Ar': 0.380894}),
    'ATLAS_EMBF_Boards_1': (2.028758, {'H': 0.017965, 'O': 0.089698, 'C': 0.169264, 'Cu': 0.313671, 'N': 0.006722, 'Ar': 0.402680}),
    'ATLAS_EMBF_Boards_2': (2.183688, {'H': 0.018006, 'O': 0.093779, 'C': 0.191758, 'Cu': 0.372867, 'N': 0.009903, 'Ar': 0.313687}),
    'ATLAS_EMBF_Boards_3': (2.338618, {'H': 0.018042, 'O': 0.097320, 'C': 0.211271, 'Cu': 0.424220, 'N': 0.012662, 'Ar': 0.236485}),
    'ATLAS_EMBF_Boards_4': (2.493536, {'H': 0.018074, 'O': 0.100420, 'C': 0.228358, 'Cu': 0.469188, 'N': 0.015078, 'Ar': 0.168882}),
    'ATLAS_EMBF_Boards_5': (2.568604, {'H': 0.018088, 'O': 0.101788, 'C': 0.235896, 'Cu': 0.489027, 'N': 0.016144, 'Ar': 0.139057}),
    'ATLAS_EMBF_Boards_6': (2.568604, {'H': 0.018088, 'O': 0.101788, 'C': 0.235896, 'Cu': 0.489027, 'N': 0.016144, 'Ar': 0.139057}),
    'ATLAS_EMBF_Boards_7': (2.568604, {'H': 0.018088, 'O': 0.101788, 'C': 0.235896, 'Cu': 0.489027, 'N': 0.016144, 'Ar': 0.139057}),
    'ATLAS_EMBF_Boards_8': (2.568604, {'H': 0.018088, 'O': 0.101788, 'C': 0.235896, 'Cu': 0.489027, 'N': 0.016144, 'Ar': 0.139057}),
    'ATLAS_EMBF_CablesOut_4': (1.396017, {'Cu': 0.000012, 'C': 0.000005, 'H': 0.000000, 'O': 0.000002, 'N': 0.000001, 'Ar': 0.999980}),
    'ATLAS_EMBF_CablesOut_5': (1.498833, {'Cu': 0.078763, 'C': 0.033364, 'H': 0.001273, 'O': 0.010100, 'N': 0.003537, 'Ar': 0.872963}),
    'ATLAS_EMBF_CablesOut_6': (1.698294, {'Cu': 0.204345, 'C': 0.086560, 'H': 0.003302, 'O': 0.026205, 'N': 0.009176, 'Ar': 0.670412}),
    'ATLAS_EMBF_CablesOut_7': (1.897754, {'Cu': 0.303528, 'C': 0.128574, 'H': 0.004904, 'O': 0.038924, 'N': 0.013630, 'Ar': 0.510438}),
    'ATLAS_EMBF_CablesOut_8': (2.097215, {'Cu': 0.383846, 'C': 0.162597, 'H': 0.006202, 'O': 0.049224, 'N': 0.017237, 'Ar': 0.380894}),
    'ATLAS_EMBF_SummingBoards': (1.607794, {'Ar': 0.843960, 'Cu': 0.156040}),
    'ATLAS_EMBF_G10Ring': (1.720003, {'H': 0.030783, 'O': 0.469777, 'C': 0.209611, 'Si': 0.289829}),
    'ATLAS_EMBF_AbsorberTips': (1.607935, {'H': 0.007184, 'O': 0.109636, 'C': 0.048919, 'Si': 0.067640, 'Fe': 0.106784, 'Ar': 0.659838}),
    'ATLAS_LAR_CABLES': (3.035179, {'Cu': 0.620000, 'C': 0.262632, 'H': 0.010018, 'O': 0.079508, 'N': 0.027842}),
}

# kind 'tube': (r0, r1, z0, z1) in mm, both_sides False means z0..z1 as given; kind 'pcon': planes (z, rmin, rmax)
ATLAS_BARREL_FRONT_VOLUMES = [
    dict(kind='tube', name='DM::Barrel::WarmInnerWall::Central', material='G4_Al', r0=1150.000, r1=1160.000, z0=-2900.000, z1=2900.000, in_field=True),  # InnerWall planes 11-12 (r 1150-1160 mm, |z| < 2956 mm); |z| < 2900 mm inside the field volume
    dict(kind='pcon', name='DM::Barrel::WarmInnerWall::EndA', material='G4_Al', planes=[(2900.000, 1150.000, 1160.000), (2956.000, 1150.000, 1160.000), (3040.000, 1150.000, 1170.000), (3040.000, 1150.000, 1208.000), (3065.000, 1150.000, 1208.000), (3065.000, 1178.000, 1208.000), (3125.000, 1178.000, 1208.000), (3125.000, 1150.000, 1208.000), (3150.000, 1150.000, 1208.000), (3175.000, 1150.000, 1228.000), (3175.000, 1215.000, 1228.000), (3385.000, 1385.000, 1401.000), (3405.000, 1400.990, 1401.000)], in_field=False),  # InnerWall planes 0-11, |z| from 2900 mm
    dict(kind='pcon', name='DM::Barrel::WarmInnerWall::EndB', material='G4_Al', planes=[(-3405.000, 1400.990, 1401.000), (-3385.000, 1385.000, 1401.000), (-3175.000, 1215.000, 1228.000), (-3175.000, 1150.000, 1228.000), (-3150.000, 1150.000, 1208.000), (-3125.000, 1150.000, 1208.000), (-3125.000, 1178.000, 1208.000), (-3065.000, 1178.000, 1208.000), (-3065.000, 1150.000, 1208.000), (-3040.000, 1150.000, 1208.000), (-3040.000, 1150.000, 1170.000), (-2956.000, 1150.000, 1160.000), (-2900.000, 1150.000, 1160.000)], in_field=False),  # InnerWall planes 12-23, |z| from 2900 mm
    dict(kind='tube', name='DM::Barrel::Solenoid::Al1', material='G4_Al', r0=1229.000, r1=1241.800, z0=-2700.000, z1=2700.000, in_field=False),  # CryoCylinders-09 cylinders 56-60
    dict(kind='tube', name='DM::Barrel::Solenoid::Cu', material='G4_Cu', r0=1241.800, r1=1244.850, z0=-2700.000, z1=2700.000, in_field=False),  # CryoCylinders-09 cylinders 56-60
    dict(kind='tube', name='DM::Barrel::Solenoid::Al2', material='G4_Al', r0=1244.850, r1=1258.650, z0=-2700.000, z1=2700.000, in_field=False),  # CryoCylinders-09 cylinders 56-60
    dict(kind='tube', name='DM::Barrel::Solenoid::G10', material='ATLAS_STD_G10', r0=1258.650, r1=1263.350, z0=-2700.000, z1=2700.000, in_field=False),  # CryoCylinders-09 cylinders 56-60
    dict(kind='tube', name='DM::Barrel::Solenoid::Al3', material='G4_Al', r0=1263.350, r1=1275.350, z0=-2840.000, z1=2840.000, in_field=False),  # CryoCylinders-09 cylinders 56-60
    dict(kind='pcon', name='DM::Barrel::ColdWall', material='G4_Al', planes=[(-3267.000, 1537.220, 1565.500), (-3101.000, 1356.720, 1385.000), (-3101.000, 1365.000, 1385.000), (-3023.000, 1365.000, 1385.000), (-2950.000, 1340.500, 1385.000), (-2728.000, 1340.500, 1385.000), (-2640.000, 1371.000, 1385.000), (-2400.000, 1369.000, 1385.000), (-1950.000, 1367.000, 1385.000), (-1120.000, 1358.900, 1385.000), (-400.000, 1356.600, 1385.000), (-75.000, 1341.000, 1385.000), (75.000, 1341.000, 1385.000), (400.000, 1356.600, 1385.000), (1120.000, 1358.900, 1385.000), (1950.000, 1367.000, 1385.000), (2400.000, 1369.000, 1385.000), (2640.000, 1371.000, 1385.000), (2728.000, 1340.500, 1385.000), (2950.000, 1340.500, 1385.000), (3023.000, 1365.000, 1385.000), (3101.000, 1365.000, 1385.000), (3101.000, 1356.720, 1385.000), (3267.000, 1537.220, 1565.500)], in_field=False),  # OuterWall planes 4-27 (CryoPcons-18)
    dict(kind='tube', name='DM::Barrel::LArBath', material='liquidArgon', r0=1385.000, r1=1410.000, z0=-3101.000, z1=3101.000, in_field=False),  # TotalLAr (CryoPcons-18) minus the PS mother
    dict(kind='tube', name='DM::Barrel::LArCentral', material='liquidArgon', r0=1410.000, r1=1500.000, z0=-3.000, z1=3.000, in_field=False),  # HalfLAr::Pos starts at z = 3 mm
    dict(kind='tube', name='DM::Barrel::PSMother::Shell::A', material='ATLAS_PS_Shell', r0=1411.400, r1=1411.800, z0=3.000, z1=3101.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Shell::B', material='ATLAS_PS_Shell', r0=1411.400, r1=1411.800, z0=-3101.000, z1=-3.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Prepreg1::A', material='ATLAS_PS_Prepreg1', r0=1412.300, r1=1413.300, z0=3.000, z1=3101.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Prepreg1::B', material='ATLAS_PS_Prepreg1', r0=1412.300, r1=1413.300, z0=-3101.000, z1=-3.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Prepreg2::A', material='ATLAS_PS_Prepreg2', r0=1426.300, r1=1430.800, z0=3.000, z1=3101.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Prepreg2::B', material='ATLAS_PS_Prepreg2', r0=1426.300, r1=1430.800, z0=-3101.000, z1=-3.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Board::A', material='ATLAS_PS_Board', r0=1430.800, r1=1433.000, z0=3.000, z1=3101.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Board::B', material='ATLAS_PS_Board', r0=1430.800, r1=1433.000, z0=-3101.000, z1=-3.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics1::A', material='ATLAS_PS_Connectics_1', r0=1433.000, r1=1438.000, z0=3.000, z1=390.250, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics1::B', material='ATLAS_PS_Connectics_1', r0=1433.000, r1=1438.000, z0=-390.250, z1=-3.000, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics2::A', material='ATLAS_PS_Connectics_2', r0=1433.000, r1=1438.000, z0=390.250, z1=777.500, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics2::B', material='ATLAS_PS_Connectics_2', r0=1433.000, r1=1438.000, z0=-777.500, z1=-390.250, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics3::A', material='ATLAS_PS_Connectics_3', r0=1433.000, r1=1438.000, z0=777.500, z1=1164.750, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics3::B', material='ATLAS_PS_Connectics_3', r0=1433.000, r1=1438.000, z0=-1164.750, z1=-777.500, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics4::A', material='ATLAS_PS_Connectics_4', r0=1433.000, r1=1438.000, z0=1164.750, z1=1552.000, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics4::B', material='ATLAS_PS_Connectics_4', r0=1433.000, r1=1438.000, z0=-1552.000, z1=-1164.750, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics5::A', material='ATLAS_PS_Connectics_5', r0=1433.000, r1=1438.000, z0=1552.000, z1=1939.250, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics5::B', material='ATLAS_PS_Connectics_5', r0=1433.000, r1=1438.000, z0=-1939.250, z1=-1552.000, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics6::A', material='ATLAS_PS_Connectics_6', r0=1433.000, r1=1438.000, z0=1939.250, z1=2326.500, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics6::B', material='ATLAS_PS_Connectics_6', r0=1433.000, r1=1438.000, z0=-2326.500, z1=-1939.250, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics7::A', material='ATLAS_PS_Connectics_7', r0=1433.000, r1=1438.000, z0=2326.500, z1=2713.750, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics7::B', material='ATLAS_PS_Connectics_7', r0=1433.000, r1=1438.000, z0=-2713.750, z1=-2326.500, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics8::A', material='ATLAS_PS_Connectics_8', r0=1433.000, r1=1438.000, z0=2713.750, z1=3101.000, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Connectics8::B', material='ATLAS_PS_Connectics_8', r0=1433.000, r1=1438.000, z0=-3101.000, z1=-2713.750, in_field=False),  # BarrelPresamplerConstruction.cxx, connectics growing with |z|
    dict(kind='tube', name='DM::Barrel::PSMother::Plate::A', material='ATLAS_PS_Plate', r0=1438.000, r1=1438.500, z0=3.000, z1=3101.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Plate::B', material='ATLAS_PS_Plate', r0=1438.000, r1=1438.500, z0=-3101.000, z1=-3.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Rails::A', material='ATLAS_PS_Rails', r0=1438.500, r1=1439.400, z0=3.000, z1=3101.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::Rails::B', material='ATLAS_PS_Rails', r0=1438.500, r1=1439.400, z0=-3101.000, z1=-3.000, in_field=False),  # BarrelPresamplerConstruction.cxx
    dict(kind='tube', name='DM::Barrel::PSMother::LAr1::A', material='liquidArgon', r0=1410.000, r1=1411.400, z0=3.000, z1=3101.000, in_field=False),  # PS mother, argon
    dict(kind='tube', name='DM::Barrel::PSMother::LAr1::B', material='liquidArgon', r0=1410.000, r1=1411.400, z0=-3101.000, z1=-3.000, in_field=False),  # PS mother, argon
    dict(kind='tube', name='DM::Barrel::PSMother::LAr2::A', material='liquidArgon', r0=1411.800, r1=1412.300, z0=3.000, z1=3101.000, in_field=False),  # PS mother, argon
    dict(kind='tube', name='DM::Barrel::PSMother::LAr2::B', material='liquidArgon', r0=1411.800, r1=1412.300, z0=-3101.000, z1=-3.000, in_field=False),  # PS mother, argon
    dict(kind='tube', name='DM::Barrel::PSMother::LAr3::A', material='liquidArgon', r0=1439.400, r1=1447.000, z0=3.000, z1=3101.000, in_field=False),  # PS mother, argon
    dict(kind='tube', name='DM::Barrel::PSMother::LAr3::B', material='liquidArgon', r0=1439.400, r1=1447.000, z0=-3101.000, z1=-3.000, in_field=False),  # PS mother, argon
    dict(kind='tube', name='DM::Barrel::EMBFront::LAr1::A', material='liquidArgon', r0=1447.000, r1=1449.180, z0=3.000, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::LAr1::B', material='liquidArgon', r0=1447.000, r1=1449.180, z0=-3097.000, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn1::A', material='liquidArgon', r0=1449.180, r1=1450.850, z0=3.000, z1=389.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn1::B', material='liquidArgon', r0=1449.180, r1=1450.850, z0=-389.750, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn2::A', material='liquidArgon', r0=1449.180, r1=1450.850, z0=389.750, z1=776.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn2::B', material='liquidArgon', r0=1449.180, r1=1450.850, z0=-776.500, z1=-389.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn3::A', material='liquidArgon', r0=1449.180, r1=1450.850, z0=776.500, z1=1163.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn3::B', material='liquidArgon', r0=1449.180, r1=1450.850, z0=-1163.250, z1=-776.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn4::A', material='ATLAS_EMBF_CablesIn_4', r0=1449.180, r1=1450.850, z0=1163.250, z1=1550.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn4::B', material='ATLAS_EMBF_CablesIn_4', r0=1449.180, r1=1450.850, z0=-1550.000, z1=-1163.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn5::A', material='ATLAS_EMBF_CablesIn_5', r0=1449.180, r1=1450.850, z0=1550.000, z1=1936.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn5::B', material='ATLAS_EMBF_CablesIn_5', r0=1449.180, r1=1450.850, z0=-1936.750, z1=-1550.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn6::A', material='ATLAS_EMBF_CablesIn_6', r0=1449.180, r1=1450.850, z0=1936.750, z1=2323.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn6::B', material='ATLAS_EMBF_CablesIn_6', r0=1449.180, r1=1450.850, z0=-2323.500, z1=-1936.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn7::A', material='ATLAS_EMBF_CablesIn_7', r0=1449.180, r1=1450.850, z0=2323.500, z1=2710.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn7::B', material='ATLAS_EMBF_CablesIn_7', r0=1449.180, r1=1450.850, z0=-2710.250, z1=-2323.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn8::A', material='ATLAS_EMBF_CablesIn_8', r0=1449.180, r1=1450.850, z0=2710.250, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesIn8::B', material='ATLAS_EMBF_CablesIn_8', r0=1449.180, r1=1450.850, z0=-3097.000, z1=-2710.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards1::A', material='ATLAS_EMBF_Boards_1', r0=1450.850, r1=1455.150, z0=3.000, z1=389.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards1::B', material='ATLAS_EMBF_Boards_1', r0=1450.850, r1=1455.150, z0=-389.750, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards2::A', material='ATLAS_EMBF_Boards_2', r0=1450.850, r1=1455.150, z0=389.750, z1=776.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards2::B', material='ATLAS_EMBF_Boards_2', r0=1450.850, r1=1455.150, z0=-776.500, z1=-389.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards3::A', material='ATLAS_EMBF_Boards_3', r0=1450.850, r1=1455.150, z0=776.500, z1=1163.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards3::B', material='ATLAS_EMBF_Boards_3', r0=1450.850, r1=1455.150, z0=-1163.250, z1=-776.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards4::A', material='ATLAS_EMBF_Boards_4', r0=1450.850, r1=1455.150, z0=1163.250, z1=1550.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards4::B', material='ATLAS_EMBF_Boards_4', r0=1450.850, r1=1455.150, z0=-1550.000, z1=-1163.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards5::A', material='ATLAS_EMBF_Boards_5', r0=1450.850, r1=1455.150, z0=1550.000, z1=1936.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards5::B', material='ATLAS_EMBF_Boards_5', r0=1450.850, r1=1455.150, z0=-1936.750, z1=-1550.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards6::A', material='ATLAS_EMBF_Boards_6', r0=1450.850, r1=1455.150, z0=1936.750, z1=2323.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards6::B', material='ATLAS_EMBF_Boards_6', r0=1450.850, r1=1455.150, z0=-2323.500, z1=-1936.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards7::A', material='ATLAS_EMBF_Boards_7', r0=1450.850, r1=1455.150, z0=2323.500, z1=2710.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards7::B', material='ATLAS_EMBF_Boards_7', r0=1450.850, r1=1455.150, z0=-2710.250, z1=-2323.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards8::A', material='ATLAS_EMBF_Boards_8', r0=1450.850, r1=1455.150, z0=2710.250, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::Boards8::B', material='ATLAS_EMBF_Boards_8', r0=1450.850, r1=1455.150, z0=-3097.000, z1=-2710.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut1::A', material='liquidArgon', r0=1455.150, r1=1456.820, z0=3.000, z1=389.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut1::B', material='liquidArgon', r0=1455.150, r1=1456.820, z0=-389.750, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut2::A', material='liquidArgon', r0=1455.150, r1=1456.820, z0=389.750, z1=776.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut2::B', material='liquidArgon', r0=1455.150, r1=1456.820, z0=-776.500, z1=-389.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut3::A', material='liquidArgon', r0=1455.150, r1=1456.820, z0=776.500, z1=1163.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut3::B', material='liquidArgon', r0=1455.150, r1=1456.820, z0=-1163.250, z1=-776.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut4::A', material='ATLAS_EMBF_CablesOut_4', r0=1455.150, r1=1456.820, z0=1163.250, z1=1550.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut4::B', material='ATLAS_EMBF_CablesOut_4', r0=1455.150, r1=1456.820, z0=-1550.000, z1=-1163.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut5::A', material='ATLAS_EMBF_CablesOut_5', r0=1455.150, r1=1456.820, z0=1550.000, z1=1936.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut5::B', material='ATLAS_EMBF_CablesOut_5', r0=1455.150, r1=1456.820, z0=-1936.750, z1=-1550.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut6::A', material='ATLAS_EMBF_CablesOut_6', r0=1455.150, r1=1456.820, z0=1936.750, z1=2323.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut6::B', material='ATLAS_EMBF_CablesOut_6', r0=1455.150, r1=1456.820, z0=-2323.500, z1=-1936.750, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut7::A', material='ATLAS_EMBF_CablesOut_7', r0=1455.150, r1=1456.820, z0=2323.500, z1=2710.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut7::B', material='ATLAS_EMBF_CablesOut_7', r0=1455.150, r1=1456.820, z0=-2710.250, z1=-2323.500, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut8::A', material='ATLAS_EMBF_CablesOut_8', r0=1455.150, r1=1456.820, z0=2710.250, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::CablesOut8::B', material='ATLAS_EMBF_CablesOut_8', r0=1455.150, r1=1456.820, z0=-3097.000, z1=-2710.250, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::LAr2::A', material='liquidArgon', r0=1456.820, r1=1460.000, z0=3.000, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::LAr2::B', material='liquidArgon', r0=1456.820, r1=1460.000, z0=-3097.000, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::SummingBoards::A', material='ATLAS_EMBF_SummingBoards', r0=1460.000, r1=1470.000, z0=3.000, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::SummingBoards::B', material='ATLAS_EMBF_SummingBoards', r0=1460.000, r1=1470.000, z0=-3097.000, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::G10Ring::A', material='ATLAS_EMBF_G10Ring', r0=1470.000, r1=1490.000, z0=3.000, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::G10Ring::B', material='ATLAS_EMBF_G10Ring', r0=1470.000, r1=1490.000, z0=-3097.000, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::AbsorberTips::A', material='ATLAS_EMBF_AbsorberTips', r0=1490.000, r1=1498.000, z0=3.000, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::AbsorberTips::B', material='ATLAS_EMBF_AbsorberTips', r0=1490.000, r1=1498.000, z0=-3097.000, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::LAr3::A', material='liquidArgon', r0=1498.000, r1=1500.000, z0=3.000, z1=3097.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='tube', name='DM::Barrel::EMBFront::LAr3::B', material='liquidArgon', r0=1498.000, r1=1500.000, z0=-3097.000, z1=-3.000, in_field=False),  # BarrelConstruction.cxx (front electronics, G10 ring, absorber tips)
    dict(kind='pcon', name='DM::Barrel::EndLAr::A', material='liquidArgon', planes=[(3097.000, 1447.000, 1500.000), (3101.000, 1447.000, 1500.000), (3101.000, 1385.001, 1500.000), (3150.000, 1438.281, 1500.000), (3205.000, 1498.085, 1500.000), (3205.000, 1498.085, 1565.500), (3266.999, 1565.500, 1565.500)], in_field=False),  # TotalLAr between the cold cone and the accordion
    dict(kind='pcon', name='DM::Barrel::EndLAr::B', material='liquidArgon', planes=[(-3266.999, 1565.500, 1565.500), (-3205.000, 1498.085, 1565.500), (-3205.000, 1498.085, 1500.000), (-3150.000, 1438.281, 1500.000), (-3101.000, 1385.001, 1500.000), (-3101.000, 1447.000, 1500.000), (-3097.000, 1447.000, 1500.000)], in_field=False),  # TotalLAr between the cold cone and the accordion
]

# solenoid field volume radius with this front (the first aluminium cylinder of the coil starts at r = 1229 mm)
ATLAS_BARREL_FRONT_FIELD_RADIUS = 1229.0

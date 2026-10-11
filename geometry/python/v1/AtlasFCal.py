"""
The ATLAS forward calorimeter (FCal) as passive material, without cells (F41 conserto 20, generated on 11/10/2026 by
Simulacao_MinBias/41_consertos_lorenzetti/c20_fcal/c20_gerar_fcal.py; do not edit by hand).

Source: the public ATLAS geometry of the GeoModel project, geometry-ATLAS-R3S-2021-03-00-00.db
(https://geomodel.web.cern.ch/atlas-geometry-data/old/): the three modules (LAr::FCAL::Module1-3), their absorbers,
electrode gaps and rods, and their position in the end-cap cryostat, without the shift of the installed end-cap
(45 mm). Each module is a tube of the mean mixture of absorber, rods and liquid argon; the liquid argon between the
modules is a tube of liquid argon. Not modelled: the cable troughs at the outer radius and the support structure.
"""

# name: (density in g/cm3, {element symbol: mass fraction})
ATLAS_FCAL_MIXTURES = {
    'ATLAS_FCAL1_MODULE': (8.096758, {'Cu': 0.979824, 'Ar': 0.020176}),
    'ATLAS_FCAL2_MODULE': (14.376014, {'Cu': 0.087091, 'Ar': 0.015865, 'W': 0.870133, 'Fe': 0.008970, 'Ni': 0.017941}),
    'ATLAS_FCAL3_MODULE': (14.088338, {'Cu': 0.084872, 'Ar': 0.018212, 'W': 0.870009, 'Fe': 0.008969, 'Ni': 0.017938}),
}

# tubes (mm); one entry per side
ATLAS_FCAL_VOLUMES = [
    {'kind': 'tube', 'name': 'DM::FCal::FCal1::A', 'material': 'ATLAS_FCAL1_MODULE', 'r0': 72.0, 'r1': 449.6, 'z0': 4668.5, 'z1': 5112.700000000001, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::Module, mean mixture'},
    {'kind': 'tube', 'name': 'DM::FCal::FCal1::B', 'material': 'ATLAS_FCAL1_MODULE', 'r0': 72.0, 'r1': 449.6, 'z0': -5112.700000000001, 'z1': -4668.5, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::Module, mean mixture'},
    {'kind': 'tube', 'name': 'DM::FCal::FCal2::A', 'material': 'ATLAS_FCAL2_MODULE', 'r0': 79.0, 'r1': 449.6, 'z0': 5128.3, 'z1': 5570.599999999999, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::Module, mean mixture'},
    {'kind': 'tube', 'name': 'DM::FCal::FCal2::B', 'material': 'ATLAS_FCAL2_MODULE', 'r0': 79.0, 'r1': 449.6, 'z0': -5570.599999999999, 'z1': -5128.3, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::Module, mean mixture'},
    {'kind': 'tube', 'name': 'DM::FCal::FCal3::A', 'material': 'ATLAS_FCAL3_MODULE', 'r0': 86.0, 'r1': 449.6, 'z0': 5602.8, 'z1': 6045.099999999999, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::Module, mean mixture'},
    {'kind': 'tube', 'name': 'DM::FCal::FCal3::B', 'material': 'ATLAS_FCAL3_MODULE', 'r0': 86.0, 'r1': 449.6, 'z0': -6045.099999999999, 'z1': -5602.8, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::Module, mean mixture'},
    {'kind': 'tube', 'name': 'DM::FCal::LAr12::A', 'material': 'liquidArgon', 'r0': 79.0, 'r1': 449.6, 'z0': 5112.700000000001, 'z1': 5128.3, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::LiquidArgonC'},
    {'kind': 'tube', 'name': 'DM::FCal::LAr12::B', 'material': 'liquidArgon', 'r0': 79.0, 'r1': 449.6, 'z0': -5128.3, 'z1': -5112.700000000001, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::LiquidArgonC'},
    {'kind': 'tube', 'name': 'DM::FCal::LAr23::A', 'material': 'liquidArgon', 'r0': 86.0, 'r1': 449.6, 'z0': 5570.599999999999, 'z1': 5602.8, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::LiquidArgonC'},
    {'kind': 'tube', 'name': 'DM::FCal::LAr23::B', 'material': 'liquidArgon', 'r0': 86.0, 'r1': 449.6, 'z0': -5602.8, 'z1': -5570.599999999999, 'both': False, 'in_field': False, 'source': 'GeoModel LAr::FCAL::LiquidArgonC'},
]

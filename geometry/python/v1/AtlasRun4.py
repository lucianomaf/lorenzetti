"""
The ATLAS Run 4 geometry of the simulation (F41 conserto 27, --atlas-run4): the geometry of the ATLAS simulation of the
reference production (ATLAS-P2-RUN4-01-01-00, Athena 23.0.14), read from its GeoModel SQLite dump made in the simulation mode
(DumpGeo with Common.ProductionStep=Simulation; geometry/scripts/atlas_run4_dumpgeo.sh).

LAr end-cap position: the ATLAS simulation places the end-cap cryostat at |z| = 6734.5 mm (LArPosition-10, LARCRYO_EC_POS), with the
alignment off. Its presampler, EMEC, HEC and FCal sit 43 to 45 mm further from the interaction point than the nominal numbers of the
readout geometry used by the ATLAS-like end-cap options (PresamplerPosition ZPOS, LArWheelCalculator, HecNominals, FCal without the
shift). shiftLArEndcap moves each of them by its own ATLAS_ENDCAP_SIM_SHIFT. The cells move with them, as in ATLAS, where the end-cap identifiers come from
local coordinates (LArG4EC PresamplerGeometry.cc: local z plus the nominal front face; LArWheelCalculator: the wheel frame; HEC and
FCal: the module frame).
"""

__all__ = ["ATLAS_ENDCAP_SIM_SHIFT", "shiftLArEndcap"]

import math
from GaugiKernel.constants import mm

# |z| of each part of the end-cap of the ATLAS simulation minus its nominal |z| in the ATLAS-like options, measured in the dump of
# ATLAS-P2-RUN4-01-01-00 (simulation mode, alignment off): presampler argon 3667-3671 against 3622-3626 mm; front face of the EMEC wheels
# 3745 against 3702 mm (LAr::EMEC::Pos::OuterWheel, 514 mm thick); front of the liquid argon of the HEC wheels 4321 and 5178 against
# 4277 and 5134 mm (LAr::HEC::LiquidArgon, HecNominals); FCal1 4713.5-5157.7 against 4668.5-5112.7 mm. By prefix of the volume names
# built by the ATLAS-like end-cap options (ECAL.py getAtlasPseCfg, EMEC.py, HEC.py, AtlasFCal.py, FCAL.py).
ATLAS_ENDCAP_SIM_SHIFT = {
    'LAr::PSE' : 45.0*mm,
    'LAr::EMEC': 43.0*mm,
    'FCal::HEC': 44.0*mm,
    'DM::HEC'  : 44.0*mm,
    'DM::FCal' : 45.0*mm,
    'LAr::FCAL': 45.0*mm,
}


def _shiftOf( name ):
    for prefix, shift in ATLAS_ENDCAP_SIM_SHIFT.items():
        if name.startswith( prefix ): return shift
    return None


def _isEndcap( name ):
    return _shiftOf( name ) is not None


def _shiftVolume( pv, dz ):
    pv.Z += dz; pv.ZMin += dz; pv.ZMax += dz
    if getattr( pv, 'ZPlanes', None ):
        pv.ZPlanes = [ z + dz for z in pv.ZPlanes ]


def shiftLArEndcap( volumes, samplings, fcal_modules ):
    """
    Moves the LAr end-cap volumes, hit makers and FCal modules away from the interaction point on each side (ATLAS_ENDCAP_SIM_SHIFT), with the fields of
    their cells that hold an absolute z: the electric focal point and the wheel boxes of the EMEC (the boundaries ZSEP and ZIW are measured
    from the focal point and move with it), the (r, z) boxes of the HEC cells, and the z subtracted from the end-cap presampler before its
    eta (EmbZShift). The eta of the centre of the FCal tiles (metadata of the hits) follows the module.
    """
    done = set()
    for pv in volumes:
        if _isEndcap( pv.Name ) and id(pv) not in done:
            _shiftVolume( pv, math.copysign( _shiftOf( pv.Name ), pv.Z ) ); done.add( id(pv) )
    for s in samplings:
        pv = s.volume()
        if not _isEndcap( pv.Name ) or id(pv) in done: continue
        sign = 1 if pv.Z > 0 else -1
        shift = _shiftOf( pv.Name )
        _shiftVolume( pv, sign * shift ); done.add( id(pv) )
        cells = getattr( s.sensitive(), 'Cells', None )
        if not cells: continue
        if 'EmecFocalShift' in cells:
            cells['EmecFocalShift'] = cells['EmecFocalShift'] + shift
            for k in ('EmecWheelZMin', 'EmecWheelZMax'):
                cells[k] = [ z + math.copysign( shift, z ) for z in cells[k] ]
        if 'BoxZMin' in cells:
            for k in ('BoxZMin', 'BoxZMax'):
                cells[k] = [ z + math.copysign( shift, z ) for z in cells[k] ]
        if cells.get('EmbMode') == 2:
            cells['EmbZShift'] = shift
    for fm in fcal_modules:
        if id(fm) in done: continue
        dz = math.copysign( _shiftOf( fm.Name ), fm.Z )
        rs = [ fm.Z / math.sinh(eta) for eta in fm.CopyEta ]
        fm.Z += dz; fm.ZMin += dz; fm.ZMax += dz
        fm.CopyEta = [ math.asinh( fm.Z / r ) for r in rs ]
        done.add( id(fm) )
